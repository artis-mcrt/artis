#include "kpkt.h"

#include <gsl/gsl_integration.h>

#include <cmath>

#include "atomic.h"
#include "grid.h"
#include "ltepop.h"
#include "macroatom.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "vectors.h"
#include "vpkt.h"

namespace kpkt {

enum coolingtype {
  COOLINGTYPE_FF = 880,
  COOLINGTYPE_FB = 881,
  COOLINGTYPE_COLLEXC = 882,
  COOLINGTYPE_COLLION = 883,
};

struct cellhistorycoolinglist {
  enum coolingtype type;
  int element;
  int ion;
  int level;
  int upperlevel;
};

static __managed__ struct cellhistorycoolinglist *coolinglist;

__host__ __device__ int get_coolinglistoffset(int element, int ion) {
  return globals::elements[element].ions[ion].coolingoffset;
}

__host__ __device__ static int get_ncoolingterms(int element, int ion) {
  return globals::elements[element].ions[ion].ncoolingterms;
}

__host__ __device__ static double get_cooling_ion_coll_exc(const int modelgridindex, const int element, const int ion,
                                                           const double T_e, const double nne) {
  double C_exc = 0.;
  const int nlevels = get_nlevels(element, ion);

  /// excitation to same ionization stage
  /// -----------------------------------
  for (int level = 0; level < nlevels; level++) {
    const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
    const double epsilon_current = epsilon(element, ion, level);
    const int nuptrans = get_nuptrans(element, ion, level);
    for (int ii = 0; ii < nuptrans; ii++) {
      const int lineindex = globals::elements[element].ions[ion].levels[level].uptrans_lineindicies[ii];
      const int upper = globals::linelist[lineindex].upperlevelindex;
      // printout("    excitation to level %d possible\n",upper);
      const double epsilon_trans = epsilon(element, ion, upper) - epsilon_current;
      const double C = nnlevel * col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans) * epsilon_trans;
      C_exc += C;
    }
  }

  return C_exc;
}

__host__ __device__ static double get_bfcoolingcoeff(int element, int ion, int level, int phixstargetindex, float T_e) {
  const int lowerindex = floor(log(T_e / MINTEMP) / T_step_log);
  if (lowerindex < TABLESIZE - 1) {
    const int upperindex = lowerindex + 1;
    const double T_lower = MINTEMP * exp(lowerindex * T_step_log);
    const double T_upper = MINTEMP * exp(upperindex * T_step_log);

    const double f_upper =
        globals::elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[upperindex];
    const double f_lower =
        globals::elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[lowerindex];

    return (f_lower + (f_upper - f_lower) / (T_upper - T_lower) * (T_e - T_lower));
  } else
    return globals::elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[TABLESIZE -
                                                                                                             1];
}

__host__ __device__ void calculate_cooling_rates(const int modelgridindex,
                                                 struct heatingcoolingrates *heatingcoolingrates)
// Calculate the cooling rates for a given cell and store them for each ion
// optionally store components (ff, bf, collisional) in heatingcoolingrates struct
{
  const float nne = grid::get_nne(modelgridindex);
  const float T_e = grid::get_Te(modelgridindex);

  double C_total = 0.;
  double C_ff_all = 0.;          /// free-free creation of rpkts
  double C_fb_all = 0.;          /// free-bound creation of rpkt
  double C_exc_all = 0.;         /// collisional excitation of macroatoms
  double C_ionization_all = 0.;  /// collisional ionisation of macroatoms
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      double C_ion = 0.;  /// all cooling for an ion
      const int nionisinglevels = get_ionisinglevels(element, ion);
      const double nncurrention = ionstagepop(modelgridindex, element, ion);

      /// ff creation of rpkt
      const int ioncharge = get_ionstage(element, ion) - 1;
      if (ioncharge > 0) {
        const double C_ff_ion = 1.426e-27 * sqrt(T_e) * pow(ioncharge, 2) * nncurrention * nne;
        C_ff_all += C_ff_ion;
        C_ion += C_ff_ion;
      }

      const double C_exc_ion = get_cooling_ion_coll_exc(modelgridindex, element, ion, T_e, nne);
      C_exc_all += C_exc_ion;
      C_ion += C_exc_ion;

      if (ion < nions - 1) {
        for (int level = 0; level < nionisinglevels; level++) {
          // printout("[debug] do_kpkt: element %d, ion %d, level %d\n",element,ion,level);
          const double epsilon_current = epsilon(element, ion, level);
          const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
          // printout("    ionisation possible\n");
          /// ionization to higher ionization stage
          /// -------------------------------------
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level);
               phixstargetindex++) {
            const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
            const double epsilon_trans = epsilon(element, ion + 1, upper) - epsilon_current;
            // printout("cooling list: col_ionization\n");
            const double C_ionization_ion_thistarget =
                nnlevel * col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans) *
                epsilon_trans;
            C_ionization_all += C_ionization_ion_thistarget;
            C_ion += C_ionization_ion_thistarget;
          }

          /// fb creation of r-pkt
          /// free bound rates are calculated from the lower ion, but associated to the higher ion
          /// --------------------
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level);
               phixstargetindex++) {
            // const int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            // const double nnupperlevel = get_levelpop(modelgridindex, element, ion + 1, upper);
            const double nnupperion = ionstagepop(modelgridindex, element, ion + 1);

            const double C_fb_ion_thistarget =
                get_bfcoolingcoeff(element, ion, level, phixstargetindex, T_e) * nnupperion * nne;
            C_fb_all += C_fb_ion_thistarget;
            C_ion += C_fb_ion_thistarget;
          }
        }
      }

      C_total += C_ion;
      grid::modelgrid[modelgridindex].cooling_contrib_ion[element][ion] = C_ion;
    }
  }
  grid::modelgrid[modelgridindex].totalcooling = C_total;

  // only used in the T_e solver and write_to_estimators file
  if (heatingcoolingrates != NULL) {
    heatingcoolingrates->cooling_collisional = C_exc_all + C_ionization_all;
    heatingcoolingrates->cooling_fb = C_fb_all;
    heatingcoolingrates->cooling_ff = C_ff_all;
  }
}

__host__ __device__ static void calculate_kpkt_rates_ion(int modelgridindex, int element, int ion, int indexionstart,
                                                         double oldcoolingsum, int tid)
/// Set up the global cooling list and determine the important entries
/// by applying a cut to the total cooling rate. Then sort the global
/// cooling list by the strength of the individual process contribution.
{
  const float nne = grid::get_nne(modelgridindex);
  const float T_e = grid::get_Te(modelgridindex);
  // double T_R = grid::get_TR(modelgridindex);
  // double W = grid::get_W(modelgridindex);

  /// calculate rates for
  // C_ff = 0.;   /// free-free creation of rpkts
  // C_fb = 0.;   /// free-bound creation of rpkt
  // C_exc = 0.;  /// collisional excitation of macroatoms
  // C_ion = 0.;  /// collisional ionisation of macroatoms
  double contrib = oldcoolingsum;
  int i = indexionstart;

  // printout("[debug] do_kpkt: element %d\n",element);
  const int nions = get_nions(element);
  const int nlevels_currention = get_nlevels(element, ion);

  const int ionisinglevels = get_ionisinglevels(element, ion);
  const double nncurrention = ionstagepop(modelgridindex, element, ion);

  /// ff creation of rpkt
  /// -------------------
  const int ioncharge = get_ionstage(element, ion) - 1;
  // printout("[debug] ioncharge %d, nncurrention %g, nne %g\n",ion,nncurrention,nne);
  if (ioncharge > 0) {
    const double C = 1.426e-27 * sqrt(T_e) * pow(ioncharge, 2) * nncurrention * nne;
    // C_ff += C;
    // C_ion += C;
    contrib += C;
    globals::cellhistory[tid].cooling_contrib[i] = contrib;

    assert_testmodeonly(coolinglist[i].type == COOLINGTYPE_FF);
    assert_testmodeonly(coolinglist[i].element == element);
    assert_testmodeonly(coolinglist[i].ion == ion);

    // if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d,
    // coolingtype %d, i %d, low
    // %d\n",contrib,oldcoolingsum,C,element,ion,-99,globals::cellhistory[tid].coolinglist[i].type,i,low);
    i++;
  }

  for (int level = 0; level < nlevels_currention; level++) {
    // printout("[debug] do_kpkt: element %d, ion %d, level %d\n",element,ion,level);
    const double epsilon_current = epsilon(element, ion, level);
    const double nnlevel = get_levelpop(modelgridindex, element, ion, level);

    const int nuptrans = get_nuptrans(element, ion, level);
    if (nuptrans > 0) {
      /// excitation to same ionization stage
      /// -----------------------------------
      for (int ii = 0; ii < nuptrans; ii++) {
        const int lineindex = globals::elements[element].ions[ion].levels[level].uptrans_lineindicies[ii];
        const int upper = globals::linelist[lineindex].upperlevelindex;
        // printout("    excitation to level %d possible\n",upper);
        const double epsilon_trans = epsilon(element, ion, upper) - epsilon_current;
        const double C = nnlevel * col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans) * epsilon_trans;
        contrib += C;
      }
      globals::cellhistory[tid].cooling_contrib[i] = contrib;

      assert_testmodeonly(coolinglist[i].type == COOLINGTYPE_COLLEXC);
      assert_testmodeonly(coolinglist[i].element == element);
      assert_testmodeonly(coolinglist[i].ion == ion);

      // if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d,
      // coolingtype %d, i %d, low
      // %d\n",contrib,oldcoolingsum,C,element,ion,level,globals::cellhistory[tid].coolinglist[i].type,i,low);
      i++;
    }

    if (ion < (nions - 1) && level < ionisinglevels)  /// check whether further ionisation stage available
    {
      // printout("    ionisation possible\n");
      /// ionization to higher ionization stage
      /// -------------------------------------
      /// for the moment we deal only with ionisations to the next ions groundlevel
      /// to allow an easy generalisation this was only made sure by col_ionization
      /// for speed issues (reduced number of calls to epsilon) it is now done also
      /// here explicitly
      const int nphixstargets = get_nphixstargets(element, ion, level);
      for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
        const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
        const double epsilon_upper = epsilon(element, ion + 1, upper);
        const double epsilon_trans = epsilon_upper - epsilon_current;
        // printout("cooling list: col_ionization\n");
        const double C = nnlevel *
                         col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans) *
                         epsilon_trans;

        contrib += C;
        globals::cellhistory[tid].cooling_contrib[i] = contrib;

        assert_testmodeonly(coolinglist[i].type == COOLINGTYPE_COLLION);
        assert_testmodeonly(coolinglist[i].element == element);
        assert_testmodeonly(coolinglist[i].ion == ion);
        // if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d,
        // coolingtype %d, i %d, low
        // %d\n",contrib,oldcoolingsum,C,element,ion,level,globals::cellhistory[tid].coolinglist[i].type,i,low);
        i++;
      }
      // C_ion += C;
      // C = 0.;
      // globals::cellhistory[tid].coolinglist[i].contribution = C;
      // }

      /// fb creation of r-pkt
      /// free bound rates are calculated from the lower ion, but associated to the higher ion
      /// --------------------
      for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++) {
        // const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
        // const double nnupperlevel = get_levelpop(modelgridindex,element,ion + 1, upper);
        const double nnupperion = ionstagepop(modelgridindex, element, ion + 1);

        const double C = get_bfcoolingcoeff(element, ion, level, phixstargetindex, T_e) * nnupperion * nne;
        contrib += C;
        globals::cellhistory[tid].cooling_contrib[i] = contrib;

        assert_testmodeonly(coolinglist[i].type == COOLINGTYPE_FB);
        assert_testmodeonly(coolinglist[i].element == element);
        assert_testmodeonly(coolinglist[i].ion == ion);

        i++;
      }
    }
  }

  assert_always(indexionstart == get_coolinglistoffset(element, ion));
  assert_always(i == indexionstart + get_ncoolingterms(element, ion));
  // we just summed up every individual cooling process. make sure it matches the stored total for the ion
  assert_always(fabs((grid::modelgrid[modelgridindex].cooling_contrib_ion[element][ion] + oldcoolingsum - contrib) /
                     contrib) < 1e-3);
}

__host__ __device__ static void set_ncoolingterms(void) {
  globals::ncoolingterms = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      int ionterms = 0;
      globals::elements[element].ions[ion].coolingoffset = globals::ncoolingterms;

      /// Ionised ions add one ff-cooling term
      if (get_ionstage(element, ion) > 1) ionterms++;
      /// Ionisinglevels below the closure ion add to bf and col ionisation
      /// All the levels add number of col excitations
      const int nlevels = get_nlevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        // if (ion < nions - 1) and (level < get_ionisinglevels(element,ion))
        if (ion < nions - 1) {
          ionterms += 2 * get_nphixstargets(element, ion, level);
        }

        if (get_nuptrans(element, ion, level) > 0) {
          ionterms++;  // level's coll. excitation cooling (all upper levels combined)
        }
      }
      globals::elements[element].ions[ion].ncoolingterms = ionterms;
      globals::ncoolingterms += ionterms;
    }
  }
}

void setup_coolinglist(void) {
  /// SET UP THE COOLING LIST
  ///======================================================
  /// Determine number of processes which allow kpkts to convert to something else.
  /// This number is given by the collisional excitations (so far determined from the oscillator strengths
  /// by the van Regemorter formula, therefore totaluptrans), the number of free-bound emissions and collisional
  /// ionisations (as long as we only deal with ionisation to the ground level this means for both of these
  /// \sum_{elements,ions}get_nlevels(element,ion) and free-free which is \sum_{elements} get_nions(element)-1

  set_ncoolingterms();
  const long mem_usage_coolinglist = globals::ncoolingterms * sizeof(struct cellhistorycoolinglist);
  coolinglist = (struct cellhistorycoolinglist *)malloc(globals::ncoolingterms * sizeof(struct cellhistorycoolinglist));
  printout("[info] mem_usage: coolinglist occupies %.3f MB\n", mem_usage_coolinglist / 1024. / 1024.);

  int i = 0;  // cooling list index
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      /// calculate rates for
      // C_ff = 0.;   /// free-free creation of rpkts
      // C_fb = 0.;   /// free-bound creation of rpkt
      // C_exc = 0.;  /// collisional excitation of macroatoms
      // C_ion = 0.;  /// collisional ionisation of macroatoms

      // printout("[debug] do_kpkt: element %d\n",element);
      const int nlevels_currention = get_nlevels(element, ion);

      const int ionisinglevels = get_ionisinglevels(element, ion);

      /// ff creation of rpkt
      /// -------------------
      const int ioncharge = get_ionstage(element, ion) - 1;
      // printout("[debug] ioncharge %d, nncurrention %g, nne %g\n",ion,nncurrention,nne);
      if (ioncharge > 0) {
        coolinglist[i].type = COOLINGTYPE_FF;
        coolinglist[i].element = element;
        coolinglist[i].ion = ion;
        coolinglist[i].level = -99;
        coolinglist[i].upperlevel = -99;
        i++;
      }

      for (int level = 0; level < nlevels_currention; level++) {
        if (get_nuptrans(element, ion, level) > 0) {
          coolinglist[i].type = COOLINGTYPE_COLLEXC;
          coolinglist[i].element = element;
          coolinglist[i].ion = ion;
          coolinglist[i].level = level;
          // upper level is not valid because this is the contribution of all upper levels combined - have to calculate
          // individually when selecting a random process
          coolinglist[i].upperlevel = -1;
          i++;
        }

        if (ion < (nions - 1) && level < ionisinglevels)  /// check whether further ionisation stage available
        {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
            coolinglist[i].type = COOLINGTYPE_COLLION;
            coolinglist[i].element = element;
            coolinglist[i].ion = ion;
            coolinglist[i].level = level;
            coolinglist[i].upperlevel = upper;
            i++;
          }

          /// fb creation of r-pkt
          /// free bound rates are calculated from the lower ion, but associated to the higher ion
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level);
               phixstargetindex++) {
            const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
            coolinglist[i].type = COOLINGTYPE_FB;
            coolinglist[i].element = element;
            coolinglist[i].ion = ion;
            coolinglist[i].level = level;
            coolinglist[i].upperlevel = upper;
            i++;
          }
        }
      }
      assert_always(i == get_coolinglistoffset(element, ion) + get_ncoolingterms(element, ion));
    }
  }
  assert_always(globals::ncoolingterms == i);  // if this doesn't match, we miscalculated the number of cooling terms
  printout("[info] read_atomicdata: number of coolingterms %d\n", globals::ncoolingterms);
}

__host__ __device__ static double sample_planck(const double T)
/// returns a randomly chosen frequency according to the Planck
/// distribution of temperature T
{
  const double nu_peak = 5.879e10 * T;
  if (nu_peak > globals::nu_max_r || nu_peak < globals::nu_min_r)
    printout("[warning] sample_planck: intensity peaks outside frequency range\n");

  const double B_peak = radfield::dbb(nu_peak, T, 1);

  double nu;
  bool endloop = false;
  // int i = 0;
  while (!endloop) {
    // i++;
    double zrand = gsl_rng_uniform(rng);
    double zrand2 = gsl_rng_uniform(rng);
    nu = globals::nu_min_r + zrand * (globals::nu_max_r - globals::nu_min_r);
    if (zrand2 * B_peak <= radfield::dbb(nu, T, 1)) endloop = true;
    // printout("[debug] sample_planck: planck_sampling %d\n", i);
  }

  return nu;
}

__host__ __device__ double do_kpkt_bb(struct packet *pkt_ptr)
/// Now routine to deal with a k-packet. Similar idea to do_gamma.
{
  // double nne = globals::cell[pkt_ptr->where].nne ;
  int cellindex = pkt_ptr->where;
  const int modelgridindex = grid::get_cell_modelgridindex(cellindex);
  const float T_e = grid::get_Te(modelgridindex);

  pkt_ptr->nu_cmf = sample_planck(T_e);
  if (!std::isfinite(pkt_ptr->nu_cmf)) {
    printout("[fatal] do_kpkt_bb: selected frequency not finite ... abort\n");
    abort();
  }
  /// and then emitt the packet randomly in the comoving frame
  emitt_rpkt(pkt_ptr);
  // printout("[debug] calculate_kappa_rpkt after kpkt to rpkt by ff\n");
  cellindex = pkt_ptr->where;
  pkt_ptr->next_trans = 0;  /// FLAG: transition history here not important, cont. process
  // if (tid == 0) k_stat_to_r_bb++;
  stats::increment(stats::COUNTER_K_STAT_TO_R_BB);
  pkt_ptr->interactions++;
  pkt_ptr->last_event = 6;
  pkt_ptr->emissiontype = -9999999;
  vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
  pkt_ptr->em_time = pkt_ptr->prop_time;
  pkt_ptr->nscatterings = 0;

  return pkt_ptr->prop_time;
}

__host__ __device__ double do_kpkt(struct packet *pkt_ptr, double t2, int nts)
/// Now routine to deal with a k-packet. Similar idea to do_gamma.
//{
//  double do_kpkt_bb(struct packet *pkt_ptr, double t1, double t2);
//  return do_kpkt_bb(pkt_ptr, t1, t2);
//}
{
  const int tid = get_thread_num();
  const double t1 = pkt_ptr->prop_time;
  const int cellindex = pkt_ptr->where;
  const int modelgridindex = grid::get_cell_modelgridindex(cellindex);

  /// don't calculate cooling rates after each cell crossings anylonger
  /// but only if we really get a kpkt and they hadn't been calculated already
  // if (globals::cellhistory[tid].totalcooling == COOLING_UNDEFINED)
  /*  int ondemand = 1;
    if (grid::modelgrid[modelgridindex].totalcooling == COOLING_UNDEFINED)
    {
      //printout("initial calculation of all cooling rates\n");
      coolingratecalccounter++;
      ondemand = 0;
      calculate_kpkt_rates(modelgridindex);
    }*/
  // printout("totalcooling %g\n",grid::modelgrid[modelgridindex].totalcooling );

  // printout("[debug] do_kpkt: propagate k-pkt\n");
  // double cut = 0.99; //1.;

  const float T_e = grid::get_Te(modelgridindex);
  double deltat = 0.;
  if (nts < globals::n_kpktdiffusion_timesteps) {
    deltat = globals::kpktdiffusion_timescale * globals::time_step[nts].width;
  }
  // double deltat = 1./(nne*1.02e-12*pow(T_e/1e4,0.843));
  // printout("kpkt diffusion time simple %g, advanced %g\n",deltat,1/(nne*1.02e-12*pow(T_e/1e4,0.843)));
  double t_current = t1 + deltat;

  if (t_current <= t2) {
    vec_scale(pkt_ptr->pos, t_current / t1);
    pkt_ptr->prop_time = t_current;

    /// Randomly select the occuring cooling process out of the important ones
    double coolingsum = 0.;
    double zrand = gsl_rng_uniform(rng);

    const double rndcool = zrand * grid::modelgrid[modelgridindex].totalcooling;
    // printout("rndcool %g totalcooling %g\n",rndcool, grid::modelgrid[modelgridindex].totalcooling);
    double oldcoolingsum = 0.;
    int element = -1;
    int ion = -1;
    for (element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (ion = 0; ion < nions; ion++) {
        oldcoolingsum = coolingsum;
        coolingsum += grid::modelgrid[modelgridindex].cooling_contrib_ion[element][ion];
        // printout("Z=%d, ionstage %d, coolingsum %g\n", get_element(element), get_ionstage(element, ion), coolingsum);
        if (coolingsum > rndcool) break;
      }
      if (coolingsum > rndcool) break;
    }
    // printout("kpkt selected Z=%d ionstage %d\n", get_element(element), get_ionstage(element, ion));

    if (element >= get_nelements() || ion >= get_nions(element)) {
      printout("do_kpkt: problem selecting a cooling process ... abort\n");
      printout("do_kpkt: modelgridindex %d element %d ion %d\n", modelgridindex, element, ion);
      printout("do_kpkt: totalcooling %g, coolingsum %g, rndcool %g\n", grid::modelgrid[modelgridindex].totalcooling,
               coolingsum, rndcool);
      printout("do_kpkt: modelgridindex %d, cellno %d, nne %g\n", modelgridindex, pkt_ptr->where,
               grid::get_nne(modelgridindex));
      for (element = 0; element < get_nelements(); element++) {
        const int nions = get_nions(element);
        for (ion = 0; ion < nions; ion++) {
          printout("do_kpkt: element %d, ion %d, coolingcontr %g\n", element, ion,
                   grid::modelgrid[modelgridindex].cooling_contrib_ion[element][ion]);
        }
      }
      abort();
    }

    // globals::debuglevel = 2;
    // printout("element %d, ion %d, coolingsum %g\n",element,ion,coolingsum);
    const int ilow = get_coolinglistoffset(element, ion);
    int low = ilow;
    int high = low + get_ncoolingterms(element, ion) - 1;
    // printout("element %d, ion %d, low %d, high %d\n",element,ion,low,high);
    if (globals::cellhistory[tid].cooling_contrib[ilow] < 0.) {
      // printout("calculate kpkt rates on demand modelgridindex %d element %d ion %d ilow %d ihigh %d oldcoolingsum
      // %g\n",
      //          modelgridindex, element, ion, ilow, high, oldcoolingsum);
      calculate_kpkt_rates_ion(modelgridindex, element, ion, ilow, oldcoolingsum, tid);
    }

    int i = -1;
    while (low <= high) {
      i = (low + high) / 2;
      if (globals::cellhistory[tid].cooling_contrib[i] >= rndcool) {
        if ((i == ilow) || ((i > 0 ? globals::cellhistory[tid].cooling_contrib[i - 1] : 0.) < rndcool))
          break;  /// found (1)
        else
          high = i - 1;
      } else
        low = i + 1;
      // else if (globals::cellhistory[tid].cooling_contrib < rndcool)
      //   low = i + 1;
      // else
      //   break; /// found (2)
    }
    // random value minus
    if (low > high) {
      printout("do_kpkt: error occured while selecting a cooling channel: low %d, high %d, i %d, rndcool %g\n", low,
               high, i, rndcool);
      printout("element %d, ion %d, offset %d, terms %d, coolingsum %g\n", element, ion,
               get_coolinglistoffset(element, ion), get_ncoolingterms(element, ion), coolingsum);
      printout("oldcoolingsum %g, coolingsum %g\n", oldcoolingsum, coolingsum);

      printout("lower %g, %g, %g\n", globals::cellhistory[tid].cooling_contrib[get_coolinglistoffset(element, ion) - 1],
               globals::cellhistory[tid].cooling_contrib[get_coolinglistoffset(element, ion)],
               globals::cellhistory[tid].cooling_contrib[get_coolinglistoffset(element, ion) + 1]);
      int finalpos = get_coolinglistoffset(element, ion) + get_ncoolingterms(element, ion) - 1;
      printout("upper %g, %g, %g\n", globals::cellhistory[tid].cooling_contrib[finalpos - 1],
               globals::cellhistory[tid].cooling_contrib[finalpos],
               globals::cellhistory[tid].cooling_contrib[finalpos + 1]);
    }

    // if (globals::debuglevel == 2) printout("do_kpkt: selected process %d, coolingsum %g\n",i,coolingsum);
    // printout("do_kpkt: selected process %d, coolingsum %g, importantcoolingterms %d, type
    // %d\n",i,coolingsum,importantcoolingterms,globals::cellhistory[tid].coolinglist[i].type);

    // printout("element Z=%d, ion_stage %d, leve %d upper %d offset %d, terms %d, coolingsum %g\n",
    //   get_element(globals::cellhistory[tid].coolinglist[i].element),
    //   get_ionstage(globals::cellhistory[tid].coolinglist[i].element, globals::cellhistory[tid].coolinglist[i].ion),
    //   globals::cellhistory[tid].coolinglist[i].level,
    //   globals::cellhistory[tid].coolinglist[i].upperlevel,
    //   get_coolinglistoffset(element,ion), get_ncoolingterms(element,ion), coolingsum);

    if (coolinglist[i].type == COOLINGTYPE_FF) {
      /// The k-packet converts directly into a r-packet by free-free-emission.
      /// Need to select the r-packets frequency and a random direction in the
      /// co-moving frame.
      // printout("[debug] do_kpkt: k-pkt -> free-free\n");
      // kdecay.to_r++;

      /// Sample the packets comoving frame frequency according to paperII 5.4.3 eq.41
      // zrand = gsl_rng_uniform(rng);   /// delivers zrand in [0,1[
      // zrand = 1. - zrand;             /// make sure that log gets a zrand in ]0,1]
      zrand = gsl_rng_uniform_pos(rng);  /// delivers zrand in ]0,1[
      pkt_ptr->nu_cmf = -KB * T_e / H * log(zrand);

      if (!std::isfinite(pkt_ptr->nu_cmf)) {
        printout("[fatal] ff cooling: selected frequency not finite ... abort\n");
        abort();
      }
      /// and then emitt the packet randomly in the comoving frame
      emitt_rpkt(pkt_ptr);
      pkt_ptr->next_trans = 0;  /// FLAG: transition history here not important, cont. process
      // if (tid == 0) k_stat_to_r_ff++;
      stats::increment(stats::COUNTER_K_STAT_TO_R_FF);

      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 6;
      pkt_ptr->emissiontype = -9999999;
      vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
      pkt_ptr->em_time = pkt_ptr->prop_time;
      pkt_ptr->nscatterings = 0;

/* call the estimator routine - generate a virtual packet */
#ifdef VPKT_ON
      int realtype = 2;
      vpkt_call_estimators(pkt_ptr, t_current, realtype);
#endif
    } else if (coolinglist[i].type == COOLINGTYPE_FB) {
      /// The k-packet converts directly into a r-packet by free-bound-emission.
      /// Need to select the r-packets frequency and a random direction in the
      /// co-moving frame.
      const int element = coolinglist[i].element;
      const int lowerion = coolinglist[i].ion;
      const int level = coolinglist[i].level;
      const int upper = coolinglist[i].upperlevel;
      // const double nu_threshold = get_phixs_threshold(element, ion, level, phixstargetindex)

      // printout("[debug] do_kpkt: k-pkt -> free-bound\n");
      // printout("[debug] do_kpkt: element  %d, ion %d, level %d, upper %d, nu_threshold
      // %g\n",element,ion,level,upper,nu_threshold);

      /// then randomly sample the packets frequency according to the continuums
      /// energy distribution and set some flags
      // zrand = gsl_rng_uniform(rng);   /// delivers zrand in [0,1[
      // zrand = 1. - zrand;   /// convert it to ]0,1]
      // pkt_ptr->nu_cmf = nu_threshold * (1 - KB*T_e/H/nu_threshold*log(zrand));
      // pkt_ptr->nu_cmf = nu_threshold * (1+sqrt(1+(4*KB*T_e/H/nu_threshold)))/2 * (1 -
      // KB*T_e/H/nu_threshold*log(zrand)); pkt_ptr->nu_cmf = nu_threshold;

      // Sample the packets comoving frame frequency according to paperII 4.2.2
      // zrand = gsl_rng_uniform(rng);
      // if (zrand < 0.5)
      { pkt_ptr->nu_cmf = select_continuum_nu(element, lowerion, level, upper, T_e); }
      // else
      // {
      //   ///Emitt like a BB
      //   pkt_ptr->nu_cmf = sample_planck(T_e);
      // }

      // printout("[debug] do_kpkt: pkt_ptr->nu_cmf %g\n",pkt_ptr->nu_cmf);

      // and then emitt the packet randomly in the comoving frame
      emitt_rpkt(pkt_ptr);

#if (TRACK_ION_STATS)
      stats::increment_ion_stats(modelgridindex, element, lowerion + 1, stats::ION_RADRECOMB_KPKT,
                                 pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf);
      const double escape_prob = get_rpkt_escape_prob(pkt_ptr, pkt_ptr->prop_time);
      stats::increment_ion_stats(modelgridindex, element, lowerion + 1, stats::ION_RADRECOMB_ESCAPED,
                                 pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf * escape_prob);
#endif

      pkt_ptr->next_trans = 0;  /// FLAG: transition history here not important, cont. process
      stats::increment(stats::COUNTER_K_STAT_TO_R_FB);
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 7;
      pkt_ptr->emissiontype = get_continuumindex(element, lowerion, level, upper);
      pkt_ptr->trueemissiontype = pkt_ptr->emissiontype;
      vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
      pkt_ptr->em_time = pkt_ptr->prop_time;
      pkt_ptr->nscatterings = 0;

// call the estimator routine - generate a virtual packet
#ifdef VPKT_ON
      int realtype = 2;
      vpkt_call_estimators(pkt_ptr, t_current, realtype);
#endif
    } else if (coolinglist[i].type == COOLINGTYPE_COLLEXC) {
      /// the k-packet activates a macro-atom due to collisional excitation
      // printout("[debug] do_kpkt: k-pkt -> collisional excitation of MA\n");
      const float nne = grid::get_nne(modelgridindex);

      // if the previous entry belongs to a different ion, it might have an invalid (uncalculated)
      // cooling_contrib, so need to use oldcoolingsum (coolingsum up to the start of the current ion)
      const double contrib_low = (i > ilow) ? globals::cellhistory[tid].cooling_contrib[i - 1] : oldcoolingsum;

      double contrib = contrib_low;
      assert_testmodeonly(coolinglist[i].element == element);
      assert_testmodeonly(coolinglist[i].ion == ion);
      const int level = coolinglist[i].level;
      const double epsilon_current = epsilon(element, ion, level);
      const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
      int upper = -1;
      // excitation to same ionization stage
      const int nuptrans = get_nuptrans(element, ion, level);
      for (int ii = 0; ii < nuptrans; ii++) {
        const int lineindex = globals::elements[element].ions[ion].levels[level].uptrans_lineindicies[ii];
        const int tmpupper = globals::linelist[lineindex].upperlevelindex;
        // printout("    excitation to level %d possible\n",upper);
        const double epsilon_trans = epsilon(element, ion, tmpupper) - epsilon_current;
        const double C = nnlevel * col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans) * epsilon_trans;
        contrib += C;
        if (contrib >= rndcool) {
          upper = tmpupper;
          break;
        }
      }

      if (upper < 0) {
        printout(
            "WARNING: Could not select an upper level. modelgridindex %d i %d element %d ion %d level %d rndcool %g "
            "contrib_low %g contrib %g (should match %g) upper %d nuptrans %d\n",
            modelgridindex, i, element, ion, level, rndcool, contrib_low, contrib,
            globals::cellhistory[tid].cooling_contrib[i], upper, nuptrans);
        abort();
      }
      assert_always(upper >= 0);

      const int element = coolinglist[i].element;
      const int ion = coolinglist[i].ion;
      // const int upper = coolinglist[i].upperlevel;
      pkt_ptr->mastate.element = element;
      pkt_ptr->mastate.ion = ion;
      pkt_ptr->mastate.level = upper;
      pkt_ptr->mastate.activatingline = -99;

      if (TRACK_ION_STATS) {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYIN_COLLEXC, pkt_ptr->e_cmf);
      }

      pkt_ptr->type = TYPE_MA;
      stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_COLLEXC);
      stats::increment(stats::COUNTER_K_STAT_TO_MA_COLLEXC);

      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 8;
      pkt_ptr->trueemissiontype = -1;  // since this is below zero, macroatom will set it
      pkt_ptr->trueemissionvelocity = -1;
    } else if (coolinglist[i].type == COOLINGTYPE_COLLION) {
      /// the k-packet activates a macro-atom due to collisional ionisation
      if (globals::debuglevel == 2) printout("[debug] do_kpkt: k-pkt -> collisional ionisation of MA\n");
      const int element = coolinglist[i].element;
      const int ion = coolinglist[i].ion + 1;
      const int upper = coolinglist[i].upperlevel;
      pkt_ptr->mastate.element = element;
      pkt_ptr->mastate.ion = ion;
      pkt_ptr->mastate.level = upper;
      pkt_ptr->mastate.activatingline = -99;

      if (TRACK_ION_STATS) {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYIN_COLLION, pkt_ptr->e_cmf);
      }

      pkt_ptr->type = TYPE_MA;
      stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_COLLION);
      stats::increment(stats::COUNTER_K_STAT_TO_MA_COLLION);
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 9;
      pkt_ptr->trueemissiontype = -1;  // since this is below zero, macroatom will set it
      pkt_ptr->trueemissionvelocity = -1;
    } else {
      printout("[fatal] do_kpkt: coolinglist.type mismatch\n");
      printout("[fatal] do_kpkt: zrand %g, grid::modelgrid[modelgridindex].totalcooling %g, coolingsum %g, i %d\n",
               zrand, grid::modelgrid[modelgridindex].totalcooling, coolingsum, i);
      printout("[fatal] do_kpkt: coolinglist[i].type %d\n", coolinglist[i].type);
      printout("[fatal] do_kpkt: pkt_ptr->where %d, mgi %d\n", pkt_ptr->where, modelgridindex);
      abort();
    }

    return pkt_ptr->prop_time;
  } else {
    vec_scale(pkt_ptr->pos, t2 / t1);
    pkt_ptr->prop_time = t2;
    return pkt_ptr->prop_time;
  }
}

/*static int compare_coolinglistentry(const void *p1, const void *p2)
/// Helper function to sort the coolinglist by the strength of the
/// individual cooling contributions.
{
  ionscoolinglist_t *a1, *a2;
  a1 = (ionscoolinglist_t *)(p1);
  a2 = (ionscoolinglist_t *)(p2);
  //printf("%d %d %d %d %g\n",a1->elementindex,a1->ionindex,a1->lowerlevelindex,a1->upperlevelindex,a1->nu);
  //printf("%d %d %d %d %g\n",a2->elementindex,a2->ionindex,a2->lowerlevelindex,a2->upperlevelindex,a2->nu);
  //printf("%g\n",a2->nu - a1->nu);
  if (a1->contribution - a2->contribution < 0)
    return 1;
  else if (a1->contribution - a2->contribution > 0)
    return -1;
  else
    return 0;
}*/

/*double get_bfcooling_direct(int element, int ion, int level, int cellnumber)
/// Returns the rate for bfheating. This can be called during packet propagation
/// or update_grid. Therefore we need to decide whether a cell history is
/// known or not.
{
  double bfcooling;

  gsl_integration_workspace *wsp;
  gslintegration_paras intparas;
  double bfcooling_integrand_gsl(double nu, void *paras);
  gsl_function F_bfcooling;
  F_bfcooling.function = &bfcooling_integrand_gsl;
  double intaccuracy = 1e-2;        /// Fractional accuracy of the integrator
  double error;
  double nu_max_phixs;

  float T_e = globals::cell[cellnumber].T_e;
  float nne = globals::cell[cellnumber].nne;
  double nnionlevel = get_groundlevelpop(cellnumber,element,ion+1);
  //upper = coolinglist[i].upperlevel;
  double nu_threshold = (epsilon(element,ion+1,0) - epsilon(element,ion,level)) / H;
  nu_max_phixs = nu_threshold * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table

  pkt_ptr->mastate.element = element;
  pkt_ptr->mastate.ion = ion;
  pkt_ptr->mastate.level = level;
  intparas.T = T_e;
  intparas.nu_edge = nu_threshold;   /// Global variable which passes the threshold to the integrator
  F_bfcooling.params = &intparas;
  gsl_integration_qag(&F_bfcooling, nu_threshold, nu_max_phixs, 0, intaccuracy, 1024, 6, wsp, &bfcooling, &error);
  bfcooling *= nnionlevel*nne*4*PI*calculate_sahafact(element,ion,level,upperionlevel,T_e,nu_threshold*H);

  return bfcooling;
}*/

}  // namespace kpkt