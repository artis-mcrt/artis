#include "kpkt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "macroatom.h"
#include "packet.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "thermalbalance.h"
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
  int level;
  int upperlevel;
};

static struct cellhistorycoolinglist *coolinglist;

static auto get_ncoolingterms_ion(int element, int ion) -> int {
  return globals::elements[element].ions[ion].ncoolingterms;
}

template <bool update_cooling_contrib_list>
static auto calculate_cooling_rates_ion(const int modelgridindex, const int element, const int ion,
                                        const int indexionstart, const int tid, double *C_ff, double *C_fb,
                                        double *C_exc, double *C_ionization) -> double
// calculate the cooling contribution list of individual levels/processes for an ion
// oldcoolingsum is the sum of lower ion (of same element or all ions of lower elements) cooling contributions
{
  const auto nne = grid::get_nne(modelgridindex);
  const auto T_e = grid::get_Te(modelgridindex);

  double C_ion = 0.;
  int i = indexionstart;  // NOLINT(misc-const-correctness)

  const int nions = get_nions(element);

  const int nionisinglevels = get_ionisinglevels(element, ion);
  const double nncurrention = get_nnion(modelgridindex, element, ion);

  /// ff creation of rpkt
  const int ioncharge = get_ionstage(element, ion) - 1;
  // printout("[debug] ioncharge %d, nncurrention %g, nne %g\n",ion,nncurrention,nne);
  if (ioncharge > 0) {
    const double C_ff_ion = 1.426e-27 * sqrt(T_e) * pow(ioncharge, 2) * nncurrention * nne;
    C_ion += C_ff_ion;

    if constexpr (update_cooling_contrib_list) {
      globals::cellhistory[tid].cooling_contrib[i] = C_ion;

      assert_testmodeonly(coolinglist[i].type == COOLINGTYPE_FF);

      // if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d,
      // coolingtype %d, i %d, low
      // %d\n",contrib,oldcoolingsum,C,element,ion,-99,globals::cellhistory[tid].coolinglist[i].type,i,low);
      i++;
    } else {
      *C_ff += C_ff_ion;
    }
  }

  /// excitation to same ionization stage
  const int nlevels = get_nlevels(element, ion);
  for (int level = 0; level < nlevels; level++) {
    // printout("[debug] do_kpkt: element %d, ion %d, level %d\n", element, ion, level);
    const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
    const double epsilon_current = epsilon(element, ion, level);
    const double statweight = stat_weight(element, ion, level);

    const int nuptrans = get_nuptrans(element, ion, level);
    if (nuptrans > 0) {
      for (int ii = 0; ii < nuptrans; ii++) {
        const int upper = globals::elements[element].ions[ion].levels[level].uptrans[ii].targetlevelindex;
        const double epsilon_trans = epsilon(element, ion, upper) - epsilon_current;
        const double C = nnlevel *
                         col_excitation_ratecoeff(T_e, nne, element, ion, level, ii, epsilon_trans, statweight) *
                         epsilon_trans;
        C_ion += C;
        if constexpr (!update_cooling_contrib_list) {
          *C_exc += C;
        }
      }
      if constexpr (update_cooling_contrib_list) {
        globals::cellhistory[tid].cooling_contrib[i] = C_ion;

        assert_testmodeonly(coolinglist[i].type == COOLINGTYPE_COLLEXC);

        i++;
      }
    }
  }

  if (ion < nions - 1) {
    const double nnupperion = get_nnion(modelgridindex, element, ion + 1);

    // ionization to higher ionization stage
    for (int level = 0; level < nionisinglevels; level++) {
      const double epsilon_current = epsilon(element, ion, level);
      const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
      const int nphixstargets = get_nphixstargets(element, ion, level);
      for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
        const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
        const double epsilon_upper = epsilon(element, ion + 1, upper);
        const double epsilon_trans = epsilon_upper - epsilon_current;
        const double C = nnlevel *
                         col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans) *
                         epsilon_trans;

        C_ion += C;
        if constexpr (update_cooling_contrib_list) {
          globals::cellhistory[tid].cooling_contrib[i] = C_ion;

          assert_testmodeonly(coolinglist[i].type == COOLINGTYPE_COLLION);
          assert_testmodeonly(coolinglist[i].level == level);
          assert_testmodeonly(coolinglist[i].upperlevel == upper);

          i++;
        } else {
          *C_ionization += C;
        }
      }
    }

    /// fb creation of r-pkt
    /// free bound rates are calculated from the lower ion, but associated to the higher ion
    for (int level = 0; level < nionisinglevels; level++) {
      const int nphixstargets = get_nphixstargets(element, ion, level);
      for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
        const double pop =
            (BFCOOLING_USELEVELPOPNOTIONPOP ? get_levelpop(modelgridindex, element, ion + 1,
                                                           get_phixsupperlevel(element, ion, level, phixstargetindex))
                                            : nnupperion);
        const double C = get_bfcoolingcoeff(element, ion, level, phixstargetindex, T_e) * pop * nne;
        C_ion += C;

        if constexpr (update_cooling_contrib_list) {
          globals::cellhistory[tid].cooling_contrib[i] = C_ion;

          assert_testmodeonly(coolinglist[i].type == COOLINGTYPE_FB);
          assert_testmodeonly(coolinglist[i].level == level);
          assert_testmodeonly(coolinglist[i].upperlevel == get_phixsupperlevel(element, ion, level, phixstargetindex));

          i++;
        } else {
          *C_fb += C;
        }
      }
    }
  }

  if constexpr (update_cooling_contrib_list) {
    assert_testmodeonly(indexionstart == get_coolinglistoffset(element, ion));
    assert_always(i == indexionstart + get_ncoolingterms_ion(element, ion));
  }

  return C_ion;
}

void calculate_cooling_rates(const int modelgridindex, struct heatingcoolingrates *heatingcoolingrates)
// Calculate the cooling rates for a given cell and store them for each ion
// optionally store components (ff, bf, collisional) in heatingcoolingrates struct
{
  double C_ff_all = 0.;          /// free-free creation of rpkts
  double C_fb_all = 0.;          /// free-bound creation of rpkt
  double C_exc_all = 0.;         /// collisional excitation of macroatoms
  double C_ionization_all = 0.;  /// collisional ionisation of macroatoms
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const double C_ion = calculate_cooling_rates_ion<false>(modelgridindex, element, ion, -1, tid, &C_ff_all,
                                                              &C_fb_all, &C_exc_all, &C_ionization_all);
      grid::modelgrid[modelgridindex].cooling_contrib_ion[element][ion] = C_ion;
    }
  }

  // this loop is made separate for future parallelisation of upper loop.
  // the ion contributions must be added in this exact order
  double C_total = 0.;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      C_total += grid::modelgrid[modelgridindex].cooling_contrib_ion[element][ion];
    }
  }
  grid::modelgrid[modelgridindex].totalcooling = C_total;

  // only used in the T_e solver and write_to_estimators file
  if (heatingcoolingrates != nullptr) {
    heatingcoolingrates->cooling_collisional = C_exc_all + C_ionization_all;
    heatingcoolingrates->cooling_fb = C_fb_all;
    heatingcoolingrates->cooling_ff = C_ff_all;
  }
}

static void set_ncoolingterms() {
  globals::ncoolingterms = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      int ionterms = 0;
      globals::elements[element].ions[ion].coolingoffset = globals::ncoolingterms;

      /// Ionised ions add one ff-cooling term
      if (get_ionstage(element, ion) > 1) {
        ionterms++;
      }
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

void setup_coolinglist() {
  /// SET UP THE COOLING LIST
  ///======================================================
  /// Determine number of processes which allow kpkts to convert to something else.
  /// This number is given by the collisional excitations (so far determined from the oscillator strengths
  /// by the van Regemorter formula, therefore totaluptrans), the number of free-bound emissions and collisional
  /// ionisations (as long as we only deal with ionisation to the ground level this means for both of these
  /// \sum_{elements,ions}get_nlevels(element,ion) and free-free which is \sum_{elements} get_nions(element)-1

  set_ncoolingterms();
  const size_t mem_usage_coolinglist = globals::ncoolingterms * sizeof(struct cellhistorycoolinglist);
  coolinglist = static_cast<struct cellhistorycoolinglist *>(
      malloc(globals::ncoolingterms * sizeof(struct cellhistorycoolinglist)));
  printout("[info] mem_usage: coolinglist occupies %.3f MB\n", mem_usage_coolinglist / 1024. / 1024.);

  int i = 0;  // cooling list index
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels_currention = get_nlevels(element, ion);

      const int nionisinglevels = get_ionisinglevels(element, ion);

      /// ff creation of rpkt
      /// -------------------
      const int ioncharge = get_ionstage(element, ion) - 1;
      // printout("[debug] ioncharge %d, nncurrention %g, nne %g\n",ion,nncurrention,nne);
      if (ioncharge > 0) {
        coolinglist[i].type = COOLINGTYPE_FF;
        coolinglist[i].level = -99;
        coolinglist[i].upperlevel = -99;
        i++;
      }

      for (int level = 0; level < nlevels_currention; level++) {
        if (get_nuptrans(element, ion, level) > 0) {
          coolinglist[i].type = COOLINGTYPE_COLLEXC;
          coolinglist[i].level = level;
          // upper level is not valid because this is the contribution of all upper levels combined - have to
          // calculate individually when selecting a random process
          coolinglist[i].upperlevel = -1;
          i++;
        }
      }

      if (ion < (nions - 1))  /// check whether further ionisation stage available
      {
        for (int level = 0; level < nionisinglevels; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
            coolinglist[i].type = COOLINGTYPE_COLLION;
            coolinglist[i].level = level;
            coolinglist[i].upperlevel = upper;
            i++;
          }
        }

        /// fb creation of r-pkt
        /// free bound rates are calculated from the lower ion, but associated to the higher ion
        for (int level = 0; level < nionisinglevels; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
            coolinglist[i].type = COOLINGTYPE_FB;
            coolinglist[i].level = level;
            coolinglist[i].upperlevel = upper;
            i++;
          }
        }
      }
      assert_always(i == get_coolinglistoffset(element, ion) + get_ncoolingterms_ion(element, ion));
    }
  }

  assert_always(globals::ncoolingterms == i);  // if this doesn't match, we miscalculated the number of cooling terms
  printout("[info] read_atomicdata: number of coolingterms %d\n", globals::ncoolingterms);
}

static auto sample_planck(const double T) -> double
/// returns a randomly chosen frequency according to the Planck
/// distribution of temperature T
{
  const double nu_peak = 5.879e10 * T;
  if (nu_peak > NU_MAX_R || nu_peak < NU_MIN_R) {
    printout("[warning] sample_planck: intensity peaks outside frequency range\n");
  }

  const double B_peak = radfield::dbb(nu_peak, T, 1);

  while (true) {
    const double nu = NU_MIN_R + rng_uniform() * (NU_MAX_R - NU_MIN_R);
    if (rng_uniform() * B_peak <= radfield::dbb(nu, T, 1)) {
      return nu;
    }
    // printout("[debug] sample_planck: planck_sampling %d\n", i);
  }
}

void do_kpkt_blackbody(struct packet *pkt_ptr)
/// handle a k-packet (e.g., in a thick cell) by emitting according to the planck function
{
  const int modelgridindex = grid::get_cell_modelgridindex(pkt_ptr->where);

  pkt_ptr->nu_cmf = sample_planck(grid::get_Te(modelgridindex));
  assert_always(std::isfinite(pkt_ptr->nu_cmf));
  /// and then emitt the packet randomly in the comoving frame
  emit_rpkt(pkt_ptr);
  // printout("[debug] calculate_chi_rpkt after kpkt to rpkt by ff\n");
  pkt_ptr->next_trans = 0;  /// FLAG: transition history here not important, cont. process
  // if (tid == 0) k_stat_to_r_bb++;
  stats::increment(stats::COUNTER_K_STAT_TO_R_BB);
  pkt_ptr->interactions++;
  pkt_ptr->last_event = 6;
  pkt_ptr->emissiontype = EMTYPE_FREEFREE;
  vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
  pkt_ptr->em_time = pkt_ptr->prop_time;
  pkt_ptr->nscatterings = 0;
}

void do_kpkt(struct packet *pkt_ptr, double t2, int nts)
/// handle a k-packet (kinetic energy of the free electrons)
{
  const int tid = get_thread_num();
  const double t1 = pkt_ptr->prop_time;
  const int modelgridindex = grid::get_cell_modelgridindex(pkt_ptr->where);

  /// don't calculate cooling rates after each cell crossings any longer
  /// but only if we really get a kpkt and they hadn't been calculated already

  // printout("[debug] do_kpkt: propagate k-pkt\n");

  const auto T_e = grid::get_Te(modelgridindex);
  double deltat = 0.;
  if (nts < globals::n_kpktdiffusion_timesteps) {
    deltat = globals::kpktdiffusion_timescale * globals::timesteps[nts].width;
  }
  // double deltat = 1. / (nne * 1.02e-12 * pow(T_e / 1e4, 0.843));
  // printout("kpkt diffusion time simple %g, advanced %g\n", deltat, 1 / (nne * 1.02e-12 * pow(T_e / 1e4, 0.843)));
  const double t_current = t1 + deltat;

  if (t_current > t2) {
    vec_scale(pkt_ptr->pos, t2 / t1);
    pkt_ptr->prop_time = t2;
  } else {
    vec_scale(pkt_ptr->pos, t_current / t1);
    pkt_ptr->prop_time = t_current;

    /// Randomly select the occuring cooling process
    double coolingsum = 0.;

    assert_always(grid::modelgrid[modelgridindex].totalcooling > 0.);
    const double rndcool_ion = rng_uniform() * grid::modelgrid[modelgridindex].totalcooling;

    int element = -1;
    int ion = -1;
    for (element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (ion = 0; ion < nions; ion++) {
        coolingsum += grid::modelgrid[modelgridindex].cooling_contrib_ion[element][ion];
        // printout("Z=%d, ionstage %d, coolingsum %g\n", get_atomicnumber(element), get_ionstage(element, ion),
        // coolingsum);
        if (coolingsum > rndcool_ion) {
          break;
        }
      }
      if (coolingsum > rndcool_ion) {
        break;
      }
    }
    // printout("kpkt selected Z=%d ionstage %d\n", get_atomicnumber(element), get_ionstage(element, ion));

    if (element >= get_nelements() || element < 0 || ion >= get_nions(element) || ion < 0) {
      printout("do_kpkt: problem selecting a cooling process ... abort\n");
      printout("do_kpkt: modelgridindex %d element %d ion %d\n", modelgridindex, element, ion);
      printout("do_kpkt: totalcooling %g, coolingsum %g, rndcool_ion %g\n",
               grid::modelgrid[modelgridindex].totalcooling, coolingsum, rndcool_ion);
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

    // printout("element %d, ion %d, coolingsum %g\n",element,ion,coolingsum);
    const int ilow = get_coolinglistoffset(element, ion);
    const int ihigh = ilow + get_ncoolingterms_ion(element, ion) - 1;
    // printout("element %d, ion %d, low %d, high %d\n",element,ion,low,high);
    if (globals::cellhistory[tid].cooling_contrib[ilow] < 0.) {
      // printout("calculate kpkt rates on demand modelgridindex %d element %d ion %d ilow %d ihigh %d
      // oldcoolingsum %g\n",
      //          modelgridindex, element, ion, ilow, high, oldcoolingsum);
      const double C_ion = calculate_cooling_rates_ion<true>(modelgridindex, element, ion, ilow, tid, nullptr, nullptr,
                                                             nullptr, nullptr);
      // we just summed up every individual cooling process. make sure it matches the stored total for the ion
      assert_always(C_ion == grid::modelgrid[modelgridindex].cooling_contrib_ion[element][ion]);
    }

    // with the ion selected, we now select a level and transition type

    const double zrand2 = rng_uniform();
    const double rndcool_ion_process = zrand2 * globals::cellhistory[tid].cooling_contrib[ihigh];

    auto *const selectedvalue =
        std::upper_bound(&globals::cellhistory[tid].cooling_contrib[ilow],
                         &globals::cellhistory[tid].cooling_contrib[ihigh + 1], rndcool_ion_process);
    const ptrdiff_t i = selectedvalue - globals::cellhistory[tid].cooling_contrib;

    if (i > ihigh) {
      printout("do_kpkt: error occured while selecting a cooling channel: low %d, high %d, i %d, rndcool %g\n", ilow,
               ihigh, i, rndcool_ion_process);
      printout("element %d, ion %d, offset %d, terms %d, coolingsum %g\n", element, ion,
               get_coolinglistoffset(element, ion), get_ncoolingterms_ion(element, ion), coolingsum);

      printout("lower %g, %g, %g\n", globals::cellhistory[tid].cooling_contrib[get_coolinglistoffset(element, ion) - 1],
               globals::cellhistory[tid].cooling_contrib[get_coolinglistoffset(element, ion)],
               globals::cellhistory[tid].cooling_contrib[get_coolinglistoffset(element, ion) + 1]);
      const int finalpos = get_coolinglistoffset(element, ion) + get_ncoolingterms_ion(element, ion) - 1;
      printout("upper %g, %g, %g\n", globals::cellhistory[tid].cooling_contrib[finalpos - 1],
               globals::cellhistory[tid].cooling_contrib[finalpos],
               globals::cellhistory[tid].cooling_contrib[finalpos + 1]);
    }

    assert_always(i <= ihigh);

    // printout("do_kpkt: selected process %d, coolingsum %g\n", i, coolingsum);
    const auto rndcoolingtype = coolinglist[i].type;
    if (rndcoolingtype == COOLINGTYPE_FF) {
      /// The k-packet converts directly into a r-packet by free-free-emission.
      /// Need to select the r-packets frequency and a random direction in the
      /// co-moving frame.
      // printout("[debug] do_kpkt: k-pkt -> free-free\n");

      /// Sample the packets comoving frame frequency according to paperII 5.4.3 eq.41
      // zrand = rng_uniform();   /// delivers zrand in [0,1[
      // zrand = 1. - zrand;             /// make sure that log gets a zrand in ]0,1]
      const double zrand = rng_uniform_pos();  /// delivers zrand in ]0,1[
      pkt_ptr->nu_cmf = -KB * T_e / H * log(zrand);

      assert_always(std::isfinite(pkt_ptr->nu_cmf));
      /// and then emitt the packet randomly in the comoving frame
      emit_rpkt(pkt_ptr);
      pkt_ptr->next_trans = 0;  /// FLAG: transition history here not important, cont. process
      stats::increment(stats::COUNTER_K_STAT_TO_R_FF);

      pkt_ptr->last_event = 6;
      pkt_ptr->emissiontype = EMTYPE_FREEFREE;
      vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
      pkt_ptr->em_time = pkt_ptr->prop_time;
      pkt_ptr->nscatterings = 0;

      vpkt_call_estimators(pkt_ptr, TYPE_KPKT);
    } else if (rndcoolingtype == COOLINGTYPE_FB) {
      /// The k-packet converts directly into a r-packet by free-bound-emission.
      /// Need to select the r-packets frequency and a random direction in the
      /// co-moving frame.
      const int lowerion = ion;
      const int lowerlevel = coolinglist[i].level;
      const int upper = coolinglist[i].upperlevel;

      /// then randomly sample the packets frequency according to the continuums
      /// energy distribution

      // Sample the packets comoving frame frequency according to paperII 4.2.2
      // const double zrand = rng_uniform();
      // if (zrand < 0.5) {
      pkt_ptr->nu_cmf = select_continuum_nu(element, lowerion, lowerlevel, upper, T_e);
      // } else {
      //   // Emit like a BB
      //   pkt_ptr->nu_cmf = sample_planck(T_e);
      // }

      // and then emitt the packet randomly in the comoving frame
      emit_rpkt(pkt_ptr);

      if constexpr (TRACK_ION_STATS) {
        stats::increment_ion_stats(modelgridindex, element, lowerion + 1, stats::ION_RADRECOMB_KPKT,
                                   pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf);
      }

      pkt_ptr->next_trans = 0;  /// FLAG: transition history here not important, cont. process
      stats::increment(stats::COUNTER_K_STAT_TO_R_FB);
      pkt_ptr->last_event = 7;
      pkt_ptr->emissiontype = get_continuumindex(element, lowerion, lowerlevel, upper);
      pkt_ptr->trueemissiontype = pkt_ptr->emissiontype;
      vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
      pkt_ptr->em_time = pkt_ptr->prop_time;
      pkt_ptr->nscatterings = 0;

      vpkt_call_estimators(pkt_ptr, TYPE_KPKT);
    } else if (rndcoolingtype == COOLINGTYPE_COLLEXC) {
      /// the k-packet activates a macro-atom due to collisional excitation
      // printout("[debug] do_kpkt: k-pkt -> collisional excitation of MA\n");
      const float nne = grid::get_nne(modelgridindex);

      // if the previous entry belongs to the same ion, then pick up the cumulative sum from
      // the previous entry, otherwise start from zero
      const double contrib_low = (i > ilow) ? globals::cellhistory[tid].cooling_contrib[i - 1] : 0.;

      double contrib = contrib_low;
      const int level = coolinglist[i].level;
      const double epsilon_current = epsilon(element, ion, level);
      const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
      const double statweight = stat_weight(element, ion, level);
      int upper = -1;
      // excitation to same ionization stage
      const int nuptrans = get_nuptrans(element, ion, level);
      const auto *const uptrans = globals::elements[element].ions[ion].levels[level].uptrans;
      for (int ii = 0; ii < nuptrans; ii++) {
        const int tmpupper = uptrans[ii].targetlevelindex;
        // printout("    excitation to level %d possible\n",upper);
        const double epsilon_trans = epsilon(element, ion, tmpupper) - epsilon_current;
        const double C = nnlevel *
                         col_excitation_ratecoeff(T_e, nne, element, ion, level, ii, epsilon_trans, statweight) *
                         epsilon_trans;
        contrib += C;
        if (contrib > rndcool_ion_process) {
          upper = tmpupper;
          break;
        }
      }

      if (upper < 0) {
        printout(
            "WARNING: Could not select an upper level. modelgridindex %d i %d element %d ion %d level %d rndcool "
            "%g "
            "contrib_low %g contrib %g (should match %g) upper %d nuptrans %d\n",
            modelgridindex, i, element, ion, level, rndcool_ion_process, contrib_low, contrib,
            globals::cellhistory[tid].cooling_contrib[i], upper, nuptrans);
        abort();
      }
      assert_always(upper >= 0);

      // const int upper = coolinglist[i].upperlevel;
      pkt_ptr->mastate.element = element;
      pkt_ptr->mastate.ion = ion;
      pkt_ptr->mastate.level = upper;
      pkt_ptr->mastate.activatingline = -99;

      if constexpr (TRACK_ION_STATS) {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYIN_COLLEXC, pkt_ptr->e_cmf);
      }

      pkt_ptr->type = TYPE_MA;
      stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_COLLEXC);
      stats::increment(stats::COUNTER_K_STAT_TO_MA_COLLEXC);

      pkt_ptr->last_event = 8;
      pkt_ptr->trueemissiontype = EMTYPE_NOTSET;
      pkt_ptr->trueemissionvelocity = -1;
    } else if (rndcoolingtype == COOLINGTYPE_COLLION) {
      /// the k-packet activates a macro-atom due to collisional ionisation
      // printout("[debug] do_kpkt: k-pkt -> collisional ionisation of MA\n");

      const int upperion = ion + 1;
      const int upper = coolinglist[i].upperlevel;
      pkt_ptr->mastate.element = element;
      pkt_ptr->mastate.ion = upperion;
      pkt_ptr->mastate.level = upper;
      pkt_ptr->mastate.activatingline = -99;

      if constexpr (TRACK_ION_STATS) {
        stats::increment_ion_stats(modelgridindex, element, upperion, stats::ION_MACROATOM_ENERGYIN_COLLION,
                                   pkt_ptr->e_cmf);
      }

      pkt_ptr->type = TYPE_MA;
      stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_COLLION);
      stats::increment(stats::COUNTER_K_STAT_TO_MA_COLLION);

      pkt_ptr->last_event = 9;
      pkt_ptr->trueemissiontype = EMTYPE_NOTSET;
      pkt_ptr->trueemissionvelocity = -1;
    } else {
      assert_testmodeonly(false);
    }

    pkt_ptr->interactions++;
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