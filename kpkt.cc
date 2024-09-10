#include "kpkt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>

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

namespace {

enum class CoolingType : std::uint8_t { FREEFREE, FREEBOUND, COLLEXC, COLLION };

struct CellCacheCoolingList {
  CoolingType type;
  int level;
  int upperlevel;
};

std::vector<CellCacheCoolingList> coolinglist;

int n_kpktdiffusion_timesteps{0};
float kpktdiffusion_timescale{0.};

// calculate the cooling contribution list of individual levels/processes for an ion
// oldcoolingsum is the sum of lower ion (of same element or all ions of lower elements) cooling contributions
template <bool update_cooling_contrib_list>
auto calculate_cooling_rates_ion(const int modelgridindex, const int element, const int ion, const int indexionstart,
                                 const int cellcacheslotid, double *const C_ff, double *const C_fb, double *const C_exc,
                                 double *const C_ionization) -> double {
  const auto nne = grid::get_nne(modelgridindex);
  const auto T_e = grid::get_Te(modelgridindex);

  if constexpr (update_cooling_contrib_list) {
    assert_always(indexionstart >= 0);
  }

  double C_ion = 0.;
  int i = indexionstart;  // NOLINT(misc-const-correctness)

  const int nionisinglevels = get_ionisinglevels(element, ion);
  const double nncurrention = get_nnion(modelgridindex, element, ion);

  // ff creation of rpkt
  const int ioncharge = get_ionstage(element, ion) - 1;
  // printout("[debug] ioncharge %d, nncurrention %g, nne %g\n",ion,nncurrention,nne);
  if (ioncharge > 0) {
    const double C_ff_ion = 1.426e-27 * sqrt(T_e) * pow(ioncharge, 2) * nncurrention * nne;
    C_ion += C_ff_ion;

    if constexpr (update_cooling_contrib_list) {
      globals::cellcache[cellcacheslotid].cooling_contrib[i] = C_ion;

      assert_testmodeonly(coolinglist[i].type == CoolingType::FREEFREE);

      i++;
    } else {
      *C_ff += C_ff_ion;
    }
  }

  // excitation to same ionization stage
  const int nlevels = get_nlevels(element, ion);
  for (int level = 0; level < nlevels; level++) {
    // printout("[debug] do_kpkt: element %d, ion %d, level %d\n", element, ion, level);
    const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
    const double epsilon_current = epsilon(element, ion, level);
    const double statweight = stat_weight(element, ion, level);

    const int nuptrans = get_nuptrans(element, ion, level);
    if (nuptrans > 0) {
      const auto *const uptranslist = get_uptranslist(element, ion, level);
      for (int ii = 0; ii < nuptrans; ii++) {
        const int upper = uptranslist[ii].targetlevelindex;
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
        globals::cellcache[cellcacheslotid].cooling_contrib[i] = C_ion;

        assert_testmodeonly(coolinglist[i].type == CoolingType::COLLEXC);

        i++;
      }
    }
  }

  if (ion < (get_nions(element) - 1)) {
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
          globals::cellcache[cellcacheslotid].cooling_contrib[i] = C_ion;

          assert_testmodeonly(coolinglist[i].type == CoolingType::COLLION);
          assert_testmodeonly(coolinglist[i].level == level);
          assert_testmodeonly(coolinglist[i].upperlevel == upper);

          i++;
        } else {
          *C_ionization += C;
        }
      }
    }

    // fb creation of r-pkt
    // free bound rates are calculated from the lower ion, but associated to the higher ion
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
          globals::cellcache[cellcacheslotid].cooling_contrib[i] = C_ion;

          assert_testmodeonly(coolinglist[i].type == CoolingType::FREEBOUND);
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

void set_ncoolingterms() {
  ncoolingterms = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      int ionterms = 0;
      globals::elements[element].ions[ion].coolingoffset = ncoolingterms;

      // Ionised ions add one ff-cooling term
      if (get_ionstage(element, ion) > 1) {
        ionterms++;
      }
      // Ionisinglevels below the closure ion add to bf and col ionisation
      // All the levels add number of col excitations
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
      ncoolingterms += ionterms;
    }
  }
}

// return a randomly chosen frequency according to the Planck distribution of temperature T using an analytic method.
// More testing of this function is needed.
auto sample_planck_analytic(const double T) -> double {
  const double nu_peak = 5.879e10 * T;
  if (nu_peak > NU_MAX_R || nu_peak < NU_MIN_R) {
    printout("[warning] sample_planck: intensity peaks outside frequency range\n");
  }

  constexpr ptrdiff_t nubins = 500;
  const auto delta_nu = (NU_MAX_R - NU_MIN_R) / (nubins - 1);
  const auto integral_total = radfield::planck_integral_analytic(T, NU_MIN_R, NU_MAX_R, false);

  const double rand_partintegral = rng_uniform() * integral_total;
  double prev_partintegral = 0.;
  double part_integral = 0.;
  double bin_nu_lower = NU_MIN_R;
  for (ptrdiff_t i = 1; i < nubins; i++) {
    bin_nu_lower = NU_MIN_R + (i - 1) * delta_nu;
    const double nu_upper = NU_MIN_R + (i * delta_nu);
    prev_partintegral = part_integral;
    part_integral = radfield::planck_integral_analytic(T, NU_MIN_R, nu_upper, false);
    if (rand_partintegral >= part_integral) {
      break;
    }
  }

  // use a linear interpolation for the frequency within the bin
  const double nuoffset = (rand_partintegral - prev_partintegral) / (part_integral - prev_partintegral) * delta_nu;

  return bin_nu_lower + nuoffset;
}

// return a randomly chosen frequency according to the Planck distribution of temperature T using a Monte Carlo method
auto sample_planck_montecarlo(const double T) -> double {
  const double nu_peak = 5.879e10 * T;
  if (nu_peak > NU_MAX_R || nu_peak < NU_MIN_R) {
    printout("[warning] sample_planck: intensity peaks outside frequency range\n");
  }

  const double B_peak = radfield::dbb(nu_peak, T, 1);

  while (true) {
    const double nu = NU_MIN_R + (rng_uniform() * (NU_MAX_R - NU_MIN_R));
    if (rng_uniform() * B_peak <= radfield::dbb(nu, T, 1)) {
      return nu;
    }
    // printout("[debug] sample_planck: planck_sampling %d\n", i);
  }
}
}  // anonymous namespace

// Calculate the cooling rates for a given cell and store them for each ion
// optionally store components (ff, bf, collisional) in heatingcoolingrates struct
void calculate_cooling_rates(const int modelgridindex, HeatingCoolingRates *heatingcoolingrates) {
  double C_ff_all = 0.;          // free-free creation of rpkts
  double C_fb_all = 0.;          // free-bound creation of rpkt
  double C_exc_all = 0.;         // collisional excitation of macroatoms
  double C_ionization_all = 0.;  // collisional ionisation of macroatoms
  for (int allionindex = 0; allionindex < get_includedions(); allionindex++) {
    const auto [element, ion] = get_ionfromuniqueionindex(allionindex);
    grid::modelgrid[modelgridindex].ion_cooling_contribs[allionindex] = calculate_cooling_rates_ion<false>(
        modelgridindex, element, ion, -1, cellcacheslotid, &C_ff_all, &C_fb_all, &C_exc_all, &C_ionization_all);
  }

  // this loop is made separate for future parallelisation of upper loop.
  // the ion contributions must be added in this exact order
  double C_total = 0.;
  for (int allionindex = 0; allionindex < get_includedions(); allionindex++) {
    C_total += grid::modelgrid[modelgridindex].ion_cooling_contribs[allionindex];
  }
  grid::modelgrid[modelgridindex].totalcooling = C_total;

  // only used in the T_e solver and write_to_estimators file
  if (heatingcoolingrates != nullptr) {
    heatingcoolingrates->cooling_collisional = C_exc_all + C_ionization_all;
    heatingcoolingrates->cooling_fb = C_fb_all;
    heatingcoolingrates->cooling_ff = C_ff_all;
  }
}

void set_kpktdiffusion(const float kpktdiffusion_timescale_in, const int n_kpktdiffusion_timesteps_in) {
  kpktdiffusion_timescale = kpktdiffusion_timescale_in;
  n_kpktdiffusion_timesteps = n_kpktdiffusion_timesteps_in;
  printout("input: kpkts diffuse %g of a time step's length for the first %d time steps\n", kpktdiffusion_timescale,
           n_kpktdiffusion_timesteps);
}

void setup_coolinglist() {
  // Determine number of processes which allow kpkts to convert to something else.
  // This number is given by the collisional excitations (so far determined from the oscillator strengths
  // by the van Regemorter formula, therefore totaluptrans), the number of free-bound emissions and collisional
  // ionisations (as long as we only deal with ionisation to the ground level this means for both of these
  // \sum_{elements,ions}get_nlevels(element,ion) and free-free which is \sum_{elements} get_nions(element)-1

  set_ncoolingterms();
  const size_t mem_usage_coolinglist = ncoolingterms * sizeof(CellCacheCoolingList);
  assert_always(ncoolingterms > 0);
  coolinglist.resize(ncoolingterms);
  printout("[info] mem_usage: coolinglist occupies %.3f MB\n", mem_usage_coolinglist / 1024. / 1024.);

  int i = 0;  // cooling list index
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels_currention = get_nlevels(element, ion);

      const int nionisinglevels = get_ionisinglevels(element, ion);

      // ff creation of rpkt
      // -------------------
      const int ioncharge = get_ionstage(element, ion) - 1;
      // printout("[debug] ioncharge %d, nncurrention %g, nne %g\n",ion,nncurrention,nne);
      if (ioncharge > 0) {
        coolinglist[i].type = CoolingType::FREEFREE;
        coolinglist[i].level = -99;
        coolinglist[i].upperlevel = -99;
        i++;
      }

      for (int level = 0; level < nlevels_currention; level++) {
        if (get_nuptrans(element, ion, level) > 0) {
          coolinglist[i].type = CoolingType::COLLEXC;
          coolinglist[i].level = level;
          // upper level is not valid because this is the contribution of all upper levels combined - have to
          // calculate individually when selecting a random process
          coolinglist[i].upperlevel = -1;
          i++;
        }
      }

      if (ion < (nions - 1))  // check whether further ionisation stage available
      {
        for (int level = 0; level < nionisinglevels; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
            coolinglist[i].type = CoolingType::COLLION;
            coolinglist[i].level = level;
            coolinglist[i].upperlevel = upper;
            i++;
          }
        }

        // fb creation of r-pkt
        // free bound rates are calculated from the lower ion, but associated to the higher ion
        for (int level = 0; level < nionisinglevels; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
            coolinglist[i].type = CoolingType::FREEBOUND;
            coolinglist[i].level = level;
            coolinglist[i].upperlevel = upper;
            i++;
          }
        }
      }
      assert_always(i == get_coolinglistoffset(element, ion) + get_ncoolingterms_ion(element, ion));
    }
  }

  assert_always(ncoolingterms == i);  // if this doesn't match, we miscalculated the number of cooling terms
  printout("[info] read_atomicdata: number of coolingterms %d\n", ncoolingterms);
}

__host__ __device__ void do_kpkt_blackbody(Packet &pkt)
// handle a k-packet (e.g., in a thick cell) by emitting according to the planck function
{
  const int modelgridindex = grid::get_cell_modelgridindex(pkt.where);

  if (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY >= 0. && grid::modelgrid[modelgridindex].thick != 1) {
    const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
    pkt.nu_cmf = sample_planck_times_expansion_opacity(nonemptymgi);
  } else {
    pkt.nu_cmf = sample_planck_montecarlo(grid::get_Te(modelgridindex));
    // TODO: is this alternative method faster or more accurate or neither?
    // pkt.nu_cmf = sample_planck_analytic(grid::get_Te(modelgridindex));
  }

  assert_always(std::isfinite(pkt.nu_cmf));
  // and then emit the packet randomly in the comoving frame
  emit_rpkt(pkt);
  // printout("[debug] calculate_chi_rpkt after kpkt to rpkt by ff\n");
  pkt.next_trans = -1;  // FLAG: transition history here not important, cont. process
  // if (tid == 0) k_stat_to_r_bb++;
  stats::increment(stats::COUNTER_K_STAT_TO_R_BB);
  stats::increment(stats::COUNTER_INTERACTIONS);
  pkt.last_event = LASTEVENT_KPKT_TO_RPKT_FFBB;
  pkt.emissiontype = EMTYPE_FREEFREE;
  pkt.em_pos = pkt.pos;
  pkt.em_time = pkt.prop_time;
  pkt.nscatterings = 0;
}

// handle a k-packet (kinetic energy of the free electrons)
__host__ __device__ void do_kpkt(Packet &pkt, const double t2, const int nts) {
  const double t1 = pkt.prop_time;
  const int modelgridindex = grid::get_cell_modelgridindex(pkt.where);

  // don't calculate cooling rates after each cell crossings any longer
  // but only if we really get a kpkt and they hadn't been calculated already

  // printout("[debug] do_kpkt: propagate k-pkt\n");

  const double deltat =
      (nts < n_kpktdiffusion_timesteps) ? kpktdiffusion_timescale * globals::timesteps[nts].width : 0.;

  const double t_current = t1 + deltat;

  if (t_current > t2) {
    pkt.pos = vec_scale(pkt.pos, t2 / t1);
    pkt.prop_time = t2;
    return;
  }
  stats::increment(stats::COUNTER_INTERACTIONS);

  pkt.pos = vec_scale(pkt.pos, t_current / t1);
  pkt.prop_time = t_current;

  assert_always(grid::modelgrid[modelgridindex].totalcooling > 0.);
  const double rndcool_ion = rng_uniform() * grid::modelgrid[modelgridindex].totalcooling;

  // Randomly select the occurring cooling process
  double coolingsum = 0.;
  int element = -1;
  int ion = -1;
  for (element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (ion = 0; ion < nions; ion++) {
      const int uniqueionindex = get_uniqueionindex(element, ion);
      coolingsum += grid::modelgrid[modelgridindex].ion_cooling_contribs[uniqueionindex];
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

  if (element >= get_nelements() || element < 0 || ion >= get_nions(element) || ion < 0) {
    printout("do_kpkt: problem selecting a cooling process ... abort\n");
    printout("do_kpkt: modelgridindex %d element %d ion %d\n", modelgridindex, element, ion);
    printout("do_kpkt: totalcooling %g, coolingsum %g, rndcool_ion %g\n", grid::modelgrid[modelgridindex].totalcooling,
             coolingsum, rndcool_ion);
    printout("do_kpkt: modelgridindex %d, cellno %d, nne %g\n", modelgridindex, pkt.where,
             grid::get_nne(modelgridindex));
    for (element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (ion = 0; ion < nions; ion++) {
        const int uniqueionindex = get_uniqueionindex(element, ion);
        printout("do_kpkt: element %d, ion %d, coolingcontr %g\n", element, ion,
                 grid::modelgrid[modelgridindex].ion_cooling_contribs[uniqueionindex]);
      }
    }
    std::abort();
  }

  const int ilow = get_coolinglistoffset(element, ion);
  const int ihigh = ilow + get_ncoolingterms_ion(element, ion) - 1;
  double C_ion_procsum = globals::cellcache[cellcacheslotid].cooling_contrib[ihigh];

  if (C_ion_procsum < 0.) {
    // printout("calculate kpkt rates on demand modelgridindex %d element %d ion %d ilow %d ihigh %d
    // oldcoolingsum %g\n",
    //          modelgridindex, element, ion, ilow, high, oldcoolingsum);
    C_ion_procsum = calculate_cooling_rates_ion<true>(modelgridindex, element, ion, ilow, cellcacheslotid, nullptr,
                                                      nullptr, nullptr, nullptr);
    assert_testmodeonly(
        (std::fabs(C_ion_procsum -
                   grid::modelgrid[modelgridindex].ion_cooling_contribs[get_uniqueionindex(element, ion)]) /
         C_ion_procsum) < 1e-3);
  }

  // with the ion selected, we now select a level and transition type

  const double rndcool_ion_process = rng_uniform() * C_ion_procsum;

  const auto i =
      std::upper_bound(globals::cellcache[cellcacheslotid].cooling_contrib + ilow,
                       globals::cellcache[cellcacheslotid].cooling_contrib + ihigh + 1, rndcool_ion_process) -
      globals::cellcache[cellcacheslotid].cooling_contrib;

  if (i > ihigh) {
    printout("do_kpkt: error occurred while selecting a cooling channel: low %d, high %d, i %td, rndcool %g\n", ilow,
             ihigh, i, rndcool_ion_process);
    printout("element %d, ion %d, offset %d, terms %d, coolingsum %g\n", element, ion,
             get_coolinglistoffset(element, ion), get_ncoolingterms_ion(element, ion), coolingsum);

    printout("lower %g, %g, %g\n",
             globals::cellcache[cellcacheslotid].cooling_contrib[get_coolinglistoffset(element, ion) - 1],
             globals::cellcache[cellcacheslotid].cooling_contrib[get_coolinglistoffset(element, ion)],
             globals::cellcache[cellcacheslotid].cooling_contrib[get_coolinglistoffset(element, ion) + 1]);
    const int finalpos = get_coolinglistoffset(element, ion) + get_ncoolingterms_ion(element, ion) - 1;
    printout("upper %g, %g, %g\n", globals::cellcache[cellcacheslotid].cooling_contrib[finalpos - 1],
             globals::cellcache[cellcacheslotid].cooling_contrib[finalpos],
             globals::cellcache[cellcacheslotid].cooling_contrib[finalpos + 1]);
  }

  assert_always(i <= ihigh);

  // printout("do_kpkt: selected process %d, coolingsum %g\n", i, coolingsum);
  const auto rndcoolingtype = coolinglist[i].type;
  const auto T_e = grid::get_Te(modelgridindex);

  if (rndcoolingtype == CoolingType::FREEFREE) {
    // The k-packet converts directly into a r-packet by free-free-emission.
    // Need to select the r-packets frequency and a random direction in the
    // co-moving frame.
    // printout("[debug] do_kpkt: k-pkt -> free-free\n");

    // Sample the packets comoving frame frequency according to paperII 5.4.3 eq.41

    const double zrand = rng_uniform_pos();  // delivers zrand in ]0,1[
    pkt.nu_cmf = -KB * T_e / H * log(zrand);

    assert_always(std::isfinite(pkt.nu_cmf));

    // and then emit the packet randomly in the comoving frame
    emit_rpkt(pkt);
    pkt.next_trans = -1;  // FLAG: transition history here not important, cont. process
    stats::increment(stats::COUNTER_K_STAT_TO_R_FF);

    pkt.last_event = LASTEVENT_KPKT_TO_RPKT_FFBB;
    pkt.emissiontype = EMTYPE_FREEFREE;
    pkt.em_pos = pkt.pos;
    pkt.em_time = pkt.prop_time;
    pkt.nscatterings = 0;
    if constexpr (VPKT_ON) {
      vpkt_call_estimators(pkt, TYPE_KPKT);
    }

  } else if (rndcoolingtype == CoolingType::FREEBOUND) {
    // The k-packet converts directly into a r-packet by free-bound-emission.
    // Need to select the r-packets frequency and a random direction in the
    // co-moving frame.
    const int lowerion = ion;
    const int lowerlevel = coolinglist[i].level;
    const int upper = coolinglist[i].upperlevel;

    // then randomly sample the packets frequency according to the continuums
    // energy distribution

    // Sample the packets comoving frame frequency according to paperII 4.2.2
    // const double zrand = rng_uniform();
    // if (zrand < 0.5) {
    pkt.nu_cmf = select_continuum_nu(element, lowerion, lowerlevel, upper, T_e);
    // } else {
    //   // Emit like a BB
    //   pkt.nu_cmf = sample_planck(T_e);
    // }

    // and then emitt the packet randomly in the comoving frame
    emit_rpkt(pkt);

    if constexpr (TRACK_ION_STATS) {
      stats::increment_ion_stats(modelgridindex, element, lowerion + 1, stats::ION_RADRECOMB_KPKT,
                                 pkt.e_cmf / H / pkt.nu_cmf);
    }

    pkt.next_trans = -1;  // FLAG: transition history here not important, cont. process
    stats::increment(stats::COUNTER_K_STAT_TO_R_FB);
    pkt.last_event = LASTEVENT_KPKT_TO_RPKT_FB;
    pkt.emissiontype = get_emtype_continuum(element, lowerion, lowerlevel, upper);
    pkt.trueemissiontype = pkt.emissiontype;
    pkt.em_pos = pkt.pos;
    pkt.em_time = pkt.prop_time;
    pkt.nscatterings = 0;

    if constexpr (VPKT_ON) {
      vpkt_call_estimators(pkt, TYPE_KPKT);
    }
  } else if (rndcoolingtype == CoolingType::COLLEXC) {
    // the k-packet activates a macro-atom due to collisional excitation
    // printout("[debug] do_kpkt: k-pkt -> collisional excitation of MA\n");
    const float nne = grid::get_nne(modelgridindex);

    // if the previous entry belongs to the same ion, then pick up the cumulative sum from
    // the previous entry, otherwise start from zero
    const double contrib_low = (i > ilow) ? globals::cellcache[cellcacheslotid].cooling_contrib[i - 1] : 0.;

    double contrib = contrib_low;
    const int level = coolinglist[i].level;
    const double epsilon_current = epsilon(element, ion, level);
    const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
    const double statweight = stat_weight(element, ion, level);
    int upper = -1;
    // excitation to same ionization stage
    const int nuptrans = get_nuptrans(element, ion, level);
    const auto *const uptranslist = get_uptranslist(element, ion, level);
    for (int ii = 0; ii < nuptrans; ii++) {
      const int tmpupper = uptranslist[ii].targetlevelindex;
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

    assert_always(upper >= 0);

    if constexpr (TRACK_ION_STATS) {
      stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_MACROATOM_ENERGYIN_COLLEXC, pkt.e_cmf);
    }

    pkt.type = TYPE_MA;
    stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_COLLEXC);
    stats::increment(stats::COUNTER_K_STAT_TO_MA_COLLEXC);

    pkt.last_event = 8;
    pkt.trueemissiontype = EMTYPE_NOTSET;
    pkt.trueemissionvelocity = -1;

    do_macroatom(pkt, {.element = element, .ion = ion, .level = upper, .activatingline = -99});
  } else if (rndcoolingtype == CoolingType::COLLION) {
    // the k-packet activates a macro-atom due to collisional ionisation
    // printout("[debug] do_kpkt: k-pkt -> collisional ionisation of MA\n");

    const int upperion = ion + 1;
    const int upper = coolinglist[i].upperlevel;

    if constexpr (TRACK_ION_STATS) {
      stats::increment_ion_stats(modelgridindex, element, upperion, stats::ION_MACROATOM_ENERGYIN_COLLION, pkt.e_cmf);
    }

    pkt.type = TYPE_MA;
    stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_COLLION);
    stats::increment(stats::COUNTER_K_STAT_TO_MA_COLLION);

    pkt.last_event = 9;
    pkt.trueemissiontype = EMTYPE_NOTSET;
    pkt.trueemissionvelocity = -1;

    do_macroatom(pkt, {.element = element, .ion = upperion, .level = upper, .activatingline = -99});
  } else if constexpr (TESTMODE) {
    assert_always(false);
  } else {
    __builtin_unreachable();
  }
}

}  // namespace kpkt
