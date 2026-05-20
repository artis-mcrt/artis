#include "kpkt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <span>
#include <utility>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "macroatom.h"
#include "mpi_logging.h"
#include "packet.h"
#include "radfield.h"
#include "random.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "stats.h"
#include "thermalbalance.h"
#include "vectors.h"
#include "vpkt.h"

namespace kpkt {

namespace {

enum class CoolingType : std::uint8_t { FREEFREE, FREEBOUND, COLLEXC, COLLION };

MPI_shared_array<const CoolingType> coolinglist_type;
MPI_shared_array<const int> coolinglist_level;
MPI_shared_array<const int> coolinglist_upperlevel;

// MK: To reduce the work imbalance between different MPI tasks I introduced a diffusion
// for kpkts, since it turned out that this work imbalance was largely dominated
// by continuous collisional interactions. By introducing a diffusion time for kpkts
// this loop is broken.
// kpktdiffusion_timescale gives the relative fraction of a time step which individual
// kpkts live.
constexpr float kpktdiffusion_timescale{0.001};

// calculate the cooling contribution list of individual levels/processes for an ion
// oldcoolingsum is the sum of lower ion (of same element or all ions of lower elements) cooling contributions
template <bool update_cellcache_contribs>
auto calculate_cooling_rates_ion(const int nonemptymgi, const int element, const int ion,
                                 std::span<double> ion_contribs, double* const C_ff, double* const C_fb,
                                 double* const C_exc, double* const C_ionisation) -> double {
  const auto nne = grid::get_nne(nonemptymgi);
  const auto T_e = grid::get_Te(nonemptymgi);

  double C_ion = 0.;
  [[maybe_unused]] int i = 0;  // NOLINT(misc-const-correctness)

  const int nionisinglevels = get_nlevels_ionising(element, ion);
  const double nncurrention = get_nnion(nonemptymgi, element, ion);

  // ff creation of rpkt
  const int ioncharge = get_ionstage(element, ion) - 1;
  if (ioncharge > 0) {
    const double C_ff_ion = 1.426e-27 * sqrt(T_e) * pow2(ioncharge) * nncurrention * nne;
    C_ion += C_ff_ion;

    if constexpr (update_cellcache_contribs) {
      ion_contribs[i] = C_ion;

      assert_testmodeonly(coolinglist_type[get_coolinglistoffset(element, ion) + i] == CoolingType::FREEFREE);

      i++;
    } else {
      *C_ff += C_ff_ion;
    }
  }

  const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
  // excitation to same ionisation stage
  const int nlevels = get_nlevels(element, ion);
  for (int level = 0; level < nlevels; level++) {
    const auto uniquelevelindex = ionuniquelevelindexstart + level;
    const double nnlevel = update_cellcache_contribs ? get_cellcache_levelpop(nonemptymgi, uniquelevelindex)
                                                     : calculate_levelpop(nonemptymgi, element, ion, level);

    const double epsilon_current = epsilon(uniquelevelindex);
    const double statweight = stat_weight(uniquelevelindex);

    const auto alltrans_startup = get_alltrans_startup(uniquelevelindex);
    const int nuptrans = get_nuptrans(uniquelevelindex);
    for (int alltransindex = alltrans_startup; alltransindex < (alltrans_startup + nuptrans); alltransindex++) {
      const int upper = globals::alltrans.targetlevelindex[alltransindex];
      const double epsilon_trans = epsilon(ionuniquelevelindexstart + upper) - epsilon_current;
      const auto upper_statweight = stat_weight(ionuniquelevelindexstart + upper);
      const double C = nnlevel *
                       col_excitation_ratecoeff(T_e, nne, upper_statweight, alltransindex, epsilon_trans, statweight) *
                       epsilon_trans;
      C_ion += C;
      if constexpr (!update_cellcache_contribs) {
        *C_exc += C;
      }
    }
    if constexpr (update_cellcache_contribs) {
      if (nuptrans > 0) {
        ion_contribs[i] = C_ion;

        assert_testmodeonly(coolinglist_type[get_coolinglistoffset(element, ion) + i] == CoolingType::COLLEXC);

        i++;
      }
    }
  }

  if (ion < (get_nions(element) - 1) && globals::nbfcontinua > 0) {
    const double nnupperion = get_nnion(nonemptymgi, element, ion + 1);

    // ionisation to higher ionisation stage
    for (int level = 0; level < nionisinglevels; level++) {
      const auto uniquelevelindex = ionuniquelevelindexstart + level;
      const double epsilon_current = epsilon(uniquelevelindex);
      const double nnlevel = update_cellcache_contribs ? get_cellcache_levelpop(nonemptymgi, uniquelevelindex)
                                                       : calculate_levelpop(nonemptymgi, element, ion, level);
      const int nphixstargets = get_nphixstargets(uniquelevelindex);
      for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
        const int upper = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
        const double epsilon_upper = epsilon(element, ion + 1, upper);
        const double epsilon_trans = epsilon_upper - epsilon_current;
        const double C = nnlevel *
                         col_ionisation_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans) *
                         epsilon_trans;

        C_ion += C;
        if constexpr (update_cellcache_contribs) {
          ion_contribs[i] = C_ion;

          assert_testmodeonly(coolinglist_type[get_coolinglistoffset(element, ion) + i] == CoolingType::COLLION);
          assert_testmodeonly(coolinglist_level[get_coolinglistoffset(element, ion) + i] == level);
          assert_testmodeonly(coolinglist_upperlevel[get_coolinglistoffset(element, ion) + i] == upper);

          i++;
        } else {
          *C_ionisation += C;
        }
      }
    }

    // fb creation of r-pkt
    for (int level = 0; level < nionisinglevels; level++) {
      const auto uniquelevelindex = ionuniquelevelindexstart + level;
      const int nphixstargets = get_nphixstargets(uniquelevelindex);
      for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
        const double pop = [&]() {
          if constexpr (BFCOOLING_USELEVELPOPNOTIONPOP) {
            const int upperlevel = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
            return update_cellcache_contribs ? get_cellcache_levelpop(nonemptymgi, element, ion + 1, upperlevel)
                                             : calculate_levelpop(nonemptymgi, element, ion + 1, upperlevel);
          }
          return nnupperion;
        }();
        const double C = get_bfcoolingcoeff(element, ion, level, phixstargetindex, T_e) * pop * nne;
        C_ion += C;

        if constexpr (update_cellcache_contribs) {
          ion_contribs[i] = C_ion;

          assert_testmodeonly(coolinglist_type[get_coolinglistoffset(element, ion) + i] == CoolingType::FREEBOUND);
          assert_testmodeonly(coolinglist_level[get_coolinglistoffset(element, ion) + i] == level);
          assert_testmodeonly(coolinglist_upperlevel[get_coolinglistoffset(element, ion) + i] ==
                              get_phixsupperlevel(uniquelevelindex, phixstargetindex));

          i++;
        } else {
          *C_fb += C;
        }
      }
    }
  }

  if constexpr (update_cellcache_contribs) {
    assert_always(i == std::ssize(ion_contribs));
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

// return a randomly chosen frequency according to the Planck distribution of temperature T using a Monte Carlo method
auto sample_planck_montecarlo(const double T) -> double {
  const double nu_peak = 5.879e10 * T;
  const double B_peak = radfield::planck(nu_peak, T);

  while (true) {
    const double nu = NU_MIN_R + (rng_uniform() * (NU_MAX_R - NU_MIN_R));
    if (rng_uniform() * B_peak <= radfield::planck(nu, T)) {
      return nu;
    }
  }
}
}  // anonymous namespace

// Calculate the cooling rates for a given cell and store them for each ion
// optionally store components (ff, bf, collisional) in heatingcoolingrates struct
void calculate_cooling_rates(const int nonemptymgi, HeatingCoolingRates* heatingcoolingrates) {
  double C_ff_all = 0.;  // free-free creation of rpkts
  double C_fb_all = 0.;  // free-bound creation of rpkt
  double C_exc_all = 0.;  // collisional excitation of macroatoms
  double C_ionisation_all = 0.;  // collisional ionisation of macroatoms
  const auto nincludedions = get_includedions();
  const auto cellioncontribs = get_cell_ion_cooling_contribs(nonemptymgi);
  double cumulative_cooling = 0.;
  for (int uniqueionindex = 0; uniqueionindex < nincludedions; uniqueionindex++) {
    const auto [element, ion] = get_ionfromuniqueionindex(uniqueionindex);
    cumulative_cooling += calculate_cooling_rates_ion<false>(nonemptymgi, element, ion, {}, &C_ff_all, &C_fb_all,
                                                             &C_exc_all, &C_ionisation_all);
    cellioncontribs[uniqueionindex] = cumulative_cooling;
  }

  // only used in the T_e solver and write_to_estimators file
  if (heatingcoolingrates != nullptr) {
    heatingcoolingrates->cooling_collisional = C_exc_all + C_ionisation_all;
    heatingcoolingrates->cooling_fb = C_fb_all;
    heatingcoolingrates->cooling_ff = C_ff_all;
  }
}

void setup_coolinglist() {
  // Determine number of processes which allow kpkts to convert to something else.
  // This number is given by the collisional excitations (so far determined from the oscillator strengths
  // by the van Regemorter formula, therefore totaluptrans), the number of free-bound emissions and collisional
  // ionisations (as long as we only deal with ionisation to the ground level this means for both of these
  // \sum_{elements,ions}get_nlevels(element,ion) and free-free which is \sum_{elements} get_nions(element)-1

  set_ncoolingterms();
  assert_always(ncoolingterms > 0);
  auto temp_coolinglist_type = MPI_shared_array<CoolingType>(ncoolingterms);
  auto temp_coolinglist_level = MPI_shared_array<int>(ncoolingterms);
  auto temp_coolinglist_upperlevel = MPI_shared_array<int>(ncoolingterms);
  const size_t mem_usage_coolinglist = ncoolingterms * (sizeof(CoolingType) + (2 * sizeof(int)));
  printlnlog("[info] mem_usage: coolinglist occupies {:.3f} MB", mem_usage_coolinglist / 1024. / 1024.);

  int i = 0;  // cooling list index
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels_currention = get_nlevels(element, ion);
      const int nionisinglevels = get_nlevels_ionising(element, ion);

      // ff creation of rpkt
      const int ioncharge = get_ionstage(element, ion) - 1;
      if (ioncharge > 0) {
        temp_coolinglist_type[i] = CoolingType::FREEFREE;
        temp_coolinglist_level[i] = -99;
        temp_coolinglist_upperlevel[i] = -99;
        i++;
      }

      for (int level = 0; level < nlevels_currention; level++) {
        if (get_nuptrans(element, ion, level) > 0) {
          temp_coolinglist_type[i] = CoolingType::COLLEXC;
          temp_coolinglist_level[i] = level;
          // upper level is not valid because this is the contribution of all upper levels combined - have to
          // calculate individually when selecting a random process
          temp_coolinglist_upperlevel[i] = -1;
          i++;
        }
      }

      // check whether further ionisation stage available
      if (ion < (nions - 1)) {
        // free-bound creation of r-pkt

        // collisional ionisation to higher ionisation stage
        for (int level = 0; level < nionisinglevels; level++) {
          const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
          const int nphixstargets = get_nphixstargets(uniquelevelindex);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            const int upper = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
            temp_coolinglist_type[i] = CoolingType::COLLION;
            temp_coolinglist_level[i] = level;
            temp_coolinglist_upperlevel[i] = upper;
            i++;
          }
        }

        // free-bound creation of r-pkt from recombination
        // this should probably be considered cooling of the higher ion, but for simplicity we store in the list of the
        // lower ion
        for (int level = 0; level < nionisinglevels; level++) {
          const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
          const int nphixstargets = get_nphixstargets(uniquelevelindex);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            const int upper = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
            temp_coolinglist_type[i] = CoolingType::FREEBOUND;
            temp_coolinglist_level[i] = level;
            temp_coolinglist_upperlevel[i] = upper;
            i++;
          }
        }
      }
      assert_always(i == get_coolinglistoffset(element, ion) + get_ncoolingterms_ion(element, ion));
    }
  }

  assert_always(ncoolingterms == i);  // if this doesn't match, we miscalculated the number of cooling terms
  printlnlog("[info] read_atomicdata: number of coolingterms {}", ncoolingterms);
  coolinglist_type = std::move(temp_coolinglist_type);
  coolinglist_level = std::move(temp_coolinglist_level);
  coolinglist_upperlevel = std::move(temp_coolinglist_upperlevel);
  MPI_Barrier_node();

  printlnlog("kpkts diffuse {:g} of a time step's length", kpktdiffusion_timescale);
}

// handle a k-packet (e.g., in a thick cell) by emitting according to the planck function
DEVICE_FUNC void do_kpkt_blackbody(Packet& pkt) {
  const auto nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);

  if (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value() && grid::thick_allcells[nonemptymgi] != 1) {
    pkt.nu_cmf = sample_planck_times_expansion_opacity(nonemptymgi);
  } else {
    pkt.nu_cmf = sample_planck_montecarlo(grid::get_Te(nonemptymgi));
  }

  assert_always(std::isfinite(pkt.nu_cmf));
  // and then emit the packet randomly in the comoving frame
  emit_rpkt(pkt);
  pkt.next_trans = -1;  // FLAG: transition history here not important, cont. process
  stats::increment(stats::Counter::K_STAT_TO_R_BB);
  stats::increment(stats::Counter::INTERACTIONS);
  pkt.emissiontype = EMTYPE_FREEFREE;
  pkt.em_pos = pkt.pos;
  pkt.em_time = static_cast<float>(pkt.prop_time);
  pkt.nscatterings = 0;
}

// handle a k-packet (kinetic energy of the free electrons)
DEVICE_FUNC void do_kpkt(Packet& pkt, const double t2, const int nts) {
  const double deltat = kpktdiffusion_timescale * globals::timesteps[nts].width;

  const double t_current = std::min(pkt.prop_time + deltat, t2);

  pkt.pos = vec_scale(pkt.pos, t_current / pkt.prop_time);
  // pkt.e_cmf *= pkt.prop_time / t_current;  // adjust energy for adiabatic losses
  pkt.prop_time = t_current;

  if (t_current >= t2) {
    return;
  }

  stats::increment(stats::Counter::INTERACTIONS);

  const auto nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);
  const std::span<const double> ion_cooling_contribs_thiscell = get_cell_ion_cooling_contribs(nonemptymgi);
  const double rndcool_ion = rng_uniform() * ion_cooling_contribs_thiscell.back();

  // Randomly select the occurring cooling process
  const int uniqueionindex = static_cast<int>(std::ranges::lower_bound(ion_cooling_contribs_thiscell, rndcool_ion) -
                                              ion_cooling_contribs_thiscell.begin());
  assert_always(uniqueionindex < get_includedions());
  const auto [element, ion] = get_ionfromuniqueionindex(uniqueionindex);

  const int ionstart = get_coolinglistoffset(element, ion);
  const int ncoolingterms_ion = get_ncoolingterms_ion(element, ion);
  assert_testmodeonly(ncoolingterms_ion > 0);
  // transition contribution list for this ion
  const std::span<double> ion_contribs =
      std::span{globals::cellcache[cellcacheslotid].cooling_contrib}.subspan(ionstart, ncoolingterms_ion);

  {
#if (defined(STDPAR_ON) || defined(_OPENMP)) && !defined(GPU_ON)
    [[maybe_unused]] ScopedMutex lock{globals::cellcache[cellcacheslotid].cooling_contrib_locks[uniqueionindex]};
#endif
    if (ion_contribs[0] < 0.) {
      calculate_cooling_rates_ion<true>(nonemptymgi, element, ion, ion_contribs, nullptr, nullptr, nullptr, nullptr);
    }
  }

  // subspan for this ion's region of the cumulative sum of cooling contributions
  const double C_ion_procsum = ion_contribs.back();
  assert_testmodeonly(C_ion_procsum > 0.);

  // with the ion selected, we now select a level and transition type

  const double rndcool_ion_process = rng_uniform() * C_ion_procsum;

  const auto ionoffset = std::ranges::upper_bound(ion_contribs, rndcool_ion_process) - ion_contribs.begin();
  assert_always(ionoffset < ncoolingterms_ion);
  const auto i = ionstart + ionoffset;

  const auto rndcoolingtype = coolinglist_type[i];
  const auto T_e = grid::get_Te(nonemptymgi);

  if (rndcoolingtype == CoolingType::FREEFREE) {
    // The k-packet converts directly into a r-packet by free-free emission.
    // Need to select the r-packets frequency and a random direction in the
    // co-moving frame.

    // Sample the packets comoving frame frequency according to paperII 5.4.3 eq.41

    pkt.nu_cmf = -KB * T_e / H * std::log(static_cast<double>(rng_uniform_pos()));

    assert_always(std::isfinite(pkt.nu_cmf));

    // and then emit the packet randomly in the comoving frame
    emit_rpkt(pkt);
    pkt.next_trans = -1;  // FLAG: transition history here not important, cont. process
    stats::increment(stats::Counter::K_STAT_TO_R_FF);

    pkt.emissiontype = EMTYPE_FREEFREE;
    pkt.em_pos = pkt.pos;
    pkt.em_time = static_cast<float>(pkt.prop_time);
    pkt.nscatterings = 0;
    if constexpr (VPKT_ON) {
      vpkt::trace_vpkts(pkt, TYPE_KPKT);
    }

  } else if (rndcoolingtype == CoolingType::FREEBOUND) {
    // The k-packet converts directly into a r-packet by free-bound-emission.
    // Need to select the r-packets frequency and a random direction in the co-moving frame.
    const int lowerion = ion;
    const int lowerlevel = coolinglist_level[i];
    const int upper = coolinglist_upperlevel[i];

    // then randomly sample the packets frequency according to the continuums energy distribution

    // Sample the packets comoving frame frequency according to paperII 4.2.2
    pkt.nu_cmf = select_continuum_nu(element, lowerion, lowerlevel, upper, T_e);

    // and then emit the packet randomly in the comoving frame
    emit_rpkt(pkt);

    pkt.next_trans = -1;  // FLAG: transition history here not important, cont. process
    stats::increment(stats::Counter::K_STAT_TO_R_FB);
    pkt.emissiontype = get_emtype_continuum(element, lowerion, lowerlevel, upper);
    pkt.trueemissiontype = pkt.emissiontype;
    pkt.em_pos = pkt.pos;
    pkt.em_time = static_cast<float>(pkt.prop_time);
    pkt.nscatterings = 0;

    if constexpr (VPKT_ON) {
      vpkt::trace_vpkts(pkt, TYPE_KPKT);
    }
  } else if (rndcoolingtype == CoolingType::COLLEXC) {
    // the k-packet activates a macro-atom due to collisional excitation
    const float nne = grid::get_nne(nonemptymgi);

    // if the previous entry belongs to the same ion, then pick up the cumulative sum from
    // the previous entry, otherwise start from zero
    const double contrib_low = (i > ionstart) ? globals::cellcache[cellcacheslotid].cooling_contrib[i - 1] : 0.;

    double contrib = contrib_low;
    const int level = coolinglist_level[i];
    const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
    const auto uniquelevelindex = ionuniquelevelindexstart + level;
    const double epsilon_current = epsilon(uniquelevelindex);
    const double nnlevel = get_cellcache_levelpop(nonemptymgi, uniquelevelindex);
    const double statweight = stat_weight(uniquelevelindex);
    int upper = -1;
    // excitation to same ionisation stage
    const auto alltrans_startup = get_alltrans_startup(uniquelevelindex);
    const int nuptrans = get_nuptrans(uniquelevelindex);
    for (int alltransindex = alltrans_startup; alltransindex < (alltrans_startup + nuptrans); alltransindex++) {
      const int tmpupper = globals::alltrans.targetlevelindex[alltransindex];
      const auto upperuniquelevelindex = ionuniquelevelindexstart + tmpupper;
      const double epsilon_trans = epsilon(upperuniquelevelindex) - epsilon_current;
      const auto upper_statweight = stat_weight(upperuniquelevelindex);
      const double C = nnlevel *
                       col_excitation_ratecoeff(T_e, nne, upper_statweight, alltransindex, epsilon_trans, statweight) *
                       epsilon_trans;
      contrib += C;
      if (contrib > rndcool_ion_process) {
        upper = tmpupper;
        break;
      }
    }

    assert_always(upper >= 0);

    stats::increment(stats::Counter::MA_STAT_ACTIVATION_COLLEXC);
    stats::increment(stats::Counter::K_STAT_TO_MA_COLLEXC);

    pkt.trueemissiontype = EMTYPE_NOTSET;
    pkt.trueem_pos = {NAN, NAN, NAN};

    do_macroatom(pkt, {.element = element, .ion = ion, .level = upper, .activatingline = -99});
  } else if (rndcoolingtype == CoolingType::COLLION) {
    // the k-packet activates a macro-atom due to collisional ionisation

    const int upperion = ion + 1;
    const int upper = coolinglist_upperlevel[i];

    stats::increment(stats::Counter::MA_STAT_ACTIVATION_COLLION);
    stats::increment(stats::Counter::K_STAT_TO_MA_COLLION);

    pkt.trueemissiontype = EMTYPE_NOTSET;
    pkt.trueem_pos = {NAN, NAN, NAN};

    do_macroatom(pkt, {.element = element, .ion = upperion, .level = upper, .activatingline = -99});
  } else if constexpr (TESTMODE) {
    assert_always(false);
  } else {
    __builtin_unreachable();
  }
}

}  // namespace kpkt
