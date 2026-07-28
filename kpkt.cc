// k-packets: packets of thermal kinetic energy. Computes the cooling rates of all processes
// (free-free, free-bound, collisional excitation and ionisation) and samples a cooling channel
// to convert a k-packet into radiation or a macro-atom excitation.

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

// Fraction of a time step that individual k-packets live before further processing. This
// diffusion time breaks up the chains of continuous collisional interactions that would
// otherwise dominate the work imbalance between MPI ranks.
constexpr float kpktdiffusion_timestep_fraction{0.001};

// Compute the collisional cooling rate of a single ion, summing the free-free, free-bound, collisional-excitation
// and collisional-ionisation contributions. Accumulates the per-process totals (C_ff/C_fb/C_exc/C_ionisation),
// records each individual term in ion_contribs, and returns the ion's total cooling rate.
template <bool update_cellcache_contribs>
auto calculate_cooling_rates_ion(const int nonemptymgi, const int element, const int ion,
                                 std::span<double> ion_contribs, double* const C_ff, double* const C_fb,
                                 double* const C_exc, double* const C_ionisation) -> double {
  const auto clumpednne = grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi);
  const auto T_e = grid::Te_allcells[nonemptymgi];

  double C_ion = 0.;
  // cursor into this ion's slice of the cooling list, only advanced when we are filling the cellcache
  [[maybe_unused]] int coolinglistindex = 0;  // NOLINT(misc-const-correctness)

  const int nionisinglevels = get_nlevels_ionising(element, ion);
  const double nncurrention = get_nnion(nonemptymgi, element, ion);

  // ff creation of rpkt
  const int ioncharge = get_ionstage(element, ion) - 1;
  if (ioncharge > 0) {
    const double C_ff_ion = 1.426e-27 * sqrt(T_e) * pow2(ioncharge) * nncurrention * clumpednne;
    C_ion += C_ff_ion;

    if constexpr (update_cellcache_contribs) {
      ion_contribs[coolinglistindex] = C_ion;

      assert_testmodeonly(coolinglist_type[get_coolinglistoffset(element, ion) + coolinglistindex] ==
                          CoolingType::FREEFREE);

      coolinglistindex++;
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
      const double C =
          nnlevel *
          col_excitation_ratecoeff(T_e, clumpednne, upper_statweight, alltransindex, epsilon_trans, statweight) *
          epsilon_trans;
      C_ion += C;
      if constexpr (!update_cellcache_contribs) {
        *C_exc += C;
      }
    }
    if constexpr (update_cellcache_contribs) {
      if (nuptrans > 0) {
        ion_contribs[coolinglistindex] = C_ion;

        assert_testmodeonly(coolinglist_type[get_coolinglistoffset(element, ion) + coolinglistindex] ==
                            CoolingType::COLLEXC);

        coolinglistindex++;
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
        const double C =
            nnlevel * col_ionisation_ratecoeff(T_e, clumpednne, element, ion, level, phixstargetindex, epsilon_trans) *
            epsilon_trans;

        C_ion += C;
        if constexpr (update_cellcache_contribs) {
          ion_contribs[coolinglistindex] = C_ion;

          assert_testmodeonly(coolinglist_type[get_coolinglistoffset(element, ion) + coolinglistindex] ==
                              CoolingType::COLLION);
          assert_testmodeonly(coolinglist_level[get_coolinglistoffset(element, ion) + coolinglistindex] == level);
          assert_testmodeonly(coolinglist_upperlevel[get_coolinglistoffset(element, ion) + coolinglistindex] == upper);

          coolinglistindex++;
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
        const double pop = [&] {
          if constexpr (BFCOOLING_USELEVELPOPNOTIONPOP) {
            const int upperlevel = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
            return update_cellcache_contribs ? get_cellcache_levelpop(nonemptymgi, element, ion + 1, upperlevel)
                                             : calculate_levelpop(nonemptymgi, element, ion + 1, upperlevel);
          }
          return nnupperion;
        }();
        const double C = get_bfcoolingcoeff(element, ion, level, phixstargetindex, T_e) * pop * clumpednne;
        C_ion += C;

        if constexpr (update_cellcache_contribs) {
          ion_contribs[coolinglistindex] = C_ion;

          assert_testmodeonly(coolinglist_type[get_coolinglistoffset(element, ion) + coolinglistindex] ==
                              CoolingType::FREEBOUND);
          assert_testmodeonly(coolinglist_level[get_coolinglistoffset(element, ion) + coolinglistindex] == level);
          assert_testmodeonly(coolinglist_upperlevel[get_coolinglistoffset(element, ion) + coolinglistindex] ==
                              get_phixsupperlevel(uniquelevelindex, phixstargetindex));

          coolinglistindex++;
        } else {
          *C_fb += C;
        }
      }
    }
  }

  if constexpr (update_cellcache_contribs) {
    assert_always(coolinglistindex == std::ssize(ion_contribs));
  }

  return C_ion;
}

// Count the total number of cooling processes (terms) over all ions and record each ion's starting offset
// (coolingoffset) into the global cooling list.
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
auto sample_planck_montecarlo(const double T, rngstate_type& rngstate) -> double {
  const double nu_peak = 5.879e10 * T;
  const double B_peak = radfield::planck(nu_peak, T);

  while (true) {
    const double nu = NU_MIN_R + (rng_uniform(rngstate) * (NU_MAX_R - NU_MIN_R));
    if (rng_uniform(rngstate) * B_peak <= radfield::planck(nu, T)) {
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
  printlnlog("[info] setup_coolinglist: number of coolingterms {}", ncoolingterms);
  coolinglist_type = std::move(temp_coolinglist_type);
  coolinglist_level = std::move(temp_coolinglist_level);
  coolinglist_upperlevel = std::move(temp_coolinglist_upperlevel);
  MPI_Barrier_node();

  printlnlog("kpkts diffuse {:g} of a time step's length", kpktdiffusion_timestep_fraction);
}

// prepopulate one ion's cooling-rate contributions into the cellcache (see header)
DEVICE_FUNC void calculate_cellcache_cooling_rates_ion(const int nonemptymgi, const int uniqueionindex) {
  const auto [element, ion] = get_ionfromuniqueionindex(uniqueionindex);
  const int ionstart = get_coolinglistoffset(element, ion);
  const int ncoolingterms_ion = get_ncoolingterms_ion(element, ion);
  const auto ion_contribs = std::span{get_cellcache(nonemptymgi).cooling_contrib}.subspan(ionstart, ncoolingterms_ion);
  calculate_cooling_rates_ion<true>(nonemptymgi, element, ion, ion_contribs, nullptr, nullptr, nullptr, nullptr);
}

// handle a k-packet (e.g., in a thick cell) by emitting according to the planck function
DEVICE_FUNC void do_kpkt_blackbody(Packet& pkt) {
  const auto nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);

  if (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value() && grid::thick_allcells[nonemptymgi] != 1) {
    pkt.nu_cmf = sample_planck_times_expansion_opacity(nonemptymgi, get_rngstate(pkt));
  } else {
    pkt.nu_cmf = sample_planck_montecarlo(grid::Te_allcells[nonemptymgi], get_rngstate(pkt));
  }

  assert_always(std::isfinite(pkt.nu_cmf));
  // and then emit the packet randomly in the comoving frame
  emit_rpkt(pkt);
  pkt.next_trans = -1;  // FLAG: transition history here not important, cont. process
  stats::increment(stats::Counter::K_STAT_TO_R_BB);
  stats::increment(stats::Counter::INTERACTIONS);
  pkt.emissiontype = EMTYPE_FREEFREE;
  // this is a thermal emission, so record it as the packet's last thermal ("true") emission
  // (emit_rpkt has just set em_pos/em_time to the current position and time)
  pkt.trueemissiontype = pkt.emissiontype;
  pkt.trueem_pos = pkt.em_pos;
  pkt.trueem_time = pkt.em_time;
  pkt.nscatterings = 0;
}

// handle a k-packet (kinetic energy of the free electrons)
DEVICE_FUNC void do_kpkt(Packet& pkt, const double t2, const int nts) {
  const double deltat = kpktdiffusion_timestep_fraction * globals::timesteps[nts].width;

  const double t_current = std::min(pkt.prop_time + deltat, t2);

  pkt.pos = vec_scale(pkt.pos, t_current / pkt.prop_time);
  pkt.e_cmf *= pkt.prop_time / t_current;  // adjust energy for adiabatic losses
  pkt.prop_time = t_current;

  if (t_current >= t2) {
    return;
  }

  stats::increment(stats::Counter::INTERACTIONS);

  const auto nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);
  const std::span<const double> ion_cooling_contribs_thiscell = get_cell_ion_cooling_contribs(nonemptymgi);
  const double rndcool_ion = rng_uniform(get_rngstate(pkt)) * ion_cooling_contribs_thiscell.back();

  // Randomly select the occurring cooling process: first ion whose cumulative cooling exceeds rndcool_ion.
  // upper_bound (not lower_bound) so that a zero-cooling ion is never selected when rndcool_ion ties its
  // cumulative value exactly (which would give a zero total and an ill-defined process selection).
  const int uniqueionindex = static_cast<int>(std::ranges::upper_bound(ion_cooling_contribs_thiscell, rndcool_ion) -
                                              ion_cooling_contribs_thiscell.begin());
  assert_always(uniqueionindex < get_includedions());
  const auto [element, ion] = get_ionfromuniqueionindex(uniqueionindex);

  const int ionstart = get_coolinglistoffset(element, ion);
  const int ncoolingterms_ion = get_ncoolingterms_ion(element, ion);
  assert_testmodeonly(ncoolingterms_ion > 0);
  // transition contribution list for this ion
  const std::span<double> ion_contribs =
      std::span{get_cellcache(nonemptymgi).cooling_contrib}.subspan(ionstart, ncoolingterms_ion);

  {
    const auto calc_cooling_if_needed = [&] {
      if (ion_contribs[0] < 0.) {
        calculate_cooling_rates_ion<true>(nonemptymgi, element, ion, ion_contribs, nullptr, nullptr, nullptr, nullptr);
      }
    };
#if (defined(STDPAR_ON) || defined(_OPENMP)) && !defined(GPU_ON)
    // Only single-slot mode shares a cache slot between threads working different cells and therefore needs
    // the per-ion mutex. Multi-slot mode pre-calculates every cell's cooling rates in cellcacheslot_populate,
    // so no locking (and no lazy calculation) is required here.
    if constexpr (cellcache_singleslot) {
      [[maybe_unused]] ScopedMutex lock{get_cellcache(nonemptymgi).cooling_contrib_locks[uniqueionindex]};
      calc_cooling_if_needed();
    } else
#endif
    {
      calc_cooling_if_needed();
    }
  }

  // subspan for this ion's region of the cumulative sum of cooling contributions
  const double C_ion_procsum = ion_contribs.back();
  assert_testmodeonly(C_ion_procsum > 0.);

  // with the ion selected, we now select a level and transition type

  const double rndcool_ion_process = rng_uniform(get_rngstate(pkt)) * C_ion_procsum;

  const auto ionoffset = std::ranges::upper_bound(ion_contribs, rndcool_ion_process) - ion_contribs.begin();
  assert_always(ionoffset < ncoolingterms_ion);
  const auto i = ionstart + ionoffset;

  const auto rndcoolingtype = coolinglist_type[i];
  const auto T_e = grid::Te_allcells[nonemptymgi];

  if (rndcoolingtype == CoolingType::FREEFREE) {
    // The k-packet converts directly into a r-packet by free-free emission.
    // Need to select the r-packets frequency and a random direction in the
    // co-moving frame.

    // Sample the packets comoving frame frequency according to paperII 5.4.3 eq.41

    pkt.nu_cmf = -KB * T_e / H * std::log(static_cast<double>(rng_uniform_pos(get_rngstate(pkt))));

    assert_always(std::isfinite(pkt.nu_cmf));

    // and then emit the packet randomly in the comoving frame
    emit_rpkt(pkt);
    pkt.next_trans = -1;  // FLAG: transition history here not important, cont. process
    stats::increment(stats::Counter::K_STAT_TO_R_FF);

    pkt.emissiontype = EMTYPE_FREEFREE;
    // this is a thermal emission, so record it as the packet's last thermal ("true") emission
    // (emit_rpkt has just set em_pos/em_time to the current position and time)
    pkt.trueemissiontype = pkt.emissiontype;
    pkt.trueem_pos = pkt.em_pos;
    pkt.trueem_time = pkt.em_time;
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
    pkt.nu_cmf = select_continuum_nu(element, lowerion, lowerlevel, upper, T_e, get_rngstate(pkt));

    // and then emit the packet randomly in the comoving frame
    emit_rpkt(pkt);

    pkt.next_trans = -1;  // FLAG: transition history here not important, cont. process
    stats::increment(stats::Counter::K_STAT_TO_R_FB);
    pkt.emissiontype = get_emtype_continuum(element, lowerion, lowerlevel, upper);
    pkt.trueemissiontype = pkt.emissiontype;
    pkt.trueem_pos = pkt.em_pos;
    pkt.trueem_time = pkt.em_time;
    pkt.nscatterings = 0;

    if constexpr (VPKT_ON) {
      vpkt::trace_vpkts(pkt, TYPE_KPKT);
    }
  } else if (rndcoolingtype == CoolingType::COLLEXC) {
    // the k-packet activates a macro-atom due to collisional excitation
    const float clumpednne = grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi);

    // if the previous entry belongs to the same ion, then pick up the cumulative sum from
    // the previous entry, otherwise start from zero
    const double contrib_low = (i > ionstart) ? get_cellcache(nonemptymgi).cooling_contrib[i - 1] : 0.;

    double contrib = contrib_low;
    const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
    const auto uniquelevelindex = ionuniquelevelindexstart + coolinglist_level[i];
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
      const double C =
          nnlevel *
          col_excitation_ratecoeff(T_e, clumpednne, upper_statweight, alltransindex, epsilon_trans, statweight) *
          epsilon_trans;
      contrib += C;
      // Strict comparison ensures that a zero target skips leading transitions with zero contribution.
      if (contrib > rndcool_ion_process) {
        upper = tmpupper;
        break;
      }
    }

    assert_always(contrib > rndcool_ion_process);

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
