// k-packets: packets of thermal kinetic energy. Computes the cooling rates of all processes
// (free-free, free-bound, collisional excitation and ionisation) and samples a cooling channel
// to convert a k-packet into radiation or a macro-atom excitation.
//
// k-packets are the thermal-pool state of the indivisible energy packet scheme of Lucy (2002),
// doi:10.1051/0004-6361:20011756; Lucy (2003), arXiv:astro-ph/0303202, which macroatom.cc implements.

#include "kpkt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
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
#include "sn3d.h"
#include "stats.h"
#include "thermalbalance.h"
#include "vectors.h"
#include "vpkt.h"

namespace kpkt {

namespace {

enum class CoolingType : std::uint8_t { FREEFREE, FREEBOUND, COLLEXC, COLLION };

MPI_shared_array<const CoolingType> coolinglist_type;
MPI_shared_array<const int> coolinglist_level;
MPI_shared_array<const int> coolinglist_phixstargetindex;

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
          col_excitation_ratecoeff(T_e, clumpednne, epsilon_trans, upper_statweight, statweight, alltransindex) *
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
          assert_testmodeonly(coolinglist_phixstargetindex[get_coolinglistoffset(element, ion) + coolinglistindex] ==
                              phixstargetindex);

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

      // The bf cooling coefficient of a (level, target) continuum is per population of the target level. When the
      // upper ion population is used instead, a level with several targets must share it among them: giving every
      // target the whole ion population would count the term's cooling once per target (with target probabilities
      // proportional to the target stat weights, each target's coefficient equals that of the whole target term).
      // The share is the LTE fraction among this level's targets, with Boltzmann factors relative to the
      // lowest-energy target so that they cannot all underflow. Single-target levels get the whole ion population.
      [[maybe_unused]] double targetweight_sum = 0.;
      [[maybe_unused]] double E_target_min = 0.;
      if constexpr (!BFCOOLING_USELEVELPOPNOTIONPOP) {
        if (nphixstargets > 1) {
          E_target_min = std::numeric_limits<double>::max();
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            E_target_min = std::min(E_target_min,
                                    epsilon(element, ion + 1, get_phixsupperlevel(uniquelevelindex, phixstargetindex)));
          }
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            const int upperlevel = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
            targetweight_sum += stat_weight(element, ion + 1, upperlevel) *
                                exp(-(epsilon(element, ion + 1, upperlevel) - E_target_min) / KB / T_e);
          }
        }
      }

      for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
        const double pop = [&] {
          if constexpr (BFCOOLING_USELEVELPOPNOTIONPOP) {
            const int upperlevel = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
            return update_cellcache_contribs ? get_cellcache_levelpop(nonemptymgi, element, ion + 1, upperlevel)
                                             : calculate_levelpop(nonemptymgi, element, ion + 1, upperlevel);
          }
          if (nphixstargets == 1) {
            return nnupperion;
          }
          const int upperlevel = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
          const double targetweight = stat_weight(element, ion + 1, upperlevel) *
                                      exp(-(epsilon(element, ion + 1, upperlevel) - E_target_min) / KB / T_e);
          return nnupperion * targetweight / targetweight_sum;
        }();
        const double C = get_bfcoolingcoeff(element, ion, level, phixstargetindex, T_e) * pop * clumpednne;
        C_ion += C;

        if constexpr (update_cellcache_contribs) {
          ion_contribs[coolinglistindex] = C_ion;

          assert_testmodeonly(coolinglist_type[get_coolinglistoffset(element, ion) + coolinglistindex] ==
                              CoolingType::FREEBOUND);
          assert_testmodeonly(coolinglist_level[get_coolinglistoffset(element, ion) + coolinglistindex] == level);
          assert_testmodeonly(coolinglist_phixstargetindex[get_coolinglistoffset(element, ion) + coolinglistindex] ==
                              phixstargetindex);

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

// Draw a frequency from the Planck distribution at temperature T, truncated to the [NU_MIN_R, NU_MAX_R] range
// that r-packets are tracked over. Uses rejection sampling against a uniform proposal, with the envelope set by
// the peak of the Planck function, so it stays valid whatever part of the curve falls inside the range.
auto sample_planck_montecarlo(const double T, rngstate_type& rngstate) -> double {
  const double nu_peak = 5.879e10 * T;  // Wien displacement law in frequency [Hz/K]
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

// Build the list of cooling processes that a k-packet can convert through.
//
// The number of processes is given by the collisional excitations (so far determined from the oscillator
// strengths by the van Regemorter formula, therefore totaluptrans), the number of free-bound emissions and
// collisional ionisations (as long as we only deal with ionisation to the ground level this means for both
// of these \sum_{elements,ions}get_nlevels(element,ion)) and free-free, which is
// \sum_{elements} get_nions(element)-1.
void setup_coolinglist() {
  set_ncoolingterms();
  assert_always(ncoolingterms > 0);
  auto temp_coolinglist_type = MPI_shared_array<CoolingType>(ncoolingterms);
  auto temp_coolinglist_level = MPI_shared_array<int>(ncoolingterms);
  auto temp_coolinglist_phixstargetindex = MPI_shared_array<int>(ncoolingterms);
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
        temp_coolinglist_phixstargetindex[i] = -99;
        i++;
      }

      for (int level = 0; level < nlevels_currention; level++) {
        if (get_nuptrans(element, ion, level) > 0) {
          temp_coolinglist_type[i] = CoolingType::COLLEXC;
          temp_coolinglist_level[i] = level;
          // a collisional excitation is bound-bound, so there is no photoionisation target. This entry is
          // the contribution of all upper levels combined, chosen individually when the process is selected
          temp_coolinglist_phixstargetindex[i] = -1;
          i++;
        }
      }

      // check whether further ionisation stage available
      if (ion < (nions - 1)) {
        // free-bound creation of r-pkt

        // collisional ionisation to higher ionisation stage
        for (int level = 0; level < nionisinglevels; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            temp_coolinglist_type[i] = CoolingType::COLLION;
            temp_coolinglist_level[i] = level;
            temp_coolinglist_phixstargetindex[i] = phixstargetindex;
            i++;
          }
        }

        // free-bound creation of r-pkt from recombination
        // this should probably be considered cooling of the higher ion, but for simplicity we store in the list of the
        // lower ion
        for (int level = 0; level < nionisinglevels; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            temp_coolinglist_type[i] = CoolingType::FREEBOUND;
            temp_coolinglist_level[i] = level;
            temp_coolinglist_phixstargetindex[i] = phixstargetindex;
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
  coolinglist_phixstargetindex = std::move(temp_coolinglist_phixstargetindex);
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

// Store and return the share of the thermal energy flux of a cell that becomes radiation (see
// NLTE_TIME_DEPENDENT_FIRST_TIMESTEP). The k-packets carry the heating rate. The thermal balance spends a part
// of it on the expansion work and on the change of the thermal energy of the gas, so the factor is
// 1 - (c_adiabatic + c_heatcapacity) / heating. This holds also when the T_e solver limits T_e and the balance
// does not close. The factor is above 1 when the gas releases stored thermal energy.
auto set_radiative_energy_factor(const int nonemptymgi, const HeatingCoolingRates& heatingcoolingrates) -> double {
  const double heating = heatingcoolingrates.get_total_heating();
  const double cooling_nonradiative = heatingcoolingrates.cooling_adiabatic + heatingcoolingrates.cooling_heatcapacity;
  double factor = 1.;
  if (heating > 0.) {
    // A negative factor means that the gas absorbs more than the heating, and a large factor means that the
    // gas releases much more than the heating. The existing k-packets cannot represent these cases.
    constexpr double max_factor = 100.;
    factor = std::clamp(1. - (cooling_nonradiative / heating), 0., max_factor);
    if (factor <= 0. || factor >= max_factor) {
      printlnlog(
          "[warning] cell {} timestep {}: k-packet energy factor limited to {:g}. heating {:g}, adiabatic cooling "
          "{:g}, heat-capacity term {:g} [erg/s/cm^3]",
          grid::get_mgi_of_nonemptymgi(nonemptymgi), globals::timestep, factor, heating,
          heatingcoolingrates.cooling_adiabatic, heatingcoolingrates.cooling_heatcapacity);
    }
  }
  if (heating <= 0. && cooling_nonradiative < 0.) {
    // without heating there are no k-packets that can carry the released thermal energy
    printlnlog(
        "[warning] cell {} timestep {}: the gas releases thermal energy at {:g} [erg/s/cm^3] but the heating is zero, "
        "so the energy cannot enter the radiation field",
        grid::get_mgi_of_nonemptymgi(nonemptymgi), globals::timestep, -cooling_nonradiative);
  }
  radiative_energy_factor_allcells[nonemptymgi] = factor;
  return factor;
}

// A cell without a thermal balance keeps the full energy of its k-packets
void reset_radiative_energy_factor(const int nonemptymgi) { radiative_energy_factor_allcells[nonemptymgi] = 1.; }

DEVICE_FUNC auto get_radiative_energy_factor(const int nonemptymgi) -> double {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < grid::get_nonempty_npts_model());
  return radiative_energy_factor_allcells[nonemptymgi];
}

// handle a k-packet (e.g., in a thick cell) by emitting according to the planck function
DEVICE_FUNC void do_kpkt_blackbody(Packet& pkt) {
  const auto nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);

  pkt.e_cmf *= get_radiative_energy_factor(nonemptymgi);

  if (RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value() &&
      grid::thick_allcells[nonemptymgi] != grid::CellThickness::THICK) {
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
  // The diffusion step stands for the k-packet to r-packet cycles that a full treatment would follow. The
  // energy loss is the adiabatic loss of the radiation in those cycles (e proportional to 1/t). It is
  // separate from the expansion work of the gas in the thermal balance.
  pkt.e_cmf *= pkt.prop_time / t_current;
  pkt.prop_time = t_current;

  if (t_current >= t2) {
    return;
  }

  stats::increment(stats::Counter::INTERACTIONS);

  const auto nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);

  // The thermal balance of the cell decides which share of the thermal energy becomes radiation. The factor
  // applies on every pass through the thermal pool, also when the macro-atom returns the packet as a
  // k-packet. The heating rate of the balance counts that return as collisional heating, so each pass is
  // one flux event. The error of this choice is of second order in (1 - factor) and proportional to the
  // share of the k-packet energy that the macro-atom returns, which is small at nebular densities.
  pkt.e_cmf *= get_radiative_energy_factor(nonemptymgi);
  const std::span<const double> ion_cooling_contribs_thiscell = get_cell_ion_cooling_contribs(nonemptymgi);
  const double rndcool_ion = rng_uniform(get_rngstate(pkt)) * ion_cooling_contribs_thiscell.back();

  // Randomly select the occurring cooling process: first ion whose cumulative cooling exceeds rndcool_ion.
  // upper_bound (not lower_bound) so that a zero-cooling ion is never selected when rndcool_ion ties its
  // cumulative value exactly (which would give a zero total and an ill-defined process selection).
  const int uniqueionindex = int_index_upperbound(ion_cooling_contribs_thiscell, rndcool_ion);
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

  const auto ionoffset = index_upperbound(ion_contribs, rndcool_ion_process);
  if (ionoffset >= ncoolingterms_ion) [[unlikely]] {
    // The ion array of the cell gave this ion a positive cooling rate, so the cumulative sum of the cell cache
    // must end above the target. A total of zero, a negative total, or a NaN means that the two arrays came
    // from different cell states or from a negative cooling term. Device code has no log file, so only the host
    // writes the details. The assertion below reports the failure on both paths.
    MY_IF_HOST(const double ion_cooling_from_ionarray =
                   ion_cooling_contribs_thiscell[uniqueionindex] -
                   ((uniqueionindex > 0) ? ion_cooling_contribs_thiscell[uniqueionindex - 1] : 0.);
               printlnlog("[error] do_kpkt: cell {} timestep {}: no cooling process found for Z={} ionstage {}. Cell "
                          "cache total {:g}, target {:g}, ion array cooling rate {:g}, {} cooling terms",
                          grid::get_mgi_of_nonemptymgi(nonemptymgi), nts, get_atomicnumber(element),
                          get_ionstage(element, ion), C_ion_procsum, rndcool_ion_process, ion_cooling_from_ionarray,
                          ncoolingterms_ion););
  }
  assert_always(ionoffset < ncoolingterms_ion);
  const auto i = ionstart + ionoffset;

  const auto rndcoolingtype = coolinglist_type[i];
  const auto T_e = grid::Te_allcells[nonemptymgi];

  if (rndcoolingtype == CoolingType::FREEFREE) {
    // The k-packet converts directly into a r-packet by free-free emission.
    // Sample the comoving-frame frequency from the free-free emissivity (paper II, Sect. 5.4.3 eq. 41). The
    // emissivity is treated as frequency-independent apart from its exp(-h nu / k T_e) factor, so nu is
    // exponentially distributed with mean k T_e / h and can be drawn by inverting the CDF.
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
    // The k-packet converts directly into a r-packet by free-bound emission.
    const int lowerion = ion;
    const int lowerlevel = coolinglist_level[i];
    const int phixstargetindex = coolinglist_phixstargetindex[i];

    // Sample the comoving-frame frequency from the energy distribution of this recombination continuum
    // (paper II, Sect. 4.2.2), i.e. by inverting the CDF of the energy-weighted emissivity over the
    // continuum's frequency range.
    pkt.nu_cmf = select_continuum_nu(element, lowerion, lowerlevel, phixstargetindex, T_e, get_rngstate(pkt));

    // and then emit the packet randomly in the comoving frame
    emit_rpkt(pkt);

    pkt.next_trans = -1;  // FLAG: transition history here not important, cont. process
    stats::increment(stats::Counter::K_STAT_TO_R_FB);
    pkt.emissiontype = get_emtype_continuum(element, lowerion, lowerlevel, phixstargetindex);
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
          col_excitation_ratecoeff(T_e, clumpednne, epsilon_trans, upper_statweight, statweight, alltransindex) *
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
    const int upper = get_phixsupperlevel(element, ion, coolinglist_level[i], coolinglist_phixstargetindex[i]);

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
