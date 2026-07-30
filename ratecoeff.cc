// Photoionisation, recombination, and bound-free heating/cooling rate coefficients:
// precalculates tables of the temperature-dependent integrals at startup and evaluates the
// radiation-field-dependent photoionisation coefficients during the simulation.

#include "ratecoeff.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <format>
#include <ios>
#include <span>
#include <sstream>
#include <string>
#include <utility>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "integrator.h"
#include "ltepop.h"
#include "macroatom.h"
#include "mpi_logging.h"
#include "radfield.h"
#include "random.h"

namespace {
constexpr double RATECOEFF_INTEGRAL_ACCURACY = 1e-3;

const double T_step_log = (std::log(MAXTEMP) - std::log(MINTEMP)) / (TABLESIZE - 1.);

const auto temperature_grid = []() {
  std::array<double, TABLESIZE + 1> grid{};
  for (auto i = 0UZ; i < grid.size(); i++) {
    grid[i] = MINTEMP * std::exp(i * T_step_log);
  }
  return grid;
}();

// Index of the first temperature grid point above the given temperature, matching
// std::ranges::upper_bound(temperature_grid, temperature). The grid is log-uniform, so the
// index can be computed directly instead of with a binary search. The correction loops
// (almost always zero iterations) keep the result identical to upper_bound under floating-point
// rounding of the analytic estimate.
[[gnu::pure]] [[nodiscard]] auto get_temperature_gridupperindex(const double temperature) -> int {
  const auto gridsize = static_cast<int>(temperature_grid.size());
  int index = std::clamp(static_cast<int>(std::log(temperature / MINTEMP) / T_step_log) + 1, 0, gridsize);
  while (index > 0 && temperature_grid[index - 1] > temperature) {
    index--;
  }
  while (index < gridsize && temperature_grid[index] <= temperature) {
    index++;
  }
  return index;
}

MPI_shared_array<const float> ion_alpha_sp;  // size is nincludedions * TABLESIZE, indexed
                                             // by (uniqueionindex * TABLESIZE) + temperatureindex

// the following spans are indexed by get_bflutindex()
MPI_shared_array<double> spontrecombcoeffs{};  // indexed by get_bflutindex()
MPI_shared_array<double> corrphotoioncoeffs{};  // for USE_LUT_PHOTOION = true
MPI_shared_array<double> bfcooling_coeffs{};

// Integrand to calculate the rate coefficient for spontaneous recombination
auto alpha_sp_integrand(const double nu_minus_nu_edge, const double nu_edge, const float T_e,
                        const std::span<const float> photoion_xs) -> double {
  const auto sigma_bf = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_minus_nu_edge + nu_edge);
  // the variable of integration has been changed from nu to nu_minus_nu_edge = nu - nu_edge
  // to get a cancellation with part of the saha factor
  return (2 / CLIGHTSQUARED) * sigma_bf * pow2(nu_edge + nu_minus_nu_edge) * exp(-HOVERKB * nu_minus_nu_edge / T_e);
}

// Energy-weighted integrand used to sample the frequency of a spontaneous free-bound emission
auto alpha_sp_E_integrand(const double nu_minus_nu_edge, const double nu_edge, const float T_e,
                          const std::span<const float> photoion_xs) -> double {
  const double nu = nu_edge + nu_minus_nu_edge;
  const auto sigma_bf = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu);
  // The omitted exp(-h nu_edge / kT_e) factor is constant across the continuum and cancels from the normalised CDF.
  return (2 / CLIGHTSQUARED) * sigma_bf * pow3(nu) / nu_edge * exp(-HOVERKB * nu_minus_nu_edge / T_e);
}

// Integrand to calculate the rate coefficient for photoionisation corrected for stimulated recombination.
auto gammacorr_integrand(const double nu, const double nu_edge, const float temperature,
                         const std::span<const float> photoion_xs) -> double {
  const auto sigma_bf = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu);

  // The correction factor for stimulated emission in gammacorr is set to its LTE value. Because the T_e dependence of
  // gammacorr is weak, this correction may be evaluated at T_R!

  // Dependence on dilution factor W is linear. This allows to set it here to
  // 1. and scale to its actual value later on.
  // Assumption T_e = T_R makes n_kappa/n_i * (n_i/n_kappa)* = 1
  return sigma_bf * (1. / H) / nu * radfield::planck(nu, temperature) * (1 - exp(-HOVERKB * nu / temperature));
}

// Integrand to precalculate the bound-free cooling rate coefficient
auto bfcooling_integrand(const double nu_minus_nu_edge, const double nu_edge, const float T_e,
                         const std::span<const float> photoion_xs) -> double {
  const float sigma_bf = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_minus_nu_edge + nu_edge);

  return sigma_bf * nu_minus_nu_edge * (2 * H / CLIGHTSQUARED) * (nu_minus_nu_edge + nu_edge) *
         (nu_minus_nu_edge + nu_edge) * exp(-HOVERKB * nu_minus_nu_edge / T_e);
}

[[gnu::pure]] [[nodiscard]] inline auto get_bflutindex(const int temperatureindex, const int uniquelevelindex,
                                                       const int phixstargetindex) -> int {
  // continuum-major layout so that the two temperature samples read by an interpolation are adjacent
  const int contindex = globals::alllevels.bflist_start[uniquelevelindex] + phixstargetindex;
  const int bflutindex = (contindex * TABLESIZE) + temperatureindex;
  assert_testmodeonly(bflutindex >= 0);
  assert_testmodeonly(bflutindex < TABLESIZE * globals::nbfcontinua);
  return bflutindex;
}

[[gnu::pure]] [[nodiscard]] inline auto get_bflutindex(const int temperatureindex, const int element, const int ion,
                                                       const int level, const int phixstargetindex) -> int {
  return get_bflutindex(temperatureindex, get_uniquelevelindex(element, ion, level), phixstargetindex);
}

void precalculate_rate_coefficient_integrals() {
  // we're writing to shared memory, so we need to synchronise
  MPI_Barrier_node();

  // Calculate the rate coefficients for each level of each ion of each element
  for (int element = 0; element < get_nelements(); element++) {
    const int atomic_number = get_atomicnumber(element);
    const int nions = get_nions(element);
    if (nions == 0) {
      continue;
    }
    printlog("Performing rate integrals for Z = {}: ion stages", atomic_number);
    for (int ion = 0; ion < nions - 1; ion++) {
      printlog(" {}", get_ionstage(element, ion));
    }
    printlnlog("");
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels_ionising(element, ion);

      for (int level = 0; level < nlevels; level++) {
        // coefficients are stored in node shared memory, so divide up the work on the node
        if ((level % globals::node_nprocs) != globals::rank_in_node) {
          continue;
        }
        const double statw_lower = stat_weight(element, ion, level);
        const int nphixstargets = get_nphixstargets(element, ion, level);
        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
          const int upperlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
          const double phixstargetprobability = get_phixsprobability(element, ion, level, phixstargetindex);
          const double statw_upper = stat_weight(element, ion + 1, upperlevel);

          const double E_threshold = get_phixs_threshold(element, ion, level, phixstargetindex);
          const double nu_threshold = E_threshold / H;
          const double nu_max_phixs =
              nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table
          // Loop over the temperature grid
          for (int temperatureindex = 0; temperatureindex < TABLESIZE; temperatureindex++) {
            const int bflutindex = get_bflutindex(temperatureindex, element, ion, level, phixstargetindex);
            double error{NAN};
            const auto temperature = static_cast<float>(temperature_grid[temperatureindex]);

            const double modified_sahafact = SAHACONST * statw_lower / statw_upper * std::pow(temperature, -1.5);
            assert_always(modified_sahafact >= 0.);
            assert_always(std::isfinite(modified_sahafact));

            assert_always(!get_phixs_table(element, ion, level).empty());
            // the threshold of the first target gives nu of the first phixstable point
            const auto photoion_xs = get_phixs_table(element, ion, level);

            // Spontaneous recombination and bf-cooling coefficient don't depend on the radiation field
            const auto alpha_sp =
                FOURPI * modified_sahafact * phixstargetprobability *
                integrator(
                    [&](const double nu_minus_nu_edge) {
                      return alpha_sp_integrand(nu_minus_nu_edge, nu_threshold, temperature, photoion_xs);
                    },
                    0, nu_max_phixs - nu_threshold, RATECOEFF_INTEGRAL_ACCURACY, &error);

            assert_always(std::isfinite(alpha_sp) && alpha_sp >= 0);
            spontrecombcoeffs[bflutindex] = alpha_sp;

            if constexpr (USE_LUT_PHOTOION) {
              auto gammacorr = integrator(
                  [&](const double nu) { return gammacorr_integrand(nu, nu_threshold, temperature, photoion_xs); },
                  nu_threshold, nu_max_phixs, RATECOEFF_INTEGRAL_ACCURACY, &error);
              gammacorr *= FOURPI * phixstargetprobability;
              assert_always(gammacorr >= 0);
              corrphotoioncoeffs[bflutindex] = gammacorr;
            }
            const auto this_bfcooling_coeff =
                FOURPI * modified_sahafact * phixstargetprobability *
                integrator(
                    [&](const double nu_minus_nu_edge) {
                      return bfcooling_integrand(nu_minus_nu_edge, nu_threshold, temperature, photoion_xs);
                    },
                    0, nu_max_phixs - nu_threshold, RATECOEFF_INTEGRAL_ACCURACY, &error);

            assert_always(std::isfinite(this_bfcooling_coeff) && this_bfcooling_coeff >= 0);
            bfcooling_coeffs[bflutindex] = this_bfcooling_coeff;
          }  // temperature loop
        }  // phixstarget loop
      }  // level loop
    }  // ion loop
  }

  MPI_Barrier_node();
}

// multiply the cross sections associated with a level by some factor and
// also update the quantities integrated from (and proportional to) the cross sections
void scale_level_phixs(const int element, const int ion, const int level, const double factor) {
  // if we store the data in node shared memory, then only one rank should update it
  if (globals::rank_in_node == 0) {
    const int nphixstargets = get_nphixstargets(element, ion, level);
    if (nphixstargets == 0) {
      return;
    }

    const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
    const auto phixsstart = globals::alllevels.phixsstart[uniquelevelindex];

    const auto phixstable =
        globals::allphixs.subspan(static_cast<ptrdiff_t>(phixsstart) * globals::NPHIXSPOINTS, globals::NPHIXSPOINTS);
    for (int n = 0; n < globals::NPHIXSPOINTS; n++) {
      phixstable[n] = static_cast<float>(phixstable[n] * factor);
    }

    for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
      for (int tempindex = 0; tempindex < TABLESIZE; tempindex++) {
        const auto bflutindex = get_bflutindex(tempindex, element, ion, level, phixstargetindex);
        spontrecombcoeffs[bflutindex] = spontrecombcoeffs[bflutindex] * factor;

        if constexpr (USE_LUT_PHOTOION) {
          corrphotoioncoeffs[bflutindex] *= factor;
        }

        bfcooling_coeffs[bflutindex] *= factor;
      }
    }
  }
}

// calibrate the recombination rates to tabulated values by scaling the photoionisation cross sections
void read_recombrate_file() {
  if (!std::filesystem::exists("recombrates.txt")) {
    printlnlog("No recombrates.txt file found. Skipping recombination rate scaling...");
    return;
  }

  printlnlog("Reading recombination rate file (recombrates.txt)...");
  auto recombrate_file = fstream_required("recombrates.txt", std::ios::in);

  const float Te_estimate = RECOMBCALIBRATION_T_ELEC;
  const double log_Te_estimate = log10(RECOMBCALIBRATION_T_ELEC);

  printlnlog("Calibrating recombination rates for a temperature of {:.1f} [K]", Te_estimate);

  struct RRCRow {
    double log_Te;
    double rrc_low_n;
    double rrc_total;
  };

  int atomicnumber = 0;
  int upperionstage = 0;
  int tablerows = 0;

  std::string line;
  while (get_noncommentline(recombrate_file, line)) {
    assert_always(std::stringstream(line) >> atomicnumber >> upperionstage >> tablerows);
    assert_always(tablerows >= 0);
    RRCRow T_highestbelow = {.log_Te = 0, .rrc_low_n = 0, .rrc_total = 0};
    RRCRow T_lowestabove = {.log_Te = 0, .rrc_low_n = 0, .rrc_total = 0};
    T_highestbelow.log_Te = -1;
    T_lowestabove.log_Te = -1;
    for (int i = 0; i < tablerows; i++) {
      RRCRow row{};
      // without these checks a truncated or malformed table would leave the previous line's values
      // in place and silently calibrate the wrong ion
      assert_always(get_noncommentline(recombrate_file, line));
      assert_always(std::stringstream(line) >> row.log_Te >> row.rrc_low_n >> row.rrc_total);

      if (row.log_Te < log_Te_estimate && row.log_Te > T_highestbelow.log_Te) {
        T_highestbelow.log_Te = row.log_Te;
        T_highestbelow.rrc_low_n = row.rrc_low_n;
        T_highestbelow.rrc_total = row.rrc_total;
      }

      if (row.log_Te > log_Te_estimate && (row.log_Te < T_lowestabove.log_Te || T_lowestabove.log_Te < 0)) {
        T_lowestabove.log_Te = row.log_Te;
        T_lowestabove.rrc_low_n = row.rrc_low_n;
        T_lowestabove.rrc_total = row.rrc_total;
      }
    }
    const int element = get_elementindex(atomicnumber);
    if (element >= 0) {
      const int ion = upperionstage - get_ionstage(element, 0);  // the index of the upper ion
      if (ion > 0 && ion < get_nions(element)) {
        printlnlog("Z={} ionstage {}->{}", atomicnumber, upperionstage, upperionstage - 1);
        assert_always(T_highestbelow.log_Te > 0);
        assert_always(T_lowestabove.log_Te > 0);

        const int nlevels = get_nlevels_ionising(element, ion - 1);

        // scale_level_phixs() writes node-shared memory (allphixs and the derived rate coefficient
        // tables) from rank_in_node == 0 only, while every rank on the node is reading those same
        // arrays in calculate_ionrecombcoeff(). Bracket the writes with node barriers so that no
        // rank can read a partially rescaled table.
        const auto scale_levels = [element, ion, nlevels](const int firstlevel, const double multiplier) {
          assert_always(std::isfinite(multiplier) && multiplier >= 0.);
          MPI_Barrier_node();
          for (int level = firstlevel; level < nlevels; level++) {
            scale_level_phixs(element, ion - 1, level, multiplier);
          }
          MPI_Barrier_node();
        };

        const double x = (log_Te_estimate - T_highestbelow.log_Te) / (T_lowestabove.log_Te - T_highestbelow.log_Te);
        const double input_rrc_low_n = std::lerp(T_highestbelow.rrc_low_n, T_lowestabove.rrc_low_n, x);
        const double input_rrc_total = std::lerp(T_highestbelow.rrc_total, T_lowestabove.rrc_total, x);

        constexpr bool assume_lte = true;
        constexpr bool per_groundmultipletpop = true;
        constexpr bool collisional_not_radiative = false;

        double rrc = calculate_ionrecombcoeff(-1, Te_estimate, element, ion, assume_lte, collisional_not_radiative,
                                              false, per_groundmultipletpop);
        printlnlog("    rrc (initial): {:10.3e} [cm^3/s]", rrc);

        if (!(rrc > 0.)) {
          // no recombination into this ion from the atomic data (e.g. the ion below has no
          // photoionisation targets), so there is nothing to calibrate and the multipliers below
          // would all be a division by zero
          printlnlog("    rrc is not positive, so skipping the recombination rate calibration for this ion");
          continue;
        }

        if (input_rrc_low_n >= 0)  // if it's < 0, ignore it
        {
          printlnlog("  input_rrc_low_n: {:10.3e}", input_rrc_low_n);

          const double phixs_multiplier = input_rrc_low_n / rrc;
          if (!(phixs_multiplier >= 0.05) || phixs_multiplier >= 2.0) {
            printlnlog("    Not scaling phixs of all levels by {:.3f} (because < 0.05 or >= 2.0)", phixs_multiplier);
          } else {
            printlnlog("    scaling phixs of all levels by {:.3f}", phixs_multiplier);

            scale_levels(0, phixs_multiplier);

            rrc = calculate_ionrecombcoeff(-1, Te_estimate, element, ion, assume_lte, collisional_not_radiative, false,
                                           per_groundmultipletpop);
            printlnlog("    rrc (after low-n scaling): {:10.3e} [cm^3/s]", rrc);
          }
        }

        // hopefully the RRC now matches the low_n value well, if it was defined
        // Next, use the superlevel recombination rates to make up the excess needed to reach the total RRC

        printlnlog("  input_rrc_total: {:10.3e} [cm^3/s]", input_rrc_total);

        if (input_rrc_total < 0) {
          // negative means no tabulated total, in the same way that a negative rrc_low_n is ignored above
          printlnlog("    input_rrc_total is negative, so not scaling to the total recombination rate");
        } else if (rrc < input_rrc_total) {
          const double rrc_superlevel = calculate_ionrecombcoeff(
              -1, Te_estimate, element, ion, assume_lte, collisional_not_radiative, true, per_groundmultipletpop);
          printlnlog("  rrc(superlevel): {:10.3e} [cm^3/s]", rrc_superlevel);

          if (rrc_superlevel > 0) {
            const double phixs_multiplier_superlevel = 1.0 + ((input_rrc_total - rrc) / rrc_superlevel);
            printlnlog("    scaling phixs of levels in the superlevel by {:.3f}", phixs_multiplier_superlevel);

            const int first_superlevel_level = get_nlevels_excited_nlte(element, ion - 1) + 1;
            scale_levels(first_superlevel_level, phixs_multiplier_superlevel);
          } else {
            printlnlog("There is no superlevel recombination, so multiplying all levels instead");
            const double phixs_multiplier = input_rrc_total / rrc;
            printlnlog("    scaling phixs of all levels by {:.3f}", phixs_multiplier);

            scale_levels(0, phixs_multiplier);
          }
        } else {
          printlnlog("    rrc {:10.3e} >= input_rrc_total {:10.3e} [cm^3/s]", rrc, input_rrc_total);
          const double phixs_multiplier = input_rrc_total / rrc;
          printlnlog("    scaling phixs of all levels by {:.3f}", phixs_multiplier);

          scale_levels(0, phixs_multiplier);
        }

        rrc = calculate_ionrecombcoeff(-1, Te_estimate, element, ion, assume_lte, collisional_not_radiative, false,
                                       per_groundmultipletpop);
        printlnlog("    rrc (final): {:10.3e} [cm^3/s]", rrc);
      }
    }
  }
}

void precalculate_ion_alpha_sp() {
  auto temp_ion_alpha_sp = MPI_shared_array<float>(get_includedions() * TABLESIZE, 0.);
  if (globals::rank_in_node == 0) {
    for (int tempindex = 0; tempindex < TABLESIZE; tempindex++) {
      const auto T_e = static_cast<float>(temperature_grid[tempindex]);
      for (int element = 0; element < get_nelements(); element++) {
        const int nions = get_nions(element) - 1;
        for (int ion = 0; ion < nions; ion++) {
          const auto uniqueionindex = get_uniqueionindex(element, ion);
          const int nionisinglevels = get_nlevels_ionising(element, ion);
          double alpha_sp = 0.;
          for (int level = 0; level < nionisinglevels; level++) {
            const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
            const auto nphixstargets = get_nphixstargets(uniquelevelindex);
            for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
              const double alpha_sp_level = get_spontrecombcoeff(uniquelevelindex, phixstargetindex, T_e);
              alpha_sp += alpha_sp_level;
            }
          }
          temp_ion_alpha_sp[(uniqueionindex * TABLESIZE) + tempindex] = static_cast<float>(alpha_sp);
        }
      }
    }
  }
  assert_always(ion_alpha_sp.empty());
  ion_alpha_sp = std::move(temp_ion_alpha_sp);
  MPI_Barrier_node();
}

// Integrand to calculate the rate coefficient for photoionisation. Corrected for stimulated recombination.
auto integrand_corrphotoioncoeff_custom_radfield(const double nu_minus_nu_edge, const double nu_edge,
                                                 const double modified_departure_ratio,
                                                 const std::span<const float> photoion_xs, const float T_e,
                                                 const int nonemptymgi) -> double {
  double corrfactor = 1. - (modified_departure_ratio * exp(-HOVERKB * nu_minus_nu_edge / T_e));
  if (corrfactor < 0) {
    corrfactor = 0.;
  }

  const float sigma_bf = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_minus_nu_edge + nu_edge);

  const double Jnu = radfield::radfield(nu_minus_nu_edge + nu_edge, nonemptymgi);

  // TODO: MK thesis page 41, use population ratios and Te?
  return (1. / H) * sigma_bf / (nu_minus_nu_edge + nu_edge) * Jnu * corrfactor;
}

auto calculate_corrphotoioncoeff_integral(const int element, const int ion, const int level, const int phixstargetindex,
                                          const int nonemptymgi, const bool use_cellcache) -> double {
  constexpr double epsrel = 1e-3;

  const auto loweruniquelevelindex = get_uniquelevelindex(element, ion, level);
  const double nu_threshold = (1. / H) * get_phixs_threshold(element, ion, level, phixstargetindex);
  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

  const auto T_e = grid::Te_allcells[nonemptymgi];

  // stimulated recombination is negative photoionisation
  const double nnlevel = use_cellcache ? get_cellcache_levelpop(nonemptymgi, loweruniquelevelindex)
                                       : calculate_levelpop(nonemptymgi, element, ion, level);
  const auto clumpednne = grid::get_nne(nonemptymgi) * grid::get_clumpfactor(nonemptymgi);

  const int upperionlevel = get_phixsupperlevel(loweruniquelevelindex, phixstargetindex);
  const auto upperuniquelevelindex = get_uniquelevelindex(element, ion + 1, upperionlevel);
  const double modified_sahafact =
      SAHACONST * stat_weight(loweruniquelevelindex) / stat_weight(upperuniquelevelindex) * std::pow(T_e, -1.5);
  const double nnupperionlevel = use_cellcache ? get_cellcache_levelpop(nonemptymgi, upperuniquelevelindex)
                                               : calculate_levelpop(nonemptymgi, element, ion + 1, upperionlevel);

  double modified_departure_ratio = nnlevel > 0. ? nnupperionlevel / nnlevel * clumpednne * modified_sahafact : 1.;
  if (!std::isfinite(modified_departure_ratio)) {
    modified_departure_ratio = 0.;
  }

  double error = 0.;
  const auto photoion_xs = get_phixs_table(loweruniquelevelindex);

  const auto gammacorr =
      FOURPI * get_phixsprobability(loweruniquelevelindex, phixstargetindex) *
      integrator(
          [&](const double nu_minus_nu_edge) {
            return integrand_corrphotoioncoeff_custom_radfield(nu_minus_nu_edge, nu_threshold, modified_departure_ratio,
                                                               photoion_xs, T_e, nonemptymgi);
          },
          0, nu_max_phixs - nu_threshold, epsrel, &error);
  assert_always(std::isfinite(gammacorr));

  return gammacorr;
}

template <typename T, typename U>
[[nodiscard]] auto lerp_or_last(const std::span<T> table, const int uniquelevelindex, const int phixstargetindex,
                                const U temperature) {
  const auto upperindex = get_temperature_gridupperindex(temperature);
  if (upperindex == 0) {
    return table[get_bflutindex(0, uniquelevelindex, phixstargetindex)];
  }
  if (upperindex < TABLESIZE) {
    const double T_lower = temperature_grid[upperindex - 1];
    const double T_upper = temperature_grid[upperindex];

    const double f_lower = table[get_bflutindex(upperindex - 1, uniquelevelindex, phixstargetindex)];
    const double f_upper = table[get_bflutindex(upperindex, uniquelevelindex, phixstargetindex)];
    return (f_lower + ((f_upper - f_lower) / (T_upper - T_lower) * (temperature - T_lower)));
  }
  return table[get_bflutindex(TABLESIZE - 1, uniquelevelindex, phixstargetindex)];
}

}  // anonymous namespace

void setup_photoion_luts() {
  assert_always(globals::nbfcontinua > 0);
  const auto lutsize = static_cast<ptrdiff_t>(TABLESIZE) * globals::nbfcontinua;
  size_t mem_usage_photoionluts = 2 * lutsize * sizeof(double);

  spontrecombcoeffs = MPI_shared_array<double>(lutsize, 0.);

  if constexpr (USE_LUT_PHOTOION) {
    corrphotoioncoeffs = MPI_shared_array<double>(lutsize, 0.);
    mem_usage_photoionluts += lutsize * sizeof(double);
  }

  bfcooling_coeffs = MPI_shared_array<double>(lutsize, 0.);

  printlnlog(
      "[info] mem_usage: lookup tables derived from photoionisation (spontrecombcoeff, bfcooling and "
      "corrphotoioncoeff if USE_LUT_PHOTOION) occupy {:.3f} MB",
      mem_usage_photoionluts / 1024. / 1024.);
}

DEVICE_FUNC auto select_continuum_nu(int element, const int lowerion, const int lower, const int upperionlevel,
                                     float T_e, rngstate_type& rngstate) -> double {
  const auto lower_uniquelevelindex = get_uniquelevelindex(element, lowerion, lower);
  const int phixstargetindex = get_phixstargetindex(lower_uniquelevelindex, upperionlevel);
  const double E_threshold = get_phixs_threshold(element, lowerion, lower, phixstargetindex);
  const double nu_threshold = (1. / H) * E_threshold;

  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

  const int npieces = globals::NPHIXSPOINTS;

  const auto photoion_xs = get_phixs_table(lower_uniquelevelindex);

  const double zrand = 1. - rng_uniform(rngstate);  // Make sure that 0 < zrand <= 1

  const double nu_range = nu_max_phixs - nu_threshold;
  const double deltanu = nu_range / npieces;
  double error{NAN};

  // integral of the energy-weighted emissivity over the whole continuum, i.e. the normalisation of the
  // distribution that we are sampling a frequency from
  const auto emissivity_integral_total = integrator<31>(
      [&](const double nu_minus_nu_edge) {
        return alpha_sp_E_integrand(nu_minus_nu_edge, nu_threshold, T_e, photoion_xs);
      },
      0., nu_range, RATECOEFF_INTEGRAL_ACCURACY, &error);

  assert_testmodeonly(std::isfinite(emissivity_integral_total));
  assert_testmodeonly(emissivity_integral_total > 0.);
  if (!(emissivity_integral_total > 0.) || !std::isfinite(emissivity_integral_total)) {
    return nu_threshold;
  }

  // emissivity_tailintegral is the part of the integral above the current piece's lower boundary, so it
  // decreases from emissivity_integral_total to zero as we step up in frequency. _prev is its value at the
  // previous piece boundary.
  double emissivity_tailintegral_prev = emissivity_integral_total;
  double emissivity_tailintegral = emissivity_integral_total;

  int i = 1;
  for (; i < npieces; i++) {
    emissivity_tailintegral_prev = emissivity_tailintegral;
    const double nu_minus_nu_edge_low = i * deltanu;

    emissivity_tailintegral = integrator<31>(
        [&](const double nu_minus_nu_edge) {
          return alpha_sp_E_integrand(nu_minus_nu_edge, nu_threshold, T_e, photoion_xs);
        },
        nu_minus_nu_edge_low, nu_range, RATECOEFF_INTEGRAL_ACCURACY, &error);

    if (zrand >= emissivity_tailintegral / emissivity_integral_total) {
      break;
    }
  }

  double nuoffset = 0.;
  if (i < npieces) {
    // the loop found the piece containing the target value, so interpolate between the remaining
    // integrals at the piece boundaries
    nuoffset = (emissivity_tailintegral != emissivity_tailintegral_prev)
                   ? ((emissivity_integral_total * zrand) - emissivity_tailintegral_prev) /
                         (emissivity_tailintegral - emissivity_tailintegral_prev) * deltanu
                   : 0.;
  } else if (emissivity_tailintegral > 0.) {
    // the loop completed without a break, so the target lies in the topmost piece, where the remaining
    // integral falls from emissivity_tailintegral at the piece's lower boundary to zero at nu_max_phixs.
    // Interpolating with the previous piece's values here would extrapolate beyond nu_max_phixs.
    nuoffset = (emissivity_tailintegral - (emissivity_integral_total * zrand)) / emissivity_tailintegral * deltanu;
  }
  const double nu_sampled = nu_threshold + ((i - 1) * deltanu) + nuoffset;

  assert_testmodeonly(std::isfinite(nu_sampled));
  assert_testmodeonly(nu_sampled >= nu_threshold);
  assert_testmodeonly(nu_sampled <= nu_max_phixs);

  return nu_sampled;
}

// Get an ion's rate coefficient for spontaneous recombination in LTE
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_ion_spontrecombcoeff(const int uniqueionindex, const float T_e)
    -> double {
  const auto upperindex = get_temperature_gridupperindex(T_e);
  if (upperindex == 0) {
    return ion_alpha_sp[uniqueionindex * TABLESIZE];
  }
  if (upperindex < TABLESIZE) {
    const double T_lower = temperature_grid[upperindex - 1];
    const double T_upper = temperature_grid[upperindex];

    const double f_lower = ion_alpha_sp[(uniqueionindex * TABLESIZE) + upperindex - 1];
    const double f_upper = ion_alpha_sp[(uniqueionindex * TABLESIZE) + upperindex];

    return f_lower + ((f_upper - f_lower) / (T_upper - T_lower) * (T_e - T_lower));
  }
  return ion_alpha_sp[(uniqueionindex * TABLESIZE) + (TABLESIZE - 1)];
}

// Return a level's rate coefficient for spontaneous recombination in LTE
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_spontrecombcoeff(const int uniquelevelindex,
                                                                  const int phixstargetindex, const float T_e)
    -> double {
  return lerp_or_last(std::span{spontrecombcoeffs}, uniquelevelindex, phixstargetindex, T_e);
}

// multiply by upper ion population (or ground population if per_groundmultipletpop is true) and nne to get a rate
auto calculate_ionrecombcoeff(const int nonemptymgi, const float T_e, const int element, const int upperion,
                              const bool assume_lte, const bool collisional_not_radiative,
                              const bool lower_superlevel_only, const bool per_groundmultipletpop) -> double {
  if (upperion <= 0) {
    return 0.;
  }

  const int lowerion = upperion - 1;
  assert_always(lowerion < get_nions(element) - 1);

  double nnupperion = 0;
  int upper_nlevels = 0;
  if (per_groundmultipletpop) {
    // assume that photoionisation of the ion below is only to the ground multiplet levels of the current ion
    upper_nlevels = get_nlevels_groundterm(element, lowerion + 1);
  } else {
    upper_nlevels = get_nlevels(element, lowerion + 1);
  }

  // population of an upper-ion level: LTE (Boltzmann relative to ground) or the full NLTE level population
  const auto calc_nnupperlevel = [&](const int upper) -> double {
    if (assume_lte) {
      const double T_exc = T_e;
      const double E_level = epsilon(element, lowerion + 1, upper);
      const double E_ground = epsilon(element, lowerion + 1, 0);
      const double nnground = (nonemptymgi >= 0) ? get_groundlevelpop(nonemptymgi, element, lowerion + 1) : 1.;

      return (nnground * stat_weight(element, lowerion + 1, upper) / stat_weight(element, lowerion + 1, 0) *
              exp(-(E_level - E_ground) / KB / T_exc));
    }
    return calculate_levelpop(nonemptymgi, element, lowerion + 1, upper);
  };

  for (int upper = 0; upper < upper_nlevels; upper++) {
    nnupperion += calc_nnupperlevel(upper);
  }

  if (nnupperion <= 0.) {
    return 0.;
  }

  // this gets divided and cancelled out in the radiative case anyway
  const auto clumpednne = (nonemptymgi >= 0) ? grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi) : 1.F;
  double alpha = 0.;
  const int maxrecombininglevel = get_maxrecombininglevel(element, lowerion + 1);
  for (int upper = 0; upper <= maxrecombininglevel; upper++) {
    const double nnupperlevel = calc_nnupperlevel(upper);
    for (int lower = 0; lower < get_nlevels(element, lowerion); lower++) {
      if (lower_superlevel_only && (!level_isinsuperlevel(element, lowerion, lower))) {
        continue;
      }

      double recomb_coeff = 0.;
      if (collisional_not_radiative) {
        const double epsilon_trans = epsilon(element, lowerion + 1, upper) - epsilon(element, lowerion, lower);
        recomb_coeff += col_recombination_ratecoeff(T_e, clumpednne, element, upperion, upper, lower, epsilon_trans);
      } else {
        recomb_coeff += rad_recombination_ratecoeff(T_e, clumpednne, element, lowerion + 1, upper, lower);
      }

      const double alpha_level = recomb_coeff / clumpednne;
      const double alpha_ion_contrib = alpha_level * nnupperlevel / nnupperion;
      alpha += alpha_ion_contrib;
    }
  }

  return alpha;
}

// Precalculates the rate coefficients for stimulated and spontaneous
// recombination and photoionisation on a given temperature grid using integration.
// NB: with the nebular approximation they only depend on T_e, T_R and W.
// W is easily factored out. For stimulated recombination we must assume
// T_e = T_R for this precalculation.
void ratecoefficients_init() {
  precalculate_rate_coefficient_integrals();

  read_recombrate_file();

  precalculate_ion_alpha_sp();
}

// Returns the (stimulated recombination corrected) photoionisation rate coefficient.
auto get_corrphotoioncoeff_ana(int element, const int ion, const int level, const int phixstargetindex, const float T_R)
    -> double {
  assert_always(USE_LUT_PHOTOION);
  const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
  return lerp_or_last(std::span{corrphotoioncoeffs}, uniquelevelindex, phixstargetindex, T_R);
}

DEVICE_FUNC auto get_bfcoolingcoeff(const int element, const int lowerion, const int lowerionlevel,
                                    const int phixstargetindex, const float T_e) -> double {
  return lerp_or_last(std::span{bfcooling_coeffs}, get_uniquelevelindex(element, lowerion, lowerionlevel),
                      phixstargetindex, T_e);
}

// Return the photoionisation rate coefficient (corrected for stimulated emission)
DEVICE_FUNC auto get_corrphotoioncoeff(const int element, const int ion, const int level, const int phixstargetindex,
                                       const int nonemptymgi, const bool use_cellcache) -> double {
  const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
  const auto allphixstargetindex = get_allphixstargetindex(uniquelevelindex, phixstargetindex);
  double gammacorr =
      use_cellcache ? get_cellcache(nonemptymgi).allphixstargets_corrphotoioncoeff[allphixstargetindex] : -1;

  if (!use_cellcache || gammacorr < 0) {
    if (DETAILED_BF_ESTIMATORS_ON && globals::timestep >= DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP) {
      gammacorr = radfield::get_bfrate_estimator(element, ion, level, phixstargetindex, nonemptymgi);
      // gammacorr will be -1 if no estimators available
    }

    if (!DETAILED_BF_ESTIMATORS_ON || gammacorr < 0) {
      if constexpr (!USE_LUT_PHOTOION) {
        gammacorr =
            calculate_corrphotoioncoeff_integral(element, ion, level, phixstargetindex, nonemptymgi, use_cellcache);
      } else {
        const double W = grid::W_allcells[nonemptymgi];
        const double T_R = grid::TR_allcells[nonemptymgi];

        gammacorr = W * lerp_or_last(std::span{corrphotoioncoeffs}, uniquelevelindex, phixstargetindex, T_R);
        const int index_in_groundlevelcontestimator = globals::alllevels.closestgroundlevelcont[uniquelevelindex];
        if (index_in_groundlevelcontestimator >= 0) {
          gammacorr *= globals::corrphotoionrenorm[(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) +
                                                   index_in_groundlevelcontestimator];
        }
      }
    }
    if (use_cellcache) {
      get_cellcache(nonemptymgi).allphixstargets_corrphotoioncoeff[allphixstargetindex] = gammacorr;
    }
  }

  return gammacorr;
}

auto iongamma_is_zero(const int nonemptymgi, const int element, const int ion) -> bool {
  const int nions = get_nions(element);
  if (ion >= nions - 1) {
    return true;
  }

  if (!elem_has_nlte_levels(element)) {
    const auto groundcontindex = get_groundcontindex(element, ion);
    if (groundcontindex < 0) {
      return true;
    }
    return (globals::gammaestimator[(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) +
                                    groundcontindex] == 0);
  }

  const auto T_e = grid::Te_allcells[nonemptymgi];
  const auto clumpednne = grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi);

  for (int level = 0; level < get_nlevels(element, ion); level++) {
    const double nnlevel = calculate_levelpop(nonemptymgi, element, ion, level);
    if (nnlevel == 0.) {
      continue;
    }
    const int nphixstargets = get_nphixstargets(element, ion, level);
    for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
      const int upperlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);

      if (nnlevel * get_corrphotoioncoeff(element, ion, level, phixstargetindex, nonemptymgi, false) > 0.) {
        return false;
      }

      const double epsilon_trans = epsilon(element, ion + 1, upperlevel) - epsilon(element, ion, level);

      if (nnlevel * col_ionisation_ratecoeff(T_e, clumpednne, element, ion, level, phixstargetindex, epsilon_trans) >
          0) {
        return false;
      }
    }
  }
  return true;
}

// ionisation rate coefficient. multiply by get_groundlevelpop to get a rate [s^-1]
auto calculate_iongamma_per_gspop(const int nonemptymgi, const int element, const int ion) -> double {
  const int nions = get_nions(element);
  if (ion >= nions - 1) {
    return 0.;
  }

  const auto T_e = grid::Te_allcells[nonemptymgi];
  const float clumpednne = grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi);

  const int nlevels_ionising = get_nlevels_ionising(element, ion);

  double ionisation_rate_rad = 0.;
  double ionisation_rate_coll = 0.;
  for (int level = 0; level < nlevels_ionising; level++) {
    const double nnlevel = calculate_levelpop(nonemptymgi, element, ion, level);
    const int nphixstargets = get_nphixstargets(element, ion, level);
    for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
      const int upperlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);

      ionisation_rate_rad += nnlevel * get_corrphotoioncoeff(element, ion, level, phixstargetindex, nonemptymgi, false);

      const double epsilon_trans = epsilon(element, ion + 1, upperlevel) - epsilon(element, ion, level);

      ionisation_rate_coll +=
          nnlevel * col_ionisation_ratecoeff(T_e, clumpednne, element, ion, level, phixstargetindex, epsilon_trans);
    }
  }
  const auto ionisation_rate = (ionisation_rate_rad + ionisation_rate_coll);
  const auto groundlevelpop = get_groundlevelpop(nonemptymgi, element, ion);
  // groundlevelpop is exactly zero only for an absent element (massfrac == 0), where ionisation_rate is
  // also zero. Avoid the resulting 0/0 = NaN being stored in the gamma estimator.
  if (groundlevelpop <= 0.) {
    return 0.;
  }
  return ionisation_rate / groundlevelpop;
}

// ionisation rate coefficient. multiply by the lower ion pop to get a rate.
// Currently only used for the estimator output file, not the simulation
auto calculate_iongamma_per_ionpop(const int nonemptymgi, const int element, const int lowerion,
                                   const bool collisional_not_radiative, const bool force_bfintegral) -> double {
  assert_always(lowerion < get_nions(element) - 1);
  // this option only makes sense for radiative ionisation
  assert_always(!collisional_not_radiative || (!force_bfintegral));

  const auto nnlowerion = get_nnion(nonemptymgi, element, lowerion);
  if (nnlowerion <= 0.) {
    return 0.;
  }

  const auto clumpednne = grid::get_clumpfactor(nonemptymgi) * grid::get_nne(nonemptymgi);
  const auto T_e = grid::Te_allcells[nonemptymgi];

  double ionisation_rate = 0.;  // rate per second
  const auto nlevels_ionising = get_nlevels_ionising(element, lowerion);
  for (int lower = 0; lower < nlevels_ionising; lower++) {
    const auto nnlowerlevel = calculate_levelpop(nonemptymgi, element, lowerion, lower);
    const auto nphixstargets = get_nphixstargets(element, lowerion, lower);
    for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
      if (collisional_not_radiative) {
        const int upper = get_phixsupperlevel(element, lowerion, lower, phixstargetindex);
        const double epsilon_trans = epsilon(element, lowerion + 1, upper) - epsilon(element, lowerion, lower);
        ionisation_rate += nnlowerlevel * col_ionisation_ratecoeff(T_e, clumpednne, element, lowerion, lower,
                                                                   phixstargetindex, epsilon_trans);
      } else if (force_bfintegral) {
        // don't use any detailed bound-free estimators, even if they are on and available
        ionisation_rate += nnlowerlevel * calculate_corrphotoioncoeff_integral(element, lowerion, lower,
                                                                               phixstargetindex, nonemptymgi, false);
      } else {
        // whatever ARTIS uses internally, maybe using detailed bound-free estimators
        ionisation_rate +=
            nnlowerlevel * get_corrphotoioncoeff(element, lowerion, lower, phixstargetindex, nonemptymgi, false);
      }
    }
  }

  return ionisation_rate / nnlowerion;
}
