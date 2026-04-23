#include "ratecoeff.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <ios>
#include <iterator>
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
#include "ltepop.h"
#include "macroatom.h"
#include "mpi_logging.h"
#include "radfield.h"
#include "random.h"

namespace {
constexpr double RATECOEFF_INTEGRAL_ACCURACY = 1e-3;

const auto temperature_grid = []() {
  const double T_step_log = (std::log(MAXTEMP) - std::log(MINTEMP)) / (TABLESIZE - 1.);
  std::array<double, TABLESIZE + 1> grid{};
  for (auto i = 0UZ; i < grid.size(); i++) {
    grid[i] = MINTEMP * std::exp(i * T_step_log);
  }
  return grid;
}();

MPI_shared_array<const float> ion_alpha_sp;  // size is nincludedions * TABLESIZE
                                             //

// the following spans are indexed by get_bflutindex()
MPI_shared_array<double> spontrecombcoeffs{};  // indexed by get_bflutindex()
MPI_shared_array<double> corrphotoioncoeffs{};  // for USE_LUT_PHOTOION = true
MPI_shared_array<double> bfcooling_coeffs{};

// Integrand to calculate the rate coefficient for spontaneous recombination
auto alpha_sp_integrand(const double nu_minus_nu_edge, const double nu_edge, const float T_e,
                        const std::span<const float> photoion_xs) -> double {
  const auto sigma_bf = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_minus_nu_edge + nu_edge);
  // the variable of integration has been changed from nu to nu_edge_minus_nu = nu - nu_edge
  // to get a cancellation with part of the saha factor
  return TWOOVERCLIGHTSQUARED * sigma_bf * pow2(nu_edge + nu_minus_nu_edge) * exp(-HOVERKB * nu_minus_nu_edge / T_e);
}

// Integrand to calculate the rate coefficient for spontaneous recombination
auto alpha_sp_E_integrand(const double nu, const double nu_edge, const float T_e,
                          const std::span<const float> photoion_xs) -> double {
  const auto sigma_bf = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu);
  return TWOOVERCLIGHTSQUARED * sigma_bf * pow3(nu) / nu_edge * exp(-HOVERKB * nu / T_e);
}

// Integrand to calculate the rate coefficient for photoionisation corrected for stimulated recombination.
auto gammacorr_integrand(const double nu, const double nu_edge, const float T, const std::span<const float> photoion_xs)
    -> double {
  const auto sigma_bf = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu);

  // The correction factor for stimulated emission in gammacorr is set to its
  // LTE value. Because the T_e dependence of gammacorr is weak, this correction
  // correction may be evaluated at T_R!

  // Dependence on dilution factor W is linear. This allows to set it here to
  // 1. and scale to its actual value later on.
  // Assumption T_e = T_R makes n_kappa/n_i * (n_i/n_kappa)* = 1
  return sigma_bf * ONEOVERH / nu * radfield::dbb(nu, T, 1) * (1 - exp(-HOVERKB * nu / T));
}

// Integrand to precalculate the bound-free cooling rate coefficient
auto bfcooling_integrand(const double nu_minus_nu_edge, const double nu_edge, const float T_e,
                         const std::span<const float> photoion_xs) -> double {
  const float sigma_bf = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_minus_nu_edge + nu_edge);

  // return sigma_bf * (1-nu_edge/nu) * TWOHOVERCLIGHTSQUARED * pow(nu,3) * exp(-HOVERKB*nu/T);
  return sigma_bf * nu_minus_nu_edge * TWOHOVERCLIGHTSQUARED * (nu_minus_nu_edge + nu_edge) *
         (nu_minus_nu_edge + nu_edge) * exp(-HOVERKB * nu_minus_nu_edge / T_e);
}

[[gnu::pure]] [[nodiscard]] inline auto get_bflutindex(const int temperatureindex, const int uniquelevelindex,
                                                       const int phixstargetindex) -> int {
  const int contindex = globals::alllevels.bflist_start[uniquelevelindex] + phixstargetindex;
  const int bflutindex = (temperatureindex * globals::nbfcontinua) + contindex;
  assert_testmodeonly(bflutindex >= 0);
  assert_testmodeonly(bflutindex <= TABLESIZE * globals::nbfcontinua);
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
    printlog("Performing rate integrals for Z = {}: ion stages", atomic_number);
    const int nions = get_nions(element) - 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ion = 0; ion < nions; ion++) {
      const int ionstage = get_ionstage(element, ion);
      printlog(" {}", ionstage);
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

          // const double E_threshold = epsilon(element,ion+1,upperlevel) - epsilon(element,ion,level);
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
              if (gammacorr < 0) {
                printlnlog("WARNING: gammacorr was negative for level {}", level);
                gammacorr = 0;
              }
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
    printlnlog("");
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

    auto phixstable = globals::allphixs.subspan(phixsstart * globals::NPHIXSPOINTS, globals::NPHIXSPOINTS);
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

  printlnlog("Calibrating recombination rates for a temperature of {:.1f} K", Te_estimate);

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
    std::stringstream(line) >> atomicnumber >> upperionstage >> tablerows;
    RRCRow T_highestbelow = {.log_Te = 0, .rrc_low_n = 0, .rrc_total = 0};
    RRCRow T_lowestabove = {.log_Te = 0, .rrc_low_n = 0, .rrc_total = 0};
    T_highestbelow.log_Te = -1;
    T_lowestabove.log_Te = -1;
    for (int i = 0; i < tablerows; i++) {
      RRCRow row{};
      get_noncommentline(recombrate_file, line);
      std::stringstream(line) >> row.log_Te >> row.rrc_low_n >> row.rrc_total;

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

        const double x = (log_Te_estimate - T_highestbelow.log_Te) / (T_lowestabove.log_Te - T_highestbelow.log_Te);
        const double input_rrc_low_n = (x * T_highestbelow.rrc_low_n) + ((1 - x) * T_lowestabove.rrc_low_n);
        const double input_rrc_total = (x * T_highestbelow.rrc_total) + ((1 - x) * T_lowestabove.rrc_total);

        constexpr bool assume_lte = true;
        constexpr bool per_groundmultipletpop = true;
        constexpr bool collisional_not_radiative = false;

        double rrc = calculate_ionrecombcoeff(-1, Te_estimate, element, ion, assume_lte, collisional_not_radiative,
                                              false, per_groundmultipletpop);
        printlnlog("              rrc: {:10.3e}", rrc);

        if (input_rrc_low_n >= 0)  // if it's < 0, ignore it
        {
          printlnlog("  input_rrc_low_n: {:10.3e}", input_rrc_low_n);

          const double phixs_multiplier = input_rrc_low_n / rrc;
          if (phixs_multiplier < 0.05 || phixs_multiplier >= 2.0) {
            printlnlog("    Not scaling phixs of all levels by {:.3f} (because < 0.05 or >= 2.0)", phixs_multiplier);
          } else {
            printlnlog("    scaling phixs of all levels by {:.3f}", phixs_multiplier);

            for (int level = 0; level < nlevels; level++) {
              scale_level_phixs(element, ion - 1, level, phixs_multiplier);
            }

            rrc = calculate_ionrecombcoeff(-1, Te_estimate, element, ion, assume_lte, collisional_not_radiative, false,
                                           per_groundmultipletpop);
            printlnlog("              rrc: {:10.3e}", rrc);
          }
        }

        // hopefully the RRC now matches the low_n value well, if it was defined
        // Next, use the superlevel recombination rates to make up the excess needed to reach the total RRC

        printlnlog("  input_rrc_total: {:10.3e}", input_rrc_total);

        if (rrc < input_rrc_total) {
          const double rrc_superlevel = calculate_ionrecombcoeff(
              -1, Te_estimate, element, ion, assume_lte, collisional_not_radiative, true, per_groundmultipletpop);
          printlnlog("  rrc(superlevel): {:10.3e}", rrc_superlevel);

          if (rrc_superlevel > 0) {
            const double phixs_multiplier_superlevel = 1.0 + ((input_rrc_total - rrc) / rrc_superlevel);
            printlnlog("    scaling phixs of levels in the superlevel by {:.3f}", phixs_multiplier_superlevel);
            assert_always(phixs_multiplier_superlevel >= 0);

            const int first_superlevel_level = get_nlevels_excited_nlte(element, ion - 1) + 1;
            for (int level = first_superlevel_level; level < nlevels; level++) {
              scale_level_phixs(element, ion - 1, level, phixs_multiplier_superlevel);
            }
          } else {
            printlnlog("There is no superlevel recombination, so multiplying all levels instead");
            const double phixs_multiplier = input_rrc_total / rrc;
            printlnlog("    scaling phixs of all levels by {:.3f}", phixs_multiplier);
            assert_always(phixs_multiplier >= 0);

            for (int level = 0; level < nlevels; level++) {
              scale_level_phixs(element, ion - 1, level, phixs_multiplier);
            }
          }
        } else {
          printlnlog("rrc >= input_rrc_total!");
          const double phixs_multiplier = input_rrc_total / rrc;
          printlnlog("    scaling phixs of all levels by {:.3f}", phixs_multiplier);
          assert_always(phixs_multiplier >= 0);

          for (int level = 0; level < nlevels; level++) {
            scale_level_phixs(element, ion - 1, level, phixs_multiplier);
          }
        }

        rrc = calculate_ionrecombcoeff(-1, Te_estimate, element, ion, assume_lte, collisional_not_radiative, false,
                                       per_groundmultipletpop);
        printlnlog("              rrc: {:10.3e}", rrc);
      }
    }
  }
}

void precalculate_ion_alpha_sp() {
  auto temp_ion_alpha_sp = MPI_shared_array<float>(get_includedions() * TABLESIZE, 0.);
  if (globals::rank_in_node == 0) {
    const auto nincludedions = get_includedions();
    for (int tempindex = 0; tempindex < TABLESIZE; tempindex++) {
      const auto T_e = static_cast<float>(temperature_grid[tempindex]);
      for (int element = 0; element < get_nelements(); element++) {
        const int nions = get_nions(element) - 1;
        for (int ion = 0; ion < nions; ion++) {
          const auto uniqueionindex = get_uniqueionindex(element, ion);
          const int nionisinglevels = get_nlevels_ionising(element, ion);
          double zeta = 0.;
          for (int level = 0; level < nionisinglevels; level++) {
            const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
            const auto nphixstargets = get_nphixstargets(uniquelevelindex);
            for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
              const double zeta_level = get_spontrecombcoeff(uniquelevelindex, phixstargetindex, T_e);
              zeta += zeta_level;
            }
          }
          temp_ion_alpha_sp[(tempindex * nincludedions) + uniqueionindex] = static_cast<float>(zeta);
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
  return ONEOVERH * sigma_bf / (nu_minus_nu_edge + nu_edge) * Jnu * corrfactor;
}

auto calculate_corrphotoioncoeff_integral(const int element, const int ion, const int level, const int phixstargetindex,
                                          const int nonemptymgi, const bool use_cellcache) -> double {
  constexpr double epsrel = 1e-3;

  const auto loweruniquelevelindex = get_uniquelevelindex(element, ion, level);
  const double nu_threshold = ONEOVERH * get_phixs_threshold(loweruniquelevelindex, phixstargetindex);
  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

  const auto T_e = grid::get_Te(nonemptymgi);

  // stimulated recombination is negative photoionisation
  const double nnlevel = use_cellcache ? get_cellcache_levelpop(nonemptymgi, loweruniquelevelindex)
                                       : calculate_levelpop(nonemptymgi, element, ion, level);
  const double nne = grid::get_nne(nonemptymgi);
  const int upperionlevel = get_phixsupperlevel(loweruniquelevelindex, phixstargetindex);
  const auto upperuniquelevelindex = get_uniquelevelindex(element, ion + 1, upperionlevel);
  const double modified_sahafact =
      SAHACONST * stat_weight(loweruniquelevelindex) / stat_weight(upperuniquelevelindex) * std::pow(T_e, -1.5);
  const double nnupperionlevel = use_cellcache ? get_cellcache_levelpop(nonemptymgi, upperuniquelevelindex)
                                               : calculate_levelpop(nonemptymgi, element, ion + 1, upperionlevel);

  double modified_departure_ratio = nnlevel > 0. ? nnupperionlevel / nnlevel * nne * modified_sahafact : 1.;
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

template <typename T>
[[nodiscard]] DEVICE_FUNC auto lerp_or_last(const std::span<T> table, const int uniquelevelindex,
                                            const int phixstargetindex, const auto temperature) {
  const auto upperindex =
      std::ranges::distance(temperature_grid.begin(), std::ranges::upper_bound(temperature_grid, temperature));
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
  size_t mem_usage_photoionluts = 2 * TABLESIZE * globals::nbfcontinua * sizeof(double);

  spontrecombcoeffs = MPI_shared_array<double>(TABLESIZE * globals::nbfcontinua, 0.);

  if constexpr (USE_LUT_PHOTOION) {
    corrphotoioncoeffs = MPI_shared_array<double>(TABLESIZE * globals::nbfcontinua, 0.);
    mem_usage_photoionluts += TABLESIZE * globals::nbfcontinua * sizeof(double);
  }

  bfcooling_coeffs = MPI_shared_array<double>(TABLESIZE * globals::nbfcontinua, 0.);

  printlnlog(
      "[info] mem_usage: lookup tables derived from photoionisation (spontrecombcoeff, bfcooling and "
      "corrphotoioncoeff if USE_LUT_PHOTOION) occupy {:.3f} MB",
      mem_usage_photoionluts / 1024. / 1024.);
}

DEVICE_FUNC auto select_continuum_nu(int element, const int lowerion, const int lower, const int upperionlevel,
                                     float T_e) -> double {
  const auto lower_uniquelevelindex = get_uniquelevelindex(element, lowerion, lower);
  const int phixstargetindex = get_phixtargetindex(lower_uniquelevelindex, upperionlevel);
  const double E_threshold = get_phixs_threshold(lower_uniquelevelindex, phixstargetindex);
  const double nu_threshold = ONEOVERH * E_threshold;

  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

  const int npieces = globals::NPHIXSPOINTS;

  const auto photoion_xs = get_phixs_table(lower_uniquelevelindex);

  const double zrand = 1. - rng_uniform();  // Make sure that 0 < zrand <= 1

  const double deltanu = (nu_max_phixs - nu_threshold) / npieces;
  double error{NAN};

  auto total_alpha_sp =
      integrator<31>([&](const double nu) { return alpha_sp_E_integrand(nu, nu_threshold, T_e, photoion_xs); },
                     nu_threshold, nu_max_phixs, RATECOEFF_INTEGRAL_ACCURACY, &error);

  double alpha_sp_old = total_alpha_sp;
  double alpha_sp = total_alpha_sp;

  int i = 1;
  for (i = 1; i < npieces; i++) {
    alpha_sp_old = alpha_sp;
    const double xlow = nu_threshold + (i * deltanu);

    // Spontaneous recombination and bf-cooling coefficient don't depend on the radiation field
    alpha_sp = integrator<31>([&](const double nu) { return alpha_sp_E_integrand(nu, nu_threshold, T_e, photoion_xs); },
                              xlow, nu_max_phixs, RATECOEFF_INTEGRAL_ACCURACY, &error);

    if (zrand >= alpha_sp / total_alpha_sp) {
      break;
    }
  }

  const double nuoffset =
      (alpha_sp != alpha_sp_old) ? ((total_alpha_sp * zrand) - alpha_sp_old) / (alpha_sp - alpha_sp_old) * deltanu : 0.;
  const double nu_lower = nu_threshold + ((i - 1) * deltanu) + nuoffset;

  assert_testmodeonly(std::isfinite(nu_lower));

  return nu_lower;
}

// Get an ion's rate coefficient for spontaneous recombination in LTE
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_ion_spontrecombcoeff(const int uniqueionindex, const float T_e)
    -> double {
  const auto upperindex =
      std::ranges::distance(temperature_grid.begin(), std::ranges::upper_bound(temperature_grid, T_e));
  const auto nincludedions = get_includedions();
  if (upperindex < TABLESIZE) {
    const double T_lower = temperature_grid[upperindex - 1];
    const double T_upper = temperature_grid[upperindex];

    const double f_lower = ion_alpha_sp[((upperindex - 1) * nincludedions) + uniqueionindex];
    const double f_upper = ion_alpha_sp[(upperindex * nincludedions) + uniqueionindex];

    return f_lower + ((f_upper - f_lower) / (T_upper - T_lower) * (T_e - T_lower));
  }
  return ion_alpha_sp[((TABLESIZE - 1) * nincludedions) + uniqueionindex];
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
    // const int nphixstargets = get_nphixstargets(element, lowerion, 0);
    // upper_nlevels = get_phixsupperlevel(element, lowerion, 0, nphixstargets - 1) + 1;

    upper_nlevels = get_nlevels_groundterm(element, lowerion + 1);
  } else {
    upper_nlevels = get_nlevels(element, lowerion + 1);
  }

  for (int upper = 0; upper < upper_nlevels; upper++) {
    double nnupperlevel{NAN};
    if (assume_lte) {
      const double T_exc = T_e;
      const double E_level = epsilon(element, lowerion + 1, upper);
      const double E_ground = epsilon(element, lowerion + 1, 0);
      const double nnground = (nonemptymgi >= 0) ? get_groundlevelpop(nonemptymgi, element, lowerion + 1) : 1.;

      nnupperlevel = (nnground * stat_weight(element, lowerion + 1, upper) / stat_weight(element, lowerion + 1, 0) *
                      exp(-(E_level - E_ground) / KB / T_exc));
    } else {
      nnupperlevel = calculate_levelpop(nonemptymgi, element, lowerion + 1, upper);
    }
    nnupperion += nnupperlevel;
  }

  if (nnupperion <= 0.) {
    return 0.;
  }

  // this gets divided and cancelled out in the radiative case anyway
  const auto nne = (nonemptymgi >= 0) ? grid::get_nne(nonemptymgi) : 1.F;
  double alpha = 0.;
  const int maxrecombininglevel = get_maxrecombininglevel(element, lowerion + 1);
  for (int upper = 0; upper <= maxrecombininglevel; upper++) {
    double nnupperlevel{NAN};
    if (assume_lte) {
      const double T_exc = T_e;
      const double E_level = epsilon(element, lowerion + 1, upper);
      const double E_ground = epsilon(element, lowerion + 1, 0);
      const double nnground = (nonemptymgi >= 0) ? get_groundlevelpop(nonemptymgi, element, lowerion + 1) : 1.;

      nnupperlevel = (nnground * stat_weight(element, lowerion + 1, upper) / stat_weight(element, lowerion + 1, 0) *
                      exp(-(E_level - E_ground) / KB / T_exc));
    } else {
      nnupperlevel = calculate_levelpop(nonemptymgi, element, lowerion + 1, upper);
    }
    for (int lower = 0; lower < get_nlevels(element, lowerion); lower++) {
      if (lower_superlevel_only && (!level_isinsuperlevel(element, lowerion, lower))) {
        continue;
      }

      double recomb_coeff = 0.;
      if (collisional_not_radiative) {
        const double epsilon_trans = epsilon(element, lowerion + 1, upper) - epsilon(element, lowerion, lower);
        recomb_coeff += col_recombination_ratecoeff(T_e, nne, element, upperion, upper, lower, epsilon_trans);
      } else {
        recomb_coeff += rad_recombination_ratecoeff(T_e, nne, element, lowerion + 1, upper, lower);
      }

      const double alpha_level = recomb_coeff / nne;
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
auto get_corrphotoioncoeff_ana(int element, const int ion, const int level, const int phixstargetindex,
                               const int nonemptymgi) -> double {
  assert_always(USE_LUT_PHOTOION);
  const double W = grid::get_W(nonemptymgi);
  const double T_R = grid::get_TR(nonemptymgi);
  const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
  return W * lerp_or_last(std::span{corrphotoioncoeffs}, uniquelevelindex, phixstargetindex, T_R);
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
  const auto allphisxtargetindex = get_allphixstargetindex(uniquelevelindex, phixstargetindex);
  double gammacorr =
      use_cellcache ? globals::cellcache[cellcacheslotid].allphixstargets_corrphotoioncoeff[allphisxtargetindex] : -1;

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
        const double W = grid::get_W(nonemptymgi);
        const double T_R = grid::get_TR(nonemptymgi);

        gammacorr = W * lerp_or_last(std::span{corrphotoioncoeffs}, uniquelevelindex, phixstargetindex, T_R);
        const int index_in_groundlevelcontestimator = globals::alllevels.closestgroundlevelcont[uniquelevelindex];
        if (index_in_groundlevelcontestimator >= 0) {
          gammacorr *= globals::corrphotoionrenorm[(nonemptymgi * globals::nbfcontinua_ground) +
                                                   index_in_groundlevelcontestimator];
        }
      }
    }
    if (use_cellcache) {
      globals::cellcache[cellcacheslotid].allphixstargets_corrphotoioncoeff[allphisxtargetindex] = gammacorr;
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
    return (globals::gammaestimator[(nonemptymgi * globals::nbfcontinua_ground) + groundcontindex] == 0);
  }

  const auto T_e = grid::get_Te(nonemptymgi);
  const auto nne = grid::get_nne(nonemptymgi);

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

      if (nnlevel * col_ionisation_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans) > 0) {
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

  const auto T_e = grid::get_Te(nonemptymgi);
  const float nne = grid::get_nne(nonemptymgi);

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
          nnlevel * col_ionisation_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);
    }
  }
  const auto ionisation_rate = (ionisation_rate_rad + ionisation_rate_coll);
  return ionisation_rate / get_groundlevelpop(nonemptymgi, element, ion);
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

  const auto nne = grid::get_nne(nonemptymgi);
  const auto T_e = grid::get_Te(nonemptymgi);

  double ionisation_rate = 0.;  // rate per second
  const auto nlevels_ionising = get_nlevels_ionising(element, lowerion);
  for (int lower = 0; lower < nlevels_ionising; lower++) {
    const auto nnlowerlevel = calculate_levelpop(nonemptymgi, element, lowerion, lower);
    const auto nphixstargets = get_nphixstargets(element, lowerion, lower);
    for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
      if (collisional_not_radiative) {
        const int upper = get_phixsupperlevel(element, lowerion, lower, phixstargetindex);
        const double epsilon_trans = epsilon(element, lowerion + 1, upper) - epsilon(element, lowerion, lower);
        ionisation_rate += nnlowerlevel * col_ionisation_ratecoeff(T_e, nne, element, lowerion, lower, phixstargetindex,
                                                                   epsilon_trans);
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
