#include "thermalbalance.h"

#pragma clang unsafe_buffer_usage begin
#if defined(USE_BOOST) && USE_BOOST
#include <boost/assert/source_location.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#else
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#endif
#pragma clang unsafe_buffer_usage end

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <span>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "sn3d.h"

namespace {

struct TeSolutionParams {
  double t_current;
  int nonemptymgi;
  HeatingCoolingRates* heatingcoolingrates;
  const std::vector<double>* bfheatingcoeffs;
};

struct BFHeatingIntegralParams {
  double nu_edge;
  int nonemptymgi;
  float T_R;
  std::span<const float> photoion_xs;
};

// Integrand to calculate the rate coefficient for bfheating.
auto integrand_bfheatingcoeff_custom_radfield(const double nu, void* const voidparas) -> double {
  const auto* const params = static_cast<const BFHeatingIntegralParams*>(voidparas);

  const int nonemptymgi = params->nonemptymgi;
  const double nu_edge = params->nu_edge;
  const float T_R = params->T_R;

  const float sigma_bf = photoionisation_crosssection_fromtable(params->photoion_xs, nu_edge, nu);

  return sigma_bf * (1 - (nu_edge / nu)) * radfield::radfield(nu, nonemptymgi) * (1 - exp(-HOVERKB * nu / T_R));
}

auto calculate_bfheatingcoeff(const int element, const int ion, const int level, const int phixstargetindex,
                              const int nonemptymgi) -> double {
  double error = 0.;
  const double epsrel = 1e-3;

  const double E_threshold = get_phixs_threshold(element, ion, level, phixstargetindex);

  const double nu_threshold = ONEOVERH * E_threshold;
  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

  const BFHeatingIntegralParams intparas = {.nu_edge = nu_threshold,
                                            .nonemptymgi = nonemptymgi,
                                            .T_R = grid::get_TR(nonemptymgi),
                                            .photoion_xs = get_phixs_table(element, ion, level)};

  double bfheating = 0.;

  const int status = integrator<integrand_bfheatingcoeff_custom_radfield>(intparas, nu_threshold, nu_max_phixs, epsrel,
                                                                          &bfheating, &error);

  if (status != 0 && (status != 18 || (error / bfheating) > 0.1)) {
#if !USE_SIMPSON_INTEGRATOR
    printlnlog(
        "bf_heating integrator gsl warning {}. modelgridindex {} Z={} ionstage {} lower {} phixstargetindex {} "
        "integral {:g} error {:g}",
        status, grid::get_mgi_of_nonemptymgi(nonemptymgi), get_atomicnumber(element), get_ionstage(element, ion), level,
        phixstargetindex, bfheating, error);
#endif
  }

  bfheating *= FOURPI * get_phixsprobability(element, ion, level, phixstargetindex);

  return bfheating;
}

auto get_heating_ion_coll_deexc(const int nonemptymgi, const int element, const int ion, const float T_e,
                                const float nne) -> double {
  double C_deexc = 0.;
  const int nlevels = get_nlevels(element, ion);
  const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
  for (int level = 0; level < nlevels; level++) {
    const auto uniquelevelindex = ionuniquelevelindexstart + level;
    const double nnlevel = calculate_levelpop(nonemptymgi, element, ion, level);
    const double epsilon_level = epsilon(uniquelevelindex);
    const double statweight = stat_weight(uniquelevelindex);

    // Collisional heating: deexcitation to same ionisation stage
    const auto alltrans_startdown = get_alltrans_startdown(uniquelevelindex);
    const int ndowntrans = get_ndowntrans(uniquelevelindex);
    for (int i = 0; i < ndowntrans; i++) {
      const auto alltransindex = alltrans_startdown + i;
      const int lower = globals::alltrans.targetlevelindex[alltransindex];
      const double epsilon_trans = epsilon_level - epsilon(ionuniquelevelindexstart + lower);
      const auto lower_statweight = stat_weight(ionuniquelevelindexstart + lower);
      const double C =
          nnlevel * col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, statweight, lower_statweight, alltransindex) *
          epsilon_trans;
      C_deexc += C;
    }
  }
  return C_deexc;
}

// Calculate the heating rates for a given cell. Results are returned via the elements of the heatingrates data
// structure.
void calculate_heating_rates(const int nonemptymgi, const float T_e, const float nne,
                             HeatingCoolingRates& heatingcoolingrates, const std::vector<double>& bfheatingcoeffs) {
  double C_deexc = 0.;

  // double C_recomb = 0.;
  double bfheating = 0.;
  double ffheating = 0.;

  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    if constexpr (DIRECT_COL_HEAT) {
      for (int ion = 0; ion < nions; ion++) {
        C_deexc += get_heating_ion_coll_deexc(nonemptymgi, element, ion, T_e, nne);
      }
    }

    // Collisional heating: recombination to lower ionisation stage (not included)

    // Bound-free heating (renormalised analytical calculation)
    for (int ion = 0; ion < nions - 1; ion++) {
      const int nbflevels = get_nlevels_ionising(element, ion);
      const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
      for (int level = 0; level < nbflevels; level++) {
        const auto uniquelevelindex = ionuniquelevelindexstart + level;
        const double nnlevel = calculate_levelpop(nonemptymgi, element, ion, level);
        bfheating += nnlevel * bfheatingcoeffs[uniquelevelindex];
      }
    }
  }

  // Free-free heating (from estimators)
  ffheating = globals::ffheatingestimator[nonemptymgi];

  if constexpr (DIRECT_COL_HEAT) {
    heatingcoolingrates.heating_collisional = C_deexc;
  } else {
    // Collisional heating (from estimators)
    heatingcoolingrates.heating_collisional = globals::colheatingestimator.at(nonemptymgi);  // C_deexc + C_recomb;
  }

  heatingcoolingrates.heating_bf = bfheating;
  heatingcoolingrates.heating_ff = ffheating;
}

// Thermal balance equation on which we have to iterate to get T_e

auto T_e_eqn_heating_minus_cooling(const double T_e, int nonemptymgi, const double t_current,
                                   HeatingCoolingRates& heatingcoolingrates, const std::vector<double>& bfheatingcoeffs)
    -> double {
  const auto fT_e = static_cast<float>(T_e);

  if constexpr (!LTEPOP_EXCITATION_USE_TJ) {
    if (std::abs((T_e / grid::get_Te(nonemptymgi)) - 1.) > 0.1) {
      grid::set_Te(nonemptymgi, fT_e);
      for (int element = 0; element < get_nelements(); element++) {
        if (!elem_has_nlte_levels(element)) {
          // recalculate the Gammas using the current level populations
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions - 1; ion++) {
            if (get_groundcontindex(element, ion) >= 0) {
              globals::gammaestimator[get_ionestimindex_nonemptymgi(nonemptymgi, element, ion)] =
                  calculate_iongamma_per_gspop(nonemptymgi, element, ion);
            }
          }
        }
      }
    }
  }

  // Set new T_e guess for the current cell and update populations
  grid::set_Te(nonemptymgi, fT_e);

  calculate_ion_balance_nne(nonemptymgi);
  const auto nne = grid::get_nne(nonemptymgi);

  // Then calculate heating and cooling rates
  kpkt::calculate_cooling_rates(nonemptymgi, &heatingcoolingrates);
  calculate_heating_rates(nonemptymgi, fT_e, nne, heatingcoolingrates, bfheatingcoeffs);

  const auto ntlepton_frac_heating = nonthermal::get_nt_frac_heating(nonemptymgi);
  const auto ntlepton_dep = nonthermal::get_deposition_rate_density(nonemptymgi);
  const auto ntalpha_frac_heating = 1.;
  const auto ntalpha_dep = heatingcoolingrates.dep_alpha;
  heatingcoolingrates.heating_dep = (ntlepton_dep * ntlepton_frac_heating) + (ntalpha_dep * ntalpha_frac_heating);
  heatingcoolingrates.dep_frac_heating =
      (ntalpha_dep > 0) ? heatingcoolingrates.heating_dep / (ntlepton_dep + ntalpha_dep) : ntlepton_frac_heating;

  // Adiabatic cooling term
  const double nntot = get_nnion_tot(nonemptymgi) + nne;
  const double p = nntot * KB * T_e;  // pressure in [erg/cm^3]
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double volumetmin = grid::get_modelcell_assocvolume_tmin(modelgridindex);
  const double dV_on_dt = 3 * volumetmin / pow(globals::tmin, 3) * pow(t_current, 2);
  const double V = volumetmin * pow(t_current / globals::tmin, 3);
  heatingcoolingrates.cooling_adiabatic = p * dV_on_dt / V;

  const double total_heating_rate = heatingcoolingrates.heating_ff + heatingcoolingrates.heating_bf +
                                    heatingcoolingrates.heating_collisional + heatingcoolingrates.heating_dep;
  const double total_coolingrate = heatingcoolingrates.cooling_ff + heatingcoolingrates.cooling_fb +
                                   heatingcoolingrates.cooling_collisional + heatingcoolingrates.cooling_adiabatic;

  return total_heating_rate - total_coolingrate;
}
auto T_e_eqn_heating_minus_cooling(const double T_e, void* const paras)  // cppcheck-suppress constParameterPointer
    -> double {
  const auto* const params = static_cast<const TeSolutionParams*>(paras);
  return T_e_eqn_heating_minus_cooling(T_e, params->nonemptymgi, params->t_current, *params->heatingcoolingrates,
                                       *params->bfheatingcoeffs);
}

}  // anonymous namespace

// depends only the radiation field - no dependence on T_e or populations
void calculate_bfheatingcoeffs(int nonemptymgi, std::vector<double>& bfheatingcoeffs) {
  bfheatingcoeffs.resize(get_includedlevels());
  const double minelfrac = 0.01;
  for (int element = 0; element < get_nelements(); element++) {
    if (grid::get_elem_abundance(nonemptymgi, element) <= minelfrac && !USE_LUT_BFHEATING) {
      printlog("skipping Z={} X={:g}, ", get_atomicnumber(element), grid::get_elem_abundance(nonemptymgi, element));
    }

    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);
      const auto levels = std::ranges::iota_view{0, nlevels};
      std::for_each(EXEC_PAR levels.begin(), levels.end(), [&](const int level) {
        double bfheatingcoeff = 0.;
        if (grid::get_elem_abundance(nonemptymgi, element) > minelfrac || USE_LUT_BFHEATING) {
          const auto nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            if constexpr (!USE_LUT_BFHEATING) {
              bfheatingcoeff += calculate_bfheatingcoeff(element, ion, level, phixstargetindex, nonemptymgi);
            } else {
              // The correction factor for stimulated emission in gammacorr is set to its
              // LTE value. Because the T_e dependence of gammacorr is weak, this correction
              // correction may be evaluated at T_R!
              const double T_R = grid::get_TR(nonemptymgi);
              const double W = grid::get_W(nonemptymgi);
              bfheatingcoeff += get_bfheatingcoeff_ana(element, ion, level, phixstargetindex, T_R, W);
            }
          }
          assert_always(std::isfinite(bfheatingcoeff));

          if constexpr (USE_LUT_BFHEATING) {
            const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
            const int index_in_groundlevelcontestimator = globals::alllevels.closestgroundlevelcont[uniquelevelindex];
            if (index_in_groundlevelcontestimator >= 0) {
              bfheatingcoeff *= globals::bfheatingestimator[(nonemptymgi * globals::nbfcontinua_ground) +
                                                            index_in_groundlevelcontestimator];
            }
          }
        }
        bfheatingcoeffs[get_uniquelevelindex(element, ion, level)] = bfheatingcoeff;
      });
    }
  }
}

void call_T_e_finder(const int nonemptymgi, const double t_current, const double T_min, const double T_max,
                     HeatingCoolingRates& heatingcoolingrates, const std::vector<double>& bfheatingcoeffs) {
  const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double T_e_old = grid::get_Te(nonemptymgi);
#if defined(USE_BOOST) && USE_BOOST
  constexpr auto method = "boost toms748_solve";
#else
  constexpr auto method = "gsl brent";
#endif
  printlog("Finding T_e in cell {} at timestep {} [{}]...", modelgridindex, globals::timestep, method);

  const auto f_T_e = [&](double T_e) -> double {
    return T_e_eqn_heating_minus_cooling(T_e, nonemptymgi, t_current, heatingcoolingrates, bfheatingcoeffs);
  };

  double thermalmin = f_T_e(T_min);
  double thermalmax = f_T_e(T_max);

  if (!std::isfinite(thermalmin) || !std::isfinite(thermalmax)) {
    printlnlog(
        "[abort request] call_T_e_finder: non-finite results in modelcell {} (T_R={:g}, W={:g}). T_e forced to be "
        "MINTEMP",
        modelgridindex, grid::get_TR(nonemptymgi), grid::get_W(nonemptymgi));
    thermalmax = thermalmin = -1;
  }

  double T_e{NAN};
  // Check whether the thermal balance equation has a root in [T_min,T_max]
  if (thermalmin * thermalmax < 0) {
    const auto maxit = 100U;
    // If it has, then solve for the root T_e
#if defined(USE_BOOST) && USE_BOOST
    {
      // use TOMS 748 solver from Boost
      boost::uintmax_t iternum = maxit;
      auto result = boost::math::tools::toms748_solve(f_T_e, T_min, T_max, ftol<TEMPERATURE_SOLVER_ACCURACY>, iternum);
      T_e = 0.5 * (result.first + result.second);
      if (iternum >= maxit) {
        printlnlog("[warning] call_T_e_finder: T_e did not converge within {} iterations. interval [{:g}, {:g}]",
                   iternum, result.first, result.second);
      } else {
        printlnlog("after {} iterations, T_e = {:g} K, interval [{:g}, {:g}]", iternum, T_e, result.first,
                   result.second);
      }
    }
#else
    {
      TeSolutionParams paras = {.t_current = t_current,
                                .nonemptymgi = nonemptymgi,
                                .heatingcoolingrates = &heatingcoolingrates,
                                .bfheatingcoeffs = &bfheatingcoeffs};

      gsl_function find_T_e_f = {.function = &T_e_eqn_heating_minus_cooling, .params = &paras};
      // one-dimensional GSL Brent root solver, bracketing type
      gsl_root_fsolver* T_e_solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

      gsl_root_fsolver_set(T_e_solver, &find_T_e_f, T_min, T_max);
      int status = 0;
      for (auto iternum = 0U; iternum < maxit; iternum++) {
        gsl_root_fsolver_iterate(T_e_solver);
        T_e = gsl_root_fsolver_root(T_e_solver);
        const double T_e_min = gsl_root_fsolver_x_lower(T_e_solver);
        const double T_e_max = gsl_root_fsolver_x_upper(T_e_solver);
        status = gsl_root_test_interval(T_e_min, T_e_max, 0, TEMPERATURE_SOLVER_ACCURACY);
        // printlnlog("iter {}, T_e interval [{:g}, {:g}], guess {:g}, status {}", iternum, T_e_min, T_e_max, T_e,
        // status);
        if (status != GSL_CONTINUE) {
          printlnlog("after {} iterations, T_e = {:g} K, interval [{:g}, {:g}]", iternum + 1, T_e, T_e_min, T_e_max);
          break;
        }
      }

      if (status == GSL_CONTINUE) {
        printlnlog("[warning] call_T_e_finder: T_e did not converge within {} iterations", maxit);
      }

      gsl_root_fsolver_free(T_e_solver);
    }
#endif
  }

  // Quick solver style: works if we can assume that there is either one or no
  // solution on [MINTEMP.MAXTEMP] (check that by doing a plot of heating-cooling vs. T_e)
  else if (thermalmax < 0) {
    // Thermal balance equation always negative ===> T_e = T_min
    T_e = MINTEMP;
    printlnlog(
        "[warning] call_T_e_finder: cooling bigger than heating at lower T_e boundary {:g} in modelcell {} "
        "(T_R={:g},W={:g}). T_e forced to be MINTEMP",
        MINTEMP, modelgridindex, grid::get_TR(nonemptymgi), grid::get_W(nonemptymgi));
  } else {
    // Thermal balance equation always negative ===> T_e = T_max
    T_e = MAXTEMP;
    printlnlog(
        "[warning] call_T_e_finder: heating bigger than cooling over the whole T_e range [{:g},{:g}] in modelcell {} "
        "(T_R={:g},W={:g}). T_e forced to be MAXTEMP",
        MINTEMP, MAXTEMP, modelgridindex, grid::get_TR(nonemptymgi), grid::get_W(nonemptymgi));
  }

  if (T_e > 2 * T_e_old) {
    T_e = 2 * T_e_old;
    printlnlog("use T_e damping in cell {}", modelgridindex);
    T_e = std::min(T_e, MAXTEMP);
  } else if (T_e < 0.5 * T_e_old) {
    T_e = 0.5 * T_e_old;
    printlnlog("use T_e damping in cell {}", modelgridindex);
    T_e = std::max(T_e, MINTEMP);
  }

  grid::set_Te(nonemptymgi, static_cast<float>(T_e));

  // this call with make sure heating/cooling rates and populations are updated for the final T_e
  // in case T_e got modified after the T_e solver finished
  f_T_e(T_e);
}
