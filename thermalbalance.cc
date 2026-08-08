// The electron thermal balance: computes the heating rates (bound-free, free-free, collisional
// deexcitation, and non-thermal) and solves heating = cooling for the electron temperature T_e
// of each cell.

#include "thermalbalance.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <span>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <cstdint>

#include "toms748.h"
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "integrator.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "mpi_logging.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "sn3d.h"

namespace {

// Integrand to calculate the rate coefficient for bfheating.
auto integrand_bfheatingcoeff(const double nu, const double nu_edge, const int nonemptymgi, const float T_R,
                              const std::span<const float> photoion_xs) -> double {
  const float sigma_bf = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu);

  return sigma_bf * (1 - (nu_edge / nu)) * radfield::radfield(nu, nonemptymgi) * (1 - exp(-HOVERKB * nu / T_R));
}

// Sum the thermal heating rate from collisional de-excitations over all levels of a single ion (the energy given
// to the thermal electron pool when bound electrons are collisionally de-excited).
auto get_heating_ion_coll_deexc(const int nonemptymgi, const int element, const int ion, const float T_e,
                                const float clumpednne) -> double {
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
          nnlevel *
          col_deexcitation_ratecoeff(T_e, clumpednne, epsilon_trans, statweight, lower_statweight, alltransindex) *
          epsilon_trans;
      C_deexc += C;
    }
  }
  return C_deexc;
}

// Calculate the heating rates for a given cell. Results are returned via the elements of the heatingrates data
// structure.
void calculate_heating_rates(const int nonemptymgi, const float T_e, const float clumpednne,
                             HeatingCoolingRates& heatingcoolingrates, const std::span<const double> bfheatingcoeffs) {
  double C_deexc = 0.;

  double bfheating = 0.;
  double ffheating = 0.;

  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    if constexpr (DIRECT_COL_HEAT) {
      for (int ion = 0; ion < nions; ion++) {
        C_deexc += get_heating_ion_coll_deexc(nonemptymgi, element, ion, T_e, clumpednne);
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
    // from Monte Carlo estimators, which accumulate collisional recombination heating as well as
    // the collisional de-excitation heating that the DIRECT_COL_HEAT branch above sums analytically
    heatingcoolingrates.heating_collisional = globals::colheatingestimator.at(nonemptymgi);
  }

  heatingcoolingrates.heating_bf = bfheating;
  heatingcoolingrates.heating_ff = ffheating;
}

// Residual (total heating minus total cooling) of the thermal balance equation, as a function of T_e.
// NB: this is not a pure function of T_e. Evaluating it stores T_e in the grid and re-solves the cell's
// ionisation balance, populations and nne at that temperature, so the cell is left in the state belonging
// to the last T_e passed in. call_T_e_finder() relies on this and re-evaluates at the final T_e.
auto T_e_eqn_heating_minus_cooling(const double T_e, int nonemptymgi, const double t_current,
                                   HeatingCoolingRates& heatingcoolingrates,
                                   const std::span<const double> bfheatingcoeffs) -> double {
  const auto fT_e = static_cast<float>(T_e);

  if constexpr (!LTEPOP_EXCITATION_USE_TJ) {
    if (std::abs((T_e / grid::Te_allcells[nonemptymgi]) - 1.) > 0.1) {
      grid::Te_allcells[nonemptymgi] = fT_e;
      for (int element = 0; element < get_nelements(); element++) {
        if (!elem_has_nlte_levels(element)) {
          // recalculate the Gammas using the current level populations
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions - 1; ion++) {
            const auto groundcontindex = get_groundcontindex(element, ion);
            if (groundcontindex >= 0) {
              globals::gammaestimator[(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) +
                                      groundcontindex] = calculate_iongamma_per_gspop(nonemptymgi, element, ion);
            }
          }
        }
      }
    }
  }

  // Set new T_e guess for the current cell and update populations
  grid::Te_allcells[nonemptymgi] = fT_e;

  calculate_ion_balance_nne(nonemptymgi);
  const auto nne = grid::get_nne(nonemptymgi);
  const float clumpfactor = grid::get_clumpfactor(nonemptymgi);

  // Then calculate heating and cooling rates
  kpkt::calculate_cooling_rates(nonemptymgi, &heatingcoolingrates);
  calculate_heating_rates(nonemptymgi, fT_e, clumpfactor * nne, heatingcoolingrates, bfheatingcoeffs);

  const auto ntlepton_frac_heating = nonthermal::get_nt_frac_heating(nonemptymgi);
  const auto ntlepton_dep = nonthermal::get_ntlepton_deposition_rate_density(nonemptymgi);
  const auto ntalpha_dep = heatingcoolingrates.dep_alpha;
  // spontaneous fission fragments deposit as pure heating via k-packets (like alpha particles)
  const auto ntspfission_dep = heatingcoolingrates.dep_spfission;
  heatingcoolingrates.heating_dep = (ntlepton_dep * ntlepton_frac_heating) + ntalpha_dep + ntspfission_dep;
  heatingcoolingrates.dep_frac_heating =
      ((ntalpha_dep + ntspfission_dep) > 0)
          ? heatingcoolingrates.heating_dep / (ntlepton_dep + ntalpha_dep + ntspfission_dep)
          : ntlepton_frac_heating;

  // Adiabatic cooling p dV/dt per unit volume. Homologous expansion gives V proportional to t^3, so this
  // reduces to 3p/t, but it is written out below in terms of the cell volume at tmin.
  const double nntot = get_nnion_tot(nonemptymgi) + nne;
  const double p = nntot * KB * T_e;  // ideal gas pressure in [erg/cm^3]
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double volumetmin = grid::get_modelcell_assocvolume_tmin(modelgridindex);
  const double dV_on_dt = 3 * volumetmin / pow3(globals::tmin) * pow2(t_current);
  const double V = volumetmin * pow3(t_current / globals::tmin);
  heatingcoolingrates.cooling_adiabatic = p * dV_on_dt / V;

  const double total_heating_rate = heatingcoolingrates.heating_ff + heatingcoolingrates.heating_bf +
                                    heatingcoolingrates.heating_collisional + heatingcoolingrates.heating_dep;
  const double total_coolingrate = heatingcoolingrates.cooling_ff + heatingcoolingrates.cooling_fb +
                                   heatingcoolingrates.cooling_collisional + heatingcoolingrates.cooling_adiabatic;

  return total_heating_rate - total_coolingrate;
}

}  // anonymous namespace

// Compute the bound-free (photoionisation) heating coefficient for one level's photoionisation target by
// integrating the heating integrand over the cross-section's frequency range for the cell's radiation field.
auto calculate_bfheatingcoeff(const int element, const int ion, const int level, const int phixstargetindex,
                              const int nonemptymgi) -> double {
  double error = 0.;
  const double epsrel = 1e-3;

  const double E_threshold = get_phixs_threshold(element, ion, level, phixstargetindex);

  const double nu_threshold = (1. / H) * E_threshold;
  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table
  const auto photoion_xs = get_phixs_table(element, ion, level);
  const auto T_R = grid::TR_allcells[nonemptymgi];
  auto bfheating = integrator(
      [&](const double nu) { return integrand_bfheatingcoeff(nu, nu_threshold, nonemptymgi, T_R, photoion_xs); },
      nu_threshold, nu_max_phixs, epsrel, &error);

  bfheating *= FOURPI * get_phixsprobability(element, ion, level, phixstargetindex);

  return bfheating;
}

// Calculate the bound-free heating coefficient of every level in a cell. These depend only on the radiation
// field, not on T_e or the populations, so they are computed once per cell and reused at every T_e that the
// temperature solver tries.
void calculate_bfheatingcoeffs(int nonemptymgi, std::span<double> bfheatingcoeffs) {
  assert_always(std::ssize(bfheatingcoeffs) == get_includedlevels());
  const double minelfrac = 0.01;
  for (int element = 0; element < get_nelements(); element++) {
    if (grid::get_elem_massfrac(nonemptymgi, element) <= minelfrac && !USE_ION_BFHEATING_ESTIMATORS) {
      printlog("skipping Z={} X={:g}, ", get_atomicnumber(element), grid::get_elem_massfrac(nonemptymgi, element));
    }

    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);
      const auto levels = std::ranges::iota_view{0, nlevels};
      std::for_each(EXEC_PAR levels.begin(), levels.end(), [&](const int level) {
        double bfheatingcoeff = 0.;
        if (grid::get_elem_massfrac(nonemptymgi, element) > minelfrac || USE_ION_BFHEATING_ESTIMATORS) {
          const auto nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            bfheatingcoeff += calculate_bfheatingcoeff(element, ion, level, phixstargetindex, nonemptymgi);
          }
          assert_always(std::isfinite(bfheatingcoeff));

          if constexpr (USE_ION_BFHEATING_ESTIMATORS) {
            const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
            const int index_in_groundlevelcontestimator = globals::alllevels.closestgroundlevelcont[uniquelevelindex];
            if (index_in_groundlevelcontestimator >= 0) {
              bfheatingcoeff *=
                  globals::bfheatingestimator[(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) +
                                              index_in_groundlevelcontestimator];
            }
          }
        }
        bfheatingcoeffs[get_uniquelevelindex(element, ion, level)] = bfheatingcoeff;
      });
    }
  }
}

// Solve the thermal-balance equation (heating = cooling) for the electron temperature T_e in a cell by
// root-finding between MINTEMP and MAXTEMP, then store the resulting T_e in the grid. If the equation has no
// sign change over that range, T_e is pinned to whichever bound the residual points towards. The change from
// the previous timestep's T_e is then damped to at most a factor of two in either direction.
void call_T_e_finder(const int nonemptymgi, const double t_current, HeatingCoolingRates& heatingcoolingrates,
                     const std::span<const double> bfheatingcoeffs) {
  const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double T_e_old = grid::Te_allcells[nonemptymgi];
  printlog("Finding T_e in cell {} at timestep {}...", modelgridindex, globals::timestep);

  const auto f_T_e = [&](double T_e) -> double {
    return T_e_eqn_heating_minus_cooling(T_e, nonemptymgi, t_current, heatingcoolingrates, bfheatingcoeffs);
  };

  const double f_T_min = f_T_e(MINTEMP);
  const double f_T_max = f_T_e(MAXTEMP);

  const bool invalid_values = (!std::isfinite(f_T_min) || !std::isfinite(f_T_max));
  if (invalid_values) {
    printlnlog(
        "[warning] call_T_e_finder: non-finite results in modelcell {} (T_R={:g}, W={:g}). T_e forced to be MINTEMP",
        modelgridindex, grid::TR_allcells[nonemptymgi], grid::W_allcells[nonemptymgi]);
  }

  double T_e{NAN};
  // a sign change over [MINTEMP, MAXTEMP] guarantees a root that the bracketing solver can find
  if (!invalid_values && f_T_min * f_T_max < 0) {
    const auto maxit = 100U;

    // TOMS 748 (Alefeld, Potra & Shi 1995, ACM Trans. Math. Softw. 21, 327, doi:10.1145/210089.210111):
    // bracketing solver with inverse cubic interpolation, so it keeps the root bracketed like bisection
    // but converges superlinearly on the smooth part of the residual
    uintmax_t iternum = maxit;
    auto result = toms748_solve(f_T_e, MINTEMP, MAXTEMP, f_T_min, f_T_max, ftol<TEMPERATURE_SOLVER_ACCURACY>, iternum);
    T_e = 0.5 * (result.first + result.second);
    if (iternum >= maxit) {
      printlnlog("[warning] call_T_e_finder: T_e did not converge within {} iterations. interval [{:g}, {:g}] [K]",
                 iternum, result.first, result.second);
    } else {
      printlnlog("after {} iterations, T_e = {:g} [K], interval [{:g}, {:g}] [K]", iternum, T_e, result.first,
                 result.second);
    }
  } else if (invalid_values || f_T_max < 0) {
    // Thermal balance equation always negative ===> T_e = T_min
    T_e = MINTEMP;
    printlnlog(
        "[warning] call_T_e_finder: cooling bigger than heating at lower T_e boundary {:g} in modelcell {} "
        "(T_R={:g},W={:g}). T_e forced to be MINTEMP",
        MINTEMP, modelgridindex, grid::TR_allcells[nonemptymgi], grid::W_allcells[nonemptymgi]);
  } else {
    // Thermal balance equation always positive ===> T_e = T_max
    T_e = MAXTEMP;
    printlnlog(
        "[warning] call_T_e_finder: heating bigger than cooling over the whole T_e range [{:g},{:g}] in modelcell {} "
        "(T_R={:g},W={:g}). T_e forced to be MAXTEMP",
        MINTEMP, MAXTEMP, modelgridindex, grid::TR_allcells[nonemptymgi], grid::W_allcells[nonemptymgi]);
  }

  if (T_e > 2 * T_e_old) {
    const double T_e_solved = T_e;
    T_e = std::min(2 * T_e_old, MAXTEMP);
    printlnlog(
        "T_e damping in cell {} at timestep {}: solver gave T_e {:g} [K], clamping the rise from the previous value "
        "{:g} [K] to {:g} [K]",
        modelgridindex, globals::timestep, T_e_solved, T_e_old, T_e);
  } else if (T_e < 0.5 * T_e_old) {
    const double T_e_solved = T_e;
    T_e = std::max(0.5 * T_e_old, MINTEMP);
    printlnlog(
        "T_e damping in cell {} at timestep {}: solver gave T_e {:g} [K], clamping the fall from the previous value "
        "{:g} [K] to {:g} [K]",
        modelgridindex, globals::timestep, T_e_solved, T_e_old, T_e);
  }

  grid::Te_allcells[nonemptymgi] = static_cast<float>(T_e);

  // this call will make sure heating/cooling rates and populations are updated for the final T_e
  // in case T_e got modified after the T_e solver finished
  f_T_e(T_e);
}
