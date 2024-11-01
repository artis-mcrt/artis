#include "thermalbalance.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <algorithm>
#include <cmath>
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
  int modelgridindex;
  HeatingCoolingRates *heatingcoolingrates;
  const std::vector<double> *bfheatingcoeffs;
};

struct BFHeatingIntegralParams {
  double nu_edge;
  int nonemptymgi;
  float T_R;
  const float *photoion_xs;
};

auto integrand_bfheatingcoeff_custom_radfield(const double nu, void *const voidparas) -> double
// Integrand to calculate the rate coefficient for bfheating using gsl integrators.
{
  const auto *const params = static_cast<const BFHeatingIntegralParams *>(voidparas);

  const int nonemptymgi = params->nonemptymgi;
  const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double nu_edge = params->nu_edge;
  const float T_R = params->T_R;
  // const double Te_TR_factor = params->Te_TR_factor; // = sqrt(T_e/T_R) * sahafac(Te) / sahafac(TR)

  const float sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge, nu);

  // const auto T_e = grid::get_Te(modelgridindex);
  // return sigma_bf * (1 - nu_edge/nu) * radfield::radfield(nu,modelgridindex) * (1 - Te_TR_factor * exp(-HOVERKB * nu
  // / T_e));

  return sigma_bf * (1 - nu_edge / nu) * radfield::radfield(nu, modelgridindex) * (1 - exp(-HOVERKB * nu / T_R));
}

auto calculate_bfheatingcoeff(const int element, const int ion, const int level, const int phixstargetindex,
                              const int modelgridindex) -> double {
  double error = 0.;
  const double epsrel = 1e-3;
  const double epsrelwarning = 1e-1;
  const double epsabs = 0.;

  // const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
  // const double E_threshold = epsilon(element, ion + 1, upperionlevel) - epsilon(element, ion, level);
  const double E_threshold = get_phixs_threshold(element, ion, level, phixstargetindex);

  const double nu_threshold = ONEOVERH * E_threshold;
  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

  // const auto T_e = grid::get_Te(modelgridindex);
  // const double T_R = grid::get_TR(nonemptymgi);
  // const double sf_Te = calculate_sahafact(element,ion,level,upperionlevel,T_e,E_threshold);
  // const double sf_TR = calculate_sahafact(element,ion,level,upperionlevel,T_R,E_threshold);
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  const BFHeatingIntegralParams intparas = {.nu_edge = nu_threshold,
                                            .nonemptymgi = nonemptymgi,
                                            .T_R = grid::get_TR(nonemptymgi),
                                            .photoion_xs = get_phixs_table(element, ion, level)};

  // intparas.Te_TR_factor = sqrt(T_e/T_R) * sf_Te / sf_TR;

  double bfheating = 0.;

  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);

  const int status = integrator<integrand_bfheatingcoeff_custom_radfield>(
      intparas, nu_threshold, nu_max_phixs, epsabs, epsrel, GSL_INTEG_GAUSS61, &bfheating, &error);

  if (status != 0 && (status != 18 || (error / bfheating) > epsrelwarning)) {
    printout(
        "bf_heating integrator gsl warning %d. modelgridindex %d Z=%d ionstage %d lower %d phixstargetindex %d "
        "integral %g error %g\n",
        status, modelgridindex, get_atomicnumber(element), get_ionstage(element, ion), level, phixstargetindex,
        bfheating, error);
  }
  gsl_set_error_handler(previous_handler);

  bfheating *= FOURPI * get_phixsprobability(element, ion, level, phixstargetindex);

  return bfheating;
}

auto get_heating_ion_coll_deexc(const int modelgridindex, const int element, const int ion, const double T_e,
                                const double nne) -> double {
  double C_deexc = 0.;
  const int nlevels = get_nlevels(element, ion);

  for (int level = 0; level < nlevels; level++) {
    const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
    const double epsilon_level = epsilon(element, ion, level);

    // Collisional heating: deexcitation to same ionization stage
    const int ndowntrans = get_ndowntrans(element, ion, level);
    const auto *const leveldowntranslist = get_downtranslist(element, ion, level);
    for (int i = 0; i < ndowntrans; i++) {
      const auto &downtransition = leveldowntranslist[i];
      const int lower = downtransition.targetlevelindex;
      const double epsilon_trans = epsilon_level - epsilon(element, ion, lower);
      const double C = nnlevel *
                       col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, element, ion, level, downtransition) *
                       epsilon_trans;
      C_deexc += C;
    }
  }
  return C_deexc;
}

// Calculate the heating rates for a given cell. Results are returned via the elements of the heatingrates data
// structure.
void calculate_heating_rates(const int modelgridindex, const double T_e, const double nne,
                             HeatingCoolingRates *heatingcoolingrates, const std::vector<double> &bfheatingcoeffs) {
  double C_deexc = 0.;

  // double C_recomb = 0.;
  double bfheating = 0.;
  double ffheating = 0.;

  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    if constexpr (DIRECT_COL_HEAT) {
      for (int ion = 0; ion < nions; ion++) {
        C_deexc += get_heating_ion_coll_deexc(modelgridindex, element, ion, T_e, nne);
      }
    }

    // Collisional heating: recombination to lower ionization stage (not included)
    // ------------------------------------------------------------

    // Bound-free heating (renormalised analytical calculation)
    // --------------------------------------------------------
    // We allow bound free-transitions only if there is a higher ionisation stage
    // left in the model atom to match the bound-free absorption in the rpkt routine.
    // There this condition is needed as we can only ionise to existing ionisation
    // stage even if there would be further ionisation stages in nature which
    // are not included in the model atom.

    for (int ion = 0; ion < nions - 1; ion++) {
      const int nbflevels = get_nlevels_ionising(element, ion);
      for (int level = 0; level < nbflevels; level++) {
        const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
        bfheating += nnlevel * bfheatingcoeffs[get_uniquelevelindex(element, ion, level)];
      }
    }
  }

  // Free-free heating (from estimators)

  const int nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  ffheating = globals::ffheatingestimator[nonemptymgi];

  if constexpr (DIRECT_COL_HEAT) {
    heatingcoolingrates->heating_collisional = C_deexc;
  } else {
    // Collisional heating (from estimators)
    heatingcoolingrates->heating_collisional = globals::colheatingestimator[nonemptymgi];  // C_deexc + C_recomb;
  }

  heatingcoolingrates->heating_bf = bfheating;
  heatingcoolingrates->heating_ff = ffheating;
}

// Thermal balance equation on which we have to iterate to get T_e
auto T_e_eqn_heating_minus_cooling(const double T_e, void *paras) -> double {
  const TeSolutionParams *const params = static_cast<TeSolutionParams *>(paras);

  const int modelgridindex = params->modelgridindex;
  const double t_current = params->t_current;
  auto *const heatingcoolingrates = params->heatingcoolingrates;

  // Set new T_e guess for the current cell and update populations
  grid::set_Te(modelgridindex, T_e);

  if constexpr (!USE_LUT_PHOTOION && !LTEPOP_EXCITATION_USE_TJ) {
    const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
    for (int element = 0; element < get_nelements(); element++) {
      if (!elem_has_nlte_levels(element)) {
        // recalculate the Gammas using the current population estimates
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions - 1; ion++) {
          globals::gammaestimator[get_ionestimindex_nonemptymgi(nonemptymgi, element, ion)] =
              calculate_iongamma_per_gspop(modelgridindex, element, ion);
        }
      }
    }
  }

  calculate_ion_balance_nne(modelgridindex);
  const auto nne = grid::get_nne(modelgridindex);

  // Then calculate heating and cooling rates
  const int nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  kpkt::calculate_cooling_rates(nonemptymgi, heatingcoolingrates);
  calculate_heating_rates(modelgridindex, T_e, nne, heatingcoolingrates, *params->bfheatingcoeffs);

  const auto ntlepton_frac_heating = nonthermal::get_nt_frac_heating(modelgridindex);
  const auto ntlepton_dep = nonthermal::get_deposition_rate_density(modelgridindex);
  const auto ntalpha_frac_heating = 1.;
  const auto ntalpha_dep = heatingcoolingrates->dep_alpha;
  heatingcoolingrates->heating_dep = ntlepton_dep * ntlepton_frac_heating + ntalpha_dep * ntalpha_frac_heating;
  heatingcoolingrates->dep_frac_heating =
      (ntalpha_dep > 0) ? heatingcoolingrates->heating_dep / (ntlepton_dep + ntalpha_dep) : ntlepton_frac_heating;

  // Adiabatic cooling term
  const double nntot = get_nnion_tot(modelgridindex) + nne;
  const double p = nntot * KB * T_e;  // pressure in [erg/cm^3]
  const double volumetmin = grid::get_modelcell_assocvolume_tmin(modelgridindex);
  const double dV = 3 * volumetmin / pow(globals::tmin, 3) * pow(t_current, 2);  // really dV/dt
  const double V = volumetmin * pow(t_current / globals::tmin, 3);
  heatingcoolingrates->cooling_adiabatic = p * dV / V;

  const double total_heating_rate = heatingcoolingrates->heating_ff + heatingcoolingrates->heating_bf +
                                    heatingcoolingrates->heating_collisional + heatingcoolingrates->heating_dep;
  const double total_coolingrate = heatingcoolingrates->cooling_ff + heatingcoolingrates->cooling_fb +
                                   heatingcoolingrates->cooling_collisional + heatingcoolingrates->cooling_adiabatic;

  return total_heating_rate - total_coolingrate;
}

}  // anonymous namespace

auto get_bfheatingcoeff_ana(const int element, const int ion, const int level, const int phixstargetindex,
                            const double T_R, const double W) -> double {
  // The correction factor for stimulated emission in gammacorr is set to its
  // LTE value. Because the T_e dependence of gammacorr is weak, this correction
  // correction may be evaluated at T_R!
  assert_always(USE_LUT_BFHEATING);
  double bfheatingcoeff = 0.;

  const int lowerindex = floor(log(T_R / MINTEMP) / T_step_log);
  if (lowerindex < TABLESIZE - 1) {
    const int upperindex = lowerindex + 1;
    const double T_lower = MINTEMP * exp(lowerindex * T_step_log);
    const double T_upper = MINTEMP * exp(upperindex * T_step_log);

    const double f_upper = globals::bfheating_coeff[get_bflutindex(upperindex, element, ion, level, phixstargetindex)];
    const double f_lower = globals::bfheating_coeff[get_bflutindex(lowerindex, element, ion, level, phixstargetindex)];

    bfheatingcoeff = (f_lower + (f_upper - f_lower) / (T_upper - T_lower) * (T_R - T_lower));
  } else {
    bfheatingcoeff = globals::bfheating_coeff[get_bflutindex(TABLESIZE - 1, element, ion, level, phixstargetindex)];
  }

  return W * bfheatingcoeff;
}

// depends only the radiation field - no dependence on T_e or populations
void calculate_bfheatingcoeffs(int nonemptymgi, std::vector<double> &bfheatingcoeffs) {
  bfheatingcoeffs.resize(get_includedlevels());
  const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double minelfrac = 0.01;
  for (int element = 0; element < get_nelements(); element++) {
    if (grid::get_elem_abundance(nonemptymgi, element) <= minelfrac && !USE_LUT_BFHEATING) {
      printout("skipping Z=%d X=%g, ", get_atomicnumber(element), grid::get_elem_abundance(nonemptymgi, element));
    }

    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        double bfheatingcoeff = 0.;
        if (grid::get_elem_abundance(nonemptymgi, element) > minelfrac || USE_LUT_BFHEATING) {
          const auto nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            if constexpr (!USE_LUT_BFHEATING) {
              bfheatingcoeff += calculate_bfheatingcoeff(element, ion, level, phixstargetindex, modelgridindex);
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
            const int index_in_groundlevelcontestimator =
                globals::elements[element].ions[ion].levels[level].closestgroundlevelcont;
            if (index_in_groundlevelcontestimator >= 0) {
              bfheatingcoeff *= globals::bfheatingestimator[(nonemptymgi * globals::nbfcontinua_ground) +
                                                            index_in_groundlevelcontestimator];
            }
          }
        }
        bfheatingcoeffs[get_uniquelevelindex(element, ion, level)] = bfheatingcoeff;
      }
    }
  }
}

void call_T_e_finder(const int nonemptymgi, const double t_current, const double T_min, const double T_max,
                     HeatingCoolingRates *heatingcoolingrates, const std::vector<double> &bfheatingcoeffs) {
  const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double T_e_old = grid::get_Te(modelgridindex);
  printout("Finding T_e in cell %d at timestep %d...", modelgridindex, globals::timestep);

  TeSolutionParams paras = {.t_current = t_current,
                            .modelgridindex = modelgridindex,
                            .heatingcoolingrates = heatingcoolingrates,
                            .bfheatingcoeffs = &bfheatingcoeffs};

  gsl_function find_T_e_f = {.function = &T_e_eqn_heating_minus_cooling, .params = &paras};

  double thermalmin = T_e_eqn_heating_minus_cooling(T_min, find_T_e_f.params);
  double thermalmax = T_e_eqn_heating_minus_cooling(T_max, find_T_e_f.params);

  if (!std::isfinite(thermalmin) || !std::isfinite(thermalmax)) {
    printout(
        "[abort request] call_T_e_finder: non-finite results in modelcell %d (T_R=%g, W=%g). T_e forced to be "
        "MINTEMP\n",
        modelgridindex, grid::get_TR(nonemptymgi), grid::get_W(nonemptymgi));
    thermalmax = thermalmin = -1;
  }

  double T_e{NAN};
  // Check whether the thermal balance equation has a root in [T_min,T_max]
  if (thermalmin * thermalmax < 0) {
    // If it has, then solve for the root T_e

    // one-dimensional gsl root solver, bracketing type
    gsl_root_fsolver *T_e_solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

    gsl_root_fsolver_set(T_e_solver, &find_T_e_f, T_min, T_max);
    const int maxit = 100;
    int status = 0;
    for (int iternum = 0; iternum < maxit; iternum++) {
      gsl_root_fsolver_iterate(T_e_solver);
      T_e = gsl_root_fsolver_root(T_e_solver);
      const double T_e_min = gsl_root_fsolver_x_lower(T_e_solver);
      const double T_e_max = gsl_root_fsolver_x_upper(T_e_solver);
      status = gsl_root_test_interval(T_e_min, T_e_max, 0, TEMPERATURE_SOLVER_ACCURACY);
      // printout("iter %d, T_e interval [%g, %g], guess %g, status %d\n", iternum, T_e_min, T_e_max, T_e, status);
      if (status != GSL_CONTINUE) {
        printout("after %d iterations, T_e = %g K, interval [%g, %g]\n", iternum + 1, T_e, T_e_min, T_e_max);
        break;
      }
    }

    if (status == GSL_CONTINUE) {
      printout("[warning] call_T_e_finder: T_e did not converge within %d iterations\n", maxit);
    }

    gsl_root_fsolver_free(T_e_solver);
  }
  // Quick solver style: works if we can assume that there is either one or no
  // solution on [MINTEMP.MAXTEMP] (check that by doing a plot of heating-cooling vs. T_e)
  else if (thermalmax < 0) {
    // Thermal balance equation always negative ===> T_e = T_min
    T_e = MINTEMP;
    printout(
        "[warning] call_T_e_finder: cooling bigger than heating at lower T_e boundary %g in modelcell %d "
        "(T_R=%g,W=%g). T_e forced to be MINTEMP\n",
        MINTEMP, modelgridindex, grid::get_TR(nonemptymgi), grid::get_W(nonemptymgi));
  } else {
    // Thermal balance equation always negative ===> T_e = T_max
    T_e = MAXTEMP;
    printout(
        "[warning] call_T_e_finder: heating bigger than cooling over the whole T_e range [%g,%g] in modelcell %d "
        "(T_R=%g,W=%g). T_e forced to be MAXTEMP\n",
        MINTEMP, MAXTEMP, modelgridindex, grid::get_TR(nonemptymgi), grid::get_W(nonemptymgi));
  }

  if (T_e > 2 * T_e_old) {
    T_e = 2 * T_e_old;
    printout("use T_e damping in cell %d\n", modelgridindex);
    T_e = std::min(T_e, MAXTEMP);
  } else if (T_e < 0.5 * T_e_old) {
    T_e = 0.5 * T_e_old;
    printout("use T_e damping in cell %d\n", modelgridindex);
    T_e = std::max(T_e, MINTEMP);
  }

  grid::set_Te(modelgridindex, T_e);

  // this call with make sure heating/cooling rates and populations are updated for the final T_e
  // in case T_e got modified after the T_e solver finished
  T_e_eqn_heating_minus_cooling(T_e, find_T_e_f.params);
}
