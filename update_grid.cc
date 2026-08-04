// The per-timestep grid update: for each nonempty model cell, normalises the Monte Carlo
// estimators collected during the previous timestep and solves for the new radiation field
// parameters, electron temperature, ion/level populations, and cooling rates.

#include "update_grid.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <ostream>
#include <print>
#include <span>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "kpkt.h"
#include "ltepop.h"
#include "mpi_logging.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "thermalbalance.h"
#include "vpkt.h"

namespace {

std::vector<HeatingCoolingRates> heatingcoolingrates_thisrankcells;

void write_to_estimators_file(std::ostream& estimators_file, const int nonemptymgi, const int timestep, const int titer,
                              const HeatingCoolingRates& heatingcoolingrates,
                              const decay::NucMassFracCoeffsSpan nuc_massfrac_coeffs) {
  // writing the estimators file is a measurable cost per cell per timestep. If the estimators output is not
  // needed, an early return here skips it.
  const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);

  const auto sys_time_start_write_estimators = std::chrono::steady_clock::now();

  const auto T_e = grid::Te_allcells[nonemptymgi];
  const auto nne = grid::get_nne(nonemptymgi);
  const auto Y_e = grid::get_electronfrac(nonemptymgi);

  std::print(estimators_file, "timestep {} modelgridindex {} titeration {} TR {:g} Te {:g} W {:g} TJ {:g}", timestep,
             mgi, titer, grid::TR_allcells[nonemptymgi], T_e, grid::W_allcells[nonemptymgi],
             grid::TJ_allcells[nonemptymgi]);
  std::println(estimators_file, " grey_depth {:g} thick {} nne {:g} Ye {:g} tdays {:7.2f}",
               grid::grey_depth_allcells[nonemptymgi], grid::thick_allcells[nonemptymgi], nne, Y_e,
               globals::timesteps[timestep].mid / DAY);

  if (globals::total_nlte_levels > 0) {
    nltepop_write_to_file(nonemptymgi, timestep);
  }

  for (int element = 0; element < get_nelements(); element++) {
    if (grid::get_elem_massfrac(nonemptymgi, element) <= 0.) {  // skip elements with no abundance
      continue;
    }

    std::print(estimators_file, "populations        Z={:2d}", get_atomicnumber(element));
    const int nions = get_nions(element);
    if (nions > 0) {
      // add spaces for missing lowest ion stages to match other elements
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
        std::print(estimators_file, "              ");
      }
    }
    double elpop = 0.;
    for (int ion = 0; ion < nions; ion++) {
      elpop += get_nnion(nonemptymgi, element, ion);
      std::print(estimators_file, "  {}: {:9.3e}", get_ionstage(element, ion), get_nnion(nonemptymgi, element, ion));
    }
    if (nions == 0) {
      elpop = grid::get_elem_numberdens(nonemptymgi, element);
    }
    std::print(estimators_file, "  SUM: {:9.3e}", elpop);

    decay::output_isotopic_densities(estimators_file, nonemptymgi, element, nuc_massfrac_coeffs);

    if (nions == 0 || elpop <= 0.) {
      // dummy element for nuclear abundances only
      continue;
    }

    std::print(estimators_file, "gamma_R            Z={:2d}", get_atomicnumber(element));
    for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
      std::print(estimators_file, "              ");
    }
    for (int ion = 0; ion < nions - 1; ion++) {
      std::print(estimators_file, "  {}: {:9.3e}", get_ionstage(element, ion),
                 calculate_iongamma_per_ionpop(nonemptymgi, element, ion, false, false));
    }
    std::println(estimators_file);

    // if we have bound-free estimators, we might want to compare how gamma_R would be if we used the radiation field
    // model instead
    if (DETAILED_BF_ESTIMATORS_ON) {
      std::print(estimators_file, "gamma_R_integral   Z={:2d}", get_atomicnumber(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
        std::print(estimators_file, "              ");
      }
      for (int ion = 0; ion < nions - 1; ion++) {
        std::print(estimators_file, "  {}: {:9.3e}", get_ionstage(element, ion),
                   calculate_iongamma_per_ionpop(nonemptymgi, element, ion, false, true));
      }
      std::println(estimators_file);
    }

    if (NT_ON) {
      std::print(estimators_file, "gamma_NT           Z={:2d}", get_atomicnumber(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
        std::print(estimators_file, "              ");
      }
      for (int ion = 0; ion < nions - 1; ion++) {
        const double Y_nt = nonthermal::nt_ionisation_ratecoeff(nonemptymgi, element, ion);
        std::print(estimators_file, "  {}: {:9.3e}", get_ionstage(element, ion), Y_nt);
      }
      std::println(estimators_file);
    }

    if (USE_LUT_PHOTOION && globals::nbfcontinua_ground > 0) {
      std::print(estimators_file, "corrphotoionrenorm Z={:2d}", get_atomicnumber(element));
      for (int ion = 0; ion < nions - 1; ion++) {
        if (get_groundcontindex(element, ion) >= 0) {
          std::print(estimators_file, "  {}: {:9.3e}", get_ionstage(element, ion),
                     globals::corrphotoionrenorm[(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) +
                                                 get_groundcontindex(element, ion)]);
        }
      }
      std::println(estimators_file);
      std::print(estimators_file, "gammaestimator     Z={:2d}", get_atomicnumber(element));
      for (int ion = 0; ion < nions - 1; ion++) {
        if (get_groundcontindex(element, ion) >= 0) {
          std::print(estimators_file, "  {}: {:9.3e}", get_ionstage(element, ion),
                     globals::gammaestimator[(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) +
                                             get_groundcontindex(element, ion)]);
        }
      }
      std::println(estimators_file);
    }
  }

  // power densities in erg / s / cm^3. 'ana' means analytical at t_mid, i.e. the rates calculated from the nuclear
  // abundances and decay data, not from Monte Carlo events
  std::print(estimators_file, "emission_ana: gamma {:11.5e} positron {:11.5e} electron {:11.5e} alpha {:11.5e}",
             heatingcoolingrates.eps_gamma_ana, heatingcoolingrates.eps_positron_ana,
             heatingcoolingrates.eps_electron_ana, heatingcoolingrates.eps_alpha_ana);
  const bool any_fission = decay::decaytype_is_used(decay::DECAYTYPE_SPONTFISSION);
  if (any_fission) {
    std::print(estimators_file, " spfission {:11.5e}", heatingcoolingrates.eps_spfission_ana);
  }
  std::println(estimators_file);

  std::print(estimators_file, "deposition: gamma {:11.5e} positron {:11.5e} electron {:11.5e} alpha {:11.5e}",
             heatingcoolingrates.dep_gamma, heatingcoolingrates.dep_positron, heatingcoolingrates.dep_electron,
             heatingcoolingrates.dep_alpha);
  if (any_fission) {
    std::print(estimators_file, " spfission {:11.5e}", heatingcoolingrates.dep_spfission);
  }
  std::println(estimators_file);

  std::print(estimators_file,
             "heating: ff {:11.5e} bf {:11.5e} coll {:11.5e}       dep {:11.5e} heating_dep/total_dep {:.3f}\n",
             heatingcoolingrates.heating_ff, heatingcoolingrates.heating_bf, heatingcoolingrates.heating_collisional,
             heatingcoolingrates.heating_dep, heatingcoolingrates.dep_frac_heating);
  std::print(estimators_file, "cooling: ff {:11.5e} fb {:11.5e} coll {:11.5e} adiabatic {:11.5e}\n",
             heatingcoolingrates.cooling_ff, heatingcoolingrates.cooling_fb, heatingcoolingrates.cooling_collisional,
             heatingcoolingrates.cooling_adiabatic);

  const auto write_estim_duration =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_write_estimators).count();
  if (write_estim_duration >= 1.) {
    printlnlog("writing estimators for timestep {} cell {} took {:.1f} seconds", timestep, mgi, write_estim_duration);
  }
}

void solve_Te_nltepops(const int nonemptymgi, const int nts, const int nts_prev,
                       HeatingCoolingRates& heatingcoolingrates) {
  const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  // bfheating coefficients are needed for the T_e solver, but they only depend on the radiation field, which is fixed
  // during the iterations below
  printlog("calculate_bfheatingcoeffs for timestep {} cell {}...", nts, mgi);
  const auto sys_time_start_calculate_bfheatingcoeffs = std::chrono::steady_clock::now();
  THREADLOCALONHOST auto bfheatingcoeffs = std::vector<double>(get_includedlevels());

  calculate_bfheatingcoeffs(nonemptymgi, bfheatingcoeffs);
  const auto bfheating_duration =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_calculate_bfheatingcoeffs)
          .count();
  printlnlog("took {:.1f} seconds", bfheating_duration);

  constexpr double nne_reltol = 0.04;
  constexpr double T_e_reltol = 0.04;
  for (int nlte_iter = 0; nlte_iter <= NLTEITER; nlte_iter++) {
    const auto sys_time_start_spencerfano = std::chrono::steady_clock::now();
    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      // SF solution depends on the ionisation balance, and weakly on nne
      nonthermal::solve_spencerfano(nonemptymgi, nts, nlte_iter);
    }
    const auto duration_solve_spencerfano =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_spencerfano).count();

    const auto sys_time_start_partfuncs_or_gamma = std::chrono::steady_clock::now();
    for (int element = 0; element < get_nelements(); element++) {
      if (!elem_has_nlte_levels(element)) {
        calculate_cellpartfuncts(nonemptymgi, element);
      }
    }
    const auto duration_solve_partfuncs_or_gamma =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_partfuncs_or_gamma).count();

    const double prev_T_e = grid::Te_allcells[nonemptymgi];
    const auto sys_time_start_Te = std::chrono::steady_clock::now();

    // Find T_e as solution for thermal balance
    call_T_e_finder(nonemptymgi, globals::timesteps[nts_prev].mid, heatingcoolingrates, bfheatingcoeffs);

    const auto duration_solve_T_e =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_Te).count();

    if (globals::total_nlte_levels == 0) {
      const auto sys_time_start_pops = std::chrono::steady_clock::now();
      calculate_ion_balance_nne(nonemptymgi);
      const auto duration_solve_pops =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_pops).count();

      printlnlog(
          "Grid solver cell {} timestep {}: time spent on: Spencer-Fano {:.1f}s, partfuncs/gamma {:.1f}s, T_e {:.1f}s, "
          "populations {:.1f}s",
          mgi, nts, duration_solve_spencerfano, duration_solve_partfuncs_or_gamma, duration_solve_T_e,
          duration_solve_pops);
      break;  // no iteration is needed without nlte pops
    }

    const double fracdiff_T_e = fabs((grid::Te_allcells[nonemptymgi] / prev_T_e) - 1);
    const auto sys_time_start_nltepops = std::chrono::steady_clock::now();
    // fractional difference between previous and current iteration's (nne or max(ground state population change))
    for (int element = 0; element < get_nelements(); element++) {
      if (get_nions(element) > 0 && elem_has_nlte_levels(element)) {
        solve_nlte_pops_element(element, nonemptymgi, nts, nlte_iter);
        calculate_cellpartfuncts(nonemptymgi, element);
      }
    }
    const auto duration_solve_nltepops =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_nltepops).count();

    const double nne_prev = grid::get_nne(nonemptymgi);
    calculate_ion_balance_nne(nonemptymgi);  // sets nne
    const auto fracdiff_nne = fabs((grid::get_nne(nonemptymgi) / nne_prev) - 1);
    printlnlog(
        "NLTE solver mgi {} timestep {} iteration {}: time spent on: Spencer-Fano {:.1f}s, T_e {:.1f}s, NLTE "
        "populations {:.1f}s",
        mgi, nts, nlte_iter, duration_solve_spencerfano, duration_solve_T_e, duration_solve_nltepops);
    printlnlog(
        "NLTE (Spencer-Fano/Te/pops) solver mgi {} timestep {} iteration {}: prev_iter nne {:e}, new nne is {:e}, "
        "fracdiff {:g}, prev T_e {:g} new T_e {:g} fracdiff {:g}",
        mgi, nts, nlte_iter, nne_prev, grid::get_nne(nonemptymgi), fracdiff_nne, prev_T_e,
        grid::Te_allcells[nonemptymgi], fracdiff_T_e);
    if (fracdiff_nne <= nne_reltol && fracdiff_T_e <= T_e_reltol) {
      printlnlog(
          "NLTE (Spencer-Fano/Te/pops) solver mgi {} timestep {} iteration {}: nne converged to tolerance {:g} <= {:g} "
          "and T_e to tolerance {:g} <= {:g}",
          mgi, nts, nlte_iter, fracdiff_nne, nne_reltol, fracdiff_T_e, T_e_reltol);
      break;
    }
    if (nlte_iter == NLTEITER) {
      printlnlog(
          "NLTE (Spencer-Fano/Te/pops) solver mgi {} timestep {} iteration {}: failed to converge... Keeping solution "
          "from last iteration",
          mgi, nts, nlte_iter);
    }
  }
}

void update_gamma_corrphotoionrenorm_bfheating_estimators(const int nonemptymgi, const double estimator_normfactor) {
  if constexpr (USE_LUT_PHOTOION) {
    const auto W = grid::W_allcells[nonemptymgi];
    const auto T_R = grid::TR_allcells[nonemptymgi];
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
        const int groundcontindex = get_groundcontindex(element, ion);
        if (groundcontindex < 0) {
          continue;
        }
        const ptrdiff_t ionestimindex =
            (static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) + groundcontindex;

        globals::gammaestimator[ionestimindex] *= estimator_normfactor / H;
#ifdef DO_TITER
        titer_average(globals::gammaestimator[ionestimindex], globals::gammaestimator_save[ionestimindex]);
#endif

        // renormalisation factor of the MC photoionisation rate estimate over the analytic rate for the cell's
        // dilute-blackbody radiation field. In cold and/or dilute cells the analytic rate can underflow to zero
        // for high-threshold continua (and the ratio can overflow for subnormal denominators). Fall back to
        // no renormalisation (factor 1) there, matching the initialisation value and the LTE/grey-cell treatment,
        // so the uncorrected analytic LUT rate is used for such continua. This drops any positive MC estimator
        // for the continuum, but no finite factor could preserve it: get_corrphotoioncoeff multiplies this factor
        // by the same analytic rate, so the product stays zero wherever that rate still underflows.
        const double gammacorr_ana = W * get_corrphotoioncoeff_ana(element, ion, 0, 0, T_R);
        // a non-finite analytic rate or MC estimator would be an upstream bug, not an underflowing denominator
        assert_always(std::isfinite(gammacorr_ana));
        assert_always(std::isfinite(globals::gammaestimator[ionestimindex]));
        const double renorm = (gammacorr_ana > 0.) ? globals::gammaestimator[ionestimindex] / gammacorr_ana
                                                   : std::numeric_limits<double>::infinity();
        if (std::isfinite(renorm)) {
          globals::corrphotoionrenorm[ionestimindex] = renorm;
        } else {
          printlnlog(
              "[warning] update_grid cell {} timestep {} Z={} ionstage {}: gammaestimator {:g} cannot be "
              "renormalised by W * corrphotoioncoeff_ana = {:g} (W {:g}, T_R {:g}). Setting corrphotoionrenorm = 1 "
              "(no renormalisation)",
              grid::get_mgi_of_nonemptymgi(nonemptymgi), globals::timestep, get_atomicnumber(element),
              get_ionstage(element, ion), globals::gammaestimator[ionestimindex], gammacorr_ana, W, T_R);
          globals::corrphotoionrenorm[ionestimindex] = 1.;
        }
      }

      // 2012-01-11. These loops should terminate here to precalculate *ALL* corrphotoionrenorm
      // values so that the values are known when required by the call to get_corrphotoioncoeff in
      // the following loops. Otherwise get_corrphotoioncoeff tries to renormalize by the closest
      // corrphotoionrenorm in frequency space which can lead to zero contributions to the total photoionisation rate!
    }
  }
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions - 1; ion++) {
      // Reuse the gammaestimator array as temporary storage of the Gamma values during
      // the remaining part of the update_grid phase. Afterwards it is reset to record
      // the next timesteps gamma estimators.
      const int groundcontindex = get_groundcontindex(element, ion);
      if (groundcontindex < 0) {
        continue;
      }
      const ptrdiff_t ionestimindex =
          (static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) + groundcontindex;

      if (!elem_has_nlte_levels(element)) {
        globals::gammaestimator[ionestimindex] = calculate_iongamma_per_gspop(nonemptymgi, element, ion);
      }

      if constexpr (USE_ION_BFHEATING_ESTIMATORS) {
        globals::bfheatingestimator[ionestimindex] *= estimator_normfactor;
#ifdef DO_TITER
        titer_average(globals::bfheatingestimator[ionestimindex], globals::bfheatingestimator_save[ionestimindex]);
#endif
        // Now convert bfheatingestimator into the bfheating renormalisation coefficient used in
        // for the remaining part of update_grid. At the start of the next update_packets, it will be reset

        // as for corrphotoionrenorm above, the analytic bf heating coefficient can be exactly zero in cells where
        // the radiation field model gives no flux above a continuum's threshold (e.g. empty radfield bins), so
        // guard the renormalisation and fall back to no renormalisation (factor 1) if the ratio is not finite.
        // A positive MC estimator over a zero analytic coefficient is dropped: calculate_bfheatingcoeffs multiplies
        // this factor by per-level coefficients from the same radiation field model, which are then zero too, so
        // no finite factor could carry the estimator's contribution through.
        const double bfheatingcoeff_ground = calculate_bfheatingcoeff(element, ion, 0, 0, nonemptymgi);
        // a non-finite analytic coefficient or MC estimator would be an upstream bug, not an underflow
        assert_always(std::isfinite(bfheatingcoeff_ground));
        assert_always(std::isfinite(globals::bfheatingestimator[ionestimindex]));
        const double bfheatingrenorm = (bfheatingcoeff_ground > 0.)
                                           ? globals::bfheatingestimator[ionestimindex] / bfheatingcoeff_ground
                                           : std::numeric_limits<double>::infinity();
        if (std::isfinite(bfheatingrenorm)) {
          globals::bfheatingestimator[ionestimindex] = bfheatingrenorm;
        } else {
          printlnlog(
              "[warning] update_grid cell {} timestep {} Z={} ionstage {}: bfheatingestimator {:g} cannot be "
              "renormalised by bfheatingcoeff_ground = {:g}. Setting bfheating renormalisation factor = 1",
              grid::get_mgi_of_nonemptymgi(nonemptymgi), globals::timestep, get_atomicnumber(element),
              get_ionstage(element, ion), globals::bfheatingestimator[ionestimindex], bfheatingcoeff_ground);
          globals::bfheatingestimator[ionestimindex] = 1.;
        }
      }
    }
  }
}

#ifdef DO_TITER
static void titer_average_estimators(const int nonemptymgi) {
  titer_average(globals::ffheatingestimator[nonemptymgi], globals::ffheatingestimator_save[nonemptymgi]);
  if constexpr (!DIRECT_COL_HEAT) {
    titer_average(globals::colheatingestimator[nonemptymgi], globals::colheatingestimator_save[nonemptymgi]);
  }
}
#endif

void update_grid_cell(const int nonemptymgi, const int nts, const int nts_prev, const int titer, const double tratmid,
                      const double deltat, HeatingCoolingRates& heatingcoolingrates,
                      const decay::AnaEmissionRateCoeffs& emissionratecoeffs) {
  const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);

  const double deltaV =
      grid::get_modelcell_assocvolume_tmin(mgi) * pow3(globals::timesteps[nts_prev].mid / globals::tmin);
  const auto sys_time_start_update_cell = std::chrono::steady_clock::now();

  const auto tmid = globals::timesteps[nts].mid;

  // Update clumping factors
  if constexpr (USE_MICROCLUMPING) {
    const double rad_vel = grid::get_modelcell_mean_radial_pos_tmin(mgi) / globals::tmin;
    const float clumpfactor = clumping_factor(tmid, rad_vel);

    grid::set_clumpfactor(nonemptymgi, clumpfactor);
  }

  nonthermal::calculate_deposition_rate_density(nonemptymgi, heatingcoolingrates, emissionratecoeffs);

  const double estimator_normfactor = 1 / deltaV / deltat / globals::nprocs;
  const double estimator_normfactor_over4pi = (1. / (4 * PI)) * estimator_normfactor;

  if (nts == globals::timestep_initial && titer == 0) {
    // For the initial timestep, temperatures have already been assigned
    // either by trapped energy release calculation, or reading from gridsave file

    if (USE_LUT_PHOTOION && !globals::simulation_continued_from_saved) {
      // Determine renormalisation factor for corrected photoionisation cross-sections
      std::ranges::fill(
          globals::corrphotoionrenorm.subspan(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground,
                                              globals::nbfcontinua_ground),
          1.);
    }

    // W == 1 indicates that this modelgrid cell was treated grey in the
    // last timestep. Therefore it has no valid Gamma estimators and must be treated in LTE at restart.
    if (grid::thick_allcells[nonemptymgi] != 1 && grid::W_allcells[nonemptymgi] == 1) {
      printlnlog(
          "force modelgrid cell {} to grey/LTE thick = 1 for update grid since existing W == 1. (will not have gamma "
          "estimators)",
          mgi);
      grid::thick_allcells[nonemptymgi] = 1;
    }

    printlnlog("mgi {} thick: {} (during grid update)", mgi, grid::thick_allcells[nonemptymgi]);

    for (int element = 0; element < get_nelements(); element++) {
      calculate_cellpartfuncts(nonemptymgi, element);
    }
    if (!globals::simulation_continued_from_saved) {
      calculate_ion_balance_nne(nonemptymgi);
    }
  } else {
    // For all other timesteps temperature corrections have to be applied

    // we have to calculate the electron density
    // and all the level populations
    // Normalise estimators and make sure that they are finite. Then update T_R and W using the estimators.
    // (This could in principle also be done for empty cells)

    const auto sys_time_start_temperature_corrections = std::chrono::steady_clock::now();

    radfield::normalise_J(nonemptymgi, estimator_normfactor_over4pi);  // this applies normalisation to the fullspec J
    // this stores the factor that will be applied later for the J bins but not fullspec J
    radfield::set_J_normfactor(nonemptymgi, estimator_normfactor_over4pi);

#ifdef DO_TITER
    radfield::titer_J(nonemptymgi);
#endif

    // lte_iteration really means either ts 0 or nts < globals::num_lte_timesteps
    if (globals::lte_iteration || grid::thick_allcells[nonemptymgi] == 1) {
      // LTE mode or grey mode (where temperature doesn't matter but is calculated anyway)

      const auto T_J = radfield::get_T_J_from_J(nonemptymgi);
      grid::TR_allcells[nonemptymgi] = T_J;
      grid::Te_allcells[nonemptymgi] = T_J;
      grid::TJ_allcells[nonemptymgi] = T_J;
      grid::W_allcells[nonemptymgi] = 1;

      if constexpr (USE_LUT_PHOTOION) {
        std::ranges::fill(
            globals::corrphotoionrenorm.subspan(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground,
                                                globals::nbfcontinua_ground),
            1.);
      }

      for (int element = 0; element < get_nelements(); element++) {
        calculate_cellpartfuncts(nonemptymgi, element);
      }
      calculate_ion_balance_nne(nonemptymgi);
    } else {
      // not lte_iteration and not a thick cell
      // non-pure-LTE timesteps with T_e from heating/cooling

      radfield::normalise_nuJ(nonemptymgi, estimator_normfactor_over4pi);

      globals::ffheatingestimator[nonemptymgi] *= estimator_normfactor;
      if constexpr (!DIRECT_COL_HEAT) {
        globals::colheatingestimator[nonemptymgi] *= estimator_normfactor;
      }

#ifdef DO_TITER
      radfield::titer_nuJ(nonemptymgi);
      titer_average_estimators(nonemptymgi);
#endif

      update_gamma_corrphotoionrenorm_bfheating_estimators(nonemptymgi, estimator_normfactor);

      // Get radiation field parameters (T_J, T_R, W, and bins if enabled) out of the
      // full-spectrum and binned J and nuJ estimators
      radfield::fit_parameters(nonemptymgi, nts);

      solve_Te_nltepops(nonemptymgi, nts, nts_prev, heatingcoolingrates);
    }
    const auto temperature_corrections_duration =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_temperature_corrections)
            .count();
    printlnlog("Temperature/NLTE solution for cell {} timestep {} took {:.1f} seconds", mgi, nts,
               temperature_corrections_duration);
  }

  // Check the temperature that the cell has ended up with for this timestep. This is done here rather
  // than where T_e is stored, because the T_e solver also stores the intermediate guesses that it
  // rejects, and those would be reported too.
  const auto T_e_cell = grid::Te_allcells[nonemptymgi];
  if (T_e_cell > 0.) {
    const double nu_peak_planck = 5.879e10 * T_e_cell;  // Wien's displacement law
    if (nu_peak_planck > NU_MAX_R || nu_peak_planck < NU_MIN_R) {
      printlnlog(
          "[warning] cell {} timestep {}: B_planck(T_e={:g} [K]) peak at {:g} [Hz] is outside the radiation field "
          "frequency range [{:g}, {:g}] [Hz]",
          mgi, nts, T_e_cell, nu_peak_planck, NU_MIN_R, NU_MAX_R);
    }
  }

  const auto nne = grid::get_nne(nonemptymgi);
  const auto propcell_width = grid::propcell_width_tmin(
      mgi, 0);  // pass mgi here because cellindex = mgi for direct map 1D or 2D, and doesn't matter for 3D grid
  const double compton_optical_depth_across_cell = SIGMA_T * nne * propcell_width * tratmid;

  if constexpr (RPKT_GREY_TYPE == RpktGreyType::JUST2022_TEMP_LANTHANIDEFRAC) {
    grid::kappagrey_allcells[nonemptymgi] = grid::calculate_cell_kappagrey(nonemptymgi);
  }
  const double radial_pos = grid::get_modelcell_mean_radial_pos_tmin(mgi) * tratmid;
  const double grey_optical_depth_across_cell =
      grid::kappagrey_allcells[nonemptymgi] * grid::get_rho(nonemptymgi) * propcell_width * tratmid;
  // cube corners will have radial pos > rmax, so clamp to 0.
  const double dist_to_obs = std::max(0., (globals::rmax * tratmid) - radial_pos);
  const auto grey_optical_depth =
      static_cast<float>(grid::kappagrey_allcells[nonemptymgi] * grid::get_rho(nonemptymgi) * dist_to_obs);
  printlnlog("modelgridindex {}, compton optical depth (/propgridcell) {:g}, grey optical depth (/propgridcell) {:g}",
             mgi, compton_optical_depth_across_cell, grey_optical_depth_across_cell);
  printlnlog("radial_pos {:g}, distance_to_obs {:g}, tau_dist {:g}", radial_pos, dist_to_obs, grey_optical_depth);

  grid::grey_depth_allcells[nonemptymgi] = grey_optical_depth;

  if ((grey_optical_depth >= globals::optical_depth_is_thick) && (nts < globals::num_grey_timesteps)) {
    printlnlog("timestep {} cell {} is treated in grey approximation (chi_grey {:g} [cm2/g], tau {:g} >= {:g})", nts,
               mgi, grid::kappagrey_allcells[nonemptymgi], grey_optical_depth, globals::optical_depth_is_thick);
    grid::thick_allcells[nonemptymgi] = 1;
  } else if (VPKT_ON && (grey_optical_depth > vpkt::optical_depth_is_thick_vpkt)) {
    grid::thick_allcells[nonemptymgi] = 2;
  } else {
    grid::thick_allcells[nonemptymgi] = 0;
  }

  if (grid::thick_allcells[nonemptymgi] == 1) {
    // cooling rates calculation can be skipped for thick cells
    // flag with negative numbers to indicate that the rates are invalid
    const auto ioncooling_contribs = kpkt::get_cell_ion_cooling_contribs(nonemptymgi);
    std::ranges::fill(ioncooling_contribs, -1.);
  } else if (globals::simulation_continued_from_saved && nts == globals::timestep_initial) {
    // cooling rates were read from the gridsave file for this timestep
    // make sure they are valid
    const auto ioncooling_contribs = kpkt::get_cell_ion_cooling_contribs(nonemptymgi);
    assert_always(ioncooling_contribs[0] >= 0.);
    assert_always(ioncooling_contribs.back() >= 0.);
  } else {
    // Cooling rates depend only on cell properties, precalculate total cooling
    // and ion contributions inside update grid and communicate between MPI tasks
    const auto sys_time_start_calc_kpkt_rates = std::chrono::steady_clock::now();

    printlog("calculating cooling_rates for timestep {} cell {}...", nts, mgi);

    // don't pass pointer to heatingcoolingrates because current populations and rates weren't used to determine T_e
    kpkt::calculate_cooling_rates(nonemptymgi, nullptr);

    const auto calc_kpkt_rates_duration =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_calc_kpkt_rates).count();
    printlnlog("took {:.1f} seconds", calc_kpkt_rates_duration);
  }

  if constexpr (RPKT_USE_EXPANSION_OPACITIES || VPKT_USE_EXPANSION_OPACITIES ||
                RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value()) {
    if (grid::thick_allcells[nonemptymgi] != 1) {
      calculate_expansion_opacities(nonemptymgi);
    }
  }

  const auto update_grid_cell_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_update_cell).count();
  if (update_grid_cell_seconds > 1.) {
    printlnlog("update_grid_cell for cell {} timestep {} took {:.1f} seconds", mgi, nts, update_grid_cell_seconds);
  }
}

}  // anonymous namespace

//  update the matter quantities in the grid cells at the start of the new timestep.
void update_grid(std::ostream& estimators_file, const int nts, const int nts_prev, const int titer,
                 const std::chrono::steady_clock::time_point real_time_start) {
  const auto my_rank = globals::my_rank;
  const auto sys_time_start_update_grid = std::chrono::steady_clock::now();
  const auto startup_elapsed_seconds =
      std::chrono::duration<double>(sys_time_start_update_grid - real_time_start).count();
  const auto tmid = globals::timesteps[nts].mid;
  printlnlog("timestep {}: time before update grid (tstart + {:.1f} seconds) simtime ts_mid {:g} [days]", nts,
             startup_elapsed_seconds, tmid / DAY);

  globals::lte_iteration = (nts < globals::num_lte_timesteps);
  printlnlog("lte_iteration {}", globals::lte_iteration ? 1 : 0);
  assert_always(globals::num_lte_timesteps > 0);  // The first time step must solve the ionisation balance in LTE

  if (NT_ON && NT_SOLVE_SPENCERFANO && (nts < globals::num_lte_timesteps + 1)) {
    printlnlog("timestep {}: skipping Spencer-Fano solutions for all cells (solver starts at timestep {})", nts,
               globals::num_lte_timesteps + 1);
  }

  const double tratmid = tmid / globals::tmin;

  // These values will not be used if nts == 0, but set them anyway
  // nts_prev is the previous timestep, unless this is timestep zero
  const double deltat = globals::timesteps[nts_prev].width;

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    radfield::normalise_bf_estimators(nts, nts_prev, titer, deltat);
  }

  const auto ndo_nonempty = grid::get_ndo_nonempty(my_rank);
  heatingcoolingrates_thisrankcells.resize(ndo_nonempty);
  std::ranges::fill(heatingcoolingrates_thisrankcells, HeatingCoolingRates{});

  const auto nstart_nonempty = grid::get_nstart_nonempty(my_rank);

  // the decay chain calculations depend only on the timestep midpoint, so the per-nuclide mass fraction
  // coefficients are computed once and reused for the abundance update, the analytic emission rates of the
  // deposition calculation, and the estimator output
  const auto nuc_massfrac_coeffs = decay::calc_nuc_massfrac_coeffs(globals::timesteps[nts].mid);
  const auto emissionratecoeffs = decay::calc_ana_emission_ratecoeffs(nuc_massfrac_coeffs);

  // for this rank's cells: update the elemental abundances (with radioactive decays to the timestep midpoint),
  // then the mass densities and the total ion number densities that depend on both
  decay::update_abundances(nstart_nonempty, ndo_nonempty, nuc_massfrac_coeffs);
  for (int nonemptymgi = nstart_nonempty; nonemptymgi < (nstart_nonempty + ndo_nonempty); nonemptymgi++) {
    const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    grid::set_rho(nonemptymgi, static_cast<float>(grid::get_rho_tmin(mgi) / pow3(tratmid)));
    grid::set_nnetot(nonemptymgi);
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int nonemptymgi = nstart_nonempty; nonemptymgi < (nstart_nonempty + ndo_nonempty); nonemptymgi++) {
    update_grid_cell(nonemptymgi, nts, nts_prev, titer, tratmid, deltat,
                     heatingcoolingrates_thisrankcells.at(nonemptymgi - nstart_nonempty), emissionratecoeffs);
  }

  // serial output of estimator data to this ranks estimator file cell by cell
  const auto nstart = grid::get_nstart(my_rank);
  const auto ndo = grid::get_ndo(my_rank);
  if (ndo > 0) {
    printlnlog("writing to estimators file timestep {} mgi {} to {} (inclusive)", nts, nstart, nstart + ndo - 1);
    for (int mgi = nstart; mgi < (nstart + ndo); mgi++) {
      const auto nonemptymgi = (grid::get_numpropcells(mgi) > 0) ? grid::get_nonemptymgi_of_mgi(mgi) : -1;
      if (nonemptymgi >= 0) {
        assert_always(nonemptymgi >= nstart_nonempty);
        assert_always(nonemptymgi < (nstart_nonempty + ndo_nonempty));
        write_to_estimators_file(estimators_file, nonemptymgi, nts, titer,
                                 heatingcoolingrates_thisrankcells.at(nonemptymgi - nstart_nonempty),
                                 nuc_massfrac_coeffs);
      } else {
        // modelgrid cells that are not represented in the simulation grid
        std::println(estimators_file, "timestep {} modelgridindex {} EMPTYCELL", nts, mgi);
      }

      std::println(estimators_file);
      estimators_file.flush();
    }
  }

  globals::max_path_step = std::min(1.e35, globals::vmax * tmid / 10.);
  printlnlog("max_path_step {:g}", globals::max_path_step);

  const auto time_update_grid_end_thisrank = std::chrono::steady_clock::now();
  const auto rank_process_time =
      std::chrono::duration<double>(time_update_grid_end_thisrank - sys_time_start_update_grid).count();
  printlnlog("timestep {}: finished update grid for rank {} (took {:.1f} seconds)", nts, my_rank, rank_process_time);

  MPI_Barrier_allranks();
  const auto time_update_grid_end_allranks = std::chrono::steady_clock::now();
  const auto rank_wait_time =
      std::chrono::duration<double>(time_update_grid_end_allranks - time_update_grid_end_thisrank).count();
  const auto rank_total_time =
      std::chrono::duration<double>(time_update_grid_end_allranks - sys_time_start_update_grid).count();
  printlnlog(
      "timestep {}: time after update grid on all processes (rank {} took {:.1f}s, waited {:.1f}s, total {:.1f}s)", nts,
      my_rank, rank_process_time, rank_wait_time, rank_total_time);
}
