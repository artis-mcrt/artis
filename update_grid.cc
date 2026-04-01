#include "update_grid.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <format>
#include <ostream>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

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
                              const HeatingCoolingRates& heatingcoolingrates) {
  // return; disable for better performance (if estimators files are not needed)
  const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);

  const auto sys_time_start_write_estimators = std::time(nullptr);

  const auto T_e = grid::get_Te(nonemptymgi);
  const auto nne = grid::get_nne(nonemptymgi);
  const auto Y_e = grid::get_electronfrac(nonemptymgi);

  estimators_file << "timestep " << timestep << " modelgridindex " << mgi << " titeration " << titer << " TR "
                  << grid::get_TR(nonemptymgi) << " Te " << T_e << " W " << grid::get_W(nonemptymgi) << " TJ "
                  << grid::get_TJ(nonemptymgi) << " grey_depth " << grid::grey_depth_allcells[nonemptymgi] << " thick "
                  << grid::thick_allcells[nonemptymgi] << " nne " << nne << " Ye " << Y_e << " tdays "
                  << std::format("{:7.2f}", globals::timesteps[timestep].mid / DAY) << '\n';

  if (globals::total_nlte_levels > 0) {
    nltepop_write_to_file(nonemptymgi, timestep);
  }

  for (int element = 0; element < get_nelements(); element++) {
    if (grid::get_elem_abundance(nonemptymgi, element) <= 0.) {  // skip elements with no abundance
      continue;
    }

    estimators_file << std::format("populations        Z={:2d}", get_atomicnumber(element));
    const int nions = get_nions(element);
    if (nions > 0) {
      // add spaces for missing lowest ion stages to match other elements
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
        estimators_file << "              ";
      }
    }
    double elpop = 0.;
    for (int ion = 0; ion < nions; ion++) {
      elpop += get_nnion(nonemptymgi, element, ion);
      estimators_file << std::format("  {}: {:9.3e}", get_ionstage(element, ion), get_nnion(nonemptymgi, element, ion));
    }
    if (nions == 0) {
      elpop = grid::get_elem_numberdens(nonemptymgi, element);
    }
    estimators_file << std::format("  SUM: {:9.3e}", elpop);

    decay::output_nuc_abundances(estimators_file, nonemptymgi, globals::timesteps[timestep].mid, element);

    if (nions == 0 || elpop <= 0.) {
      // dummy element for nuclear abundances only
      continue;
    }

    estimators_file << std::format("gamma_R            Z={:2d}", get_atomicnumber(element));
    for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
      estimators_file << "              ";
    }
    for (int ion = 0; ion < nions - 1; ion++) {
      estimators_file << std::format("  {}: {:9.3e}", get_ionstage(element, ion),
                                     calculate_iongamma_per_ionpop(nonemptymgi, element, ion, false, false));
    }
    estimators_file << '\n';

    // if we have bound-free estimators, we might want to compare how gamma_R would be if we used the radiation field
    // model instead
    if (DETAILED_BF_ESTIMATORS_ON) {
      estimators_file << std::format("gamma_R_integral   Z={:2d}", get_atomicnumber(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
        estimators_file << "              ";
      }
      for (int ion = 0; ion < nions - 1; ion++) {
        estimators_file << std::format("  {}: {:9.3e}", get_ionstage(element, ion),
                                       calculate_iongamma_per_ionpop(nonemptymgi, element, ion, false, true));
      }
      estimators_file << '\n';
    }

    if (NT_ON) {
      estimators_file << std::format("gamma_NT           Z={:2d}", get_atomicnumber(element));
      for (int ionstage = 1; ionstage < get_ionstage(element, 0); ionstage++) {
        estimators_file << "              ";
      }
      for (int ion = 0; ion < nions - 1; ion++) {
        const double Y_nt = nonthermal::nt_ionisation_ratecoeff(nonemptymgi, element, ion);
        estimators_file << std::format("  {}: {:9.3e}", get_ionstage(element, ion), Y_nt);
      }
      estimators_file << '\n';
    }

    if (USE_LUT_PHOTOION && globals::nbfcontinua_ground > 0) {
      estimators_file << std::format("corrphotoionrenorm Z={:2d}", get_atomicnumber(element));
      for (int ion = 0; ion < nions - 1; ion++) {
        if (get_groundcontindex(element, ion) >= 0) {
          estimators_file << std::format(
              "  {}: {:9.3e}", get_ionstage(element, ion),
              globals::corrphotoionrenorm[get_ionestimindex_nonemptymgi(nonemptymgi, element, ion)]);
        }
      }
      estimators_file << '\n';
      estimators_file << std::format("gammaestimator     Z={:2d}", get_atomicnumber(element));
      for (int ion = 0; ion < nions - 1; ion++) {
        if (get_groundcontindex(element, ion) >= 0) {
          estimators_file << std::format(
              "  {}: {:9.3e}", get_ionstage(element, ion),
              globals::gammaestimator[get_ionestimindex_nonemptymgi(nonemptymgi, element, ion)]);
        }
      }
      estimators_file << '\n';
    }
  }

  // power densities in erg / s / cm^3. 'ana' means analytical at t_mid, i.e. the rates calculated from the nuclear
  // abundances and decay data, not from Monte Carlo events
  estimators_file << std::format("emission_ana: gamma {:11.5e} positron {:11.5e} electron {:11.5e} alpha {:11.5e}",
                                 heatingcoolingrates.eps_gamma_ana, heatingcoolingrates.eps_positron_ana,
                                 heatingcoolingrates.eps_electron_ana, heatingcoolingrates.eps_alpha_ana);
  const bool any_fission = decay::decaytype_is_used(decay::DECAYTYPE_SPONTFISSION);
  if (any_fission) {
    estimators_file << std::format(" spfission {:11.5e}", heatingcoolingrates.eps_spfission_ana);
  }
  estimators_file << '\n';

  estimators_file << std::format("deposition: gamma {:11.5e} positron {:11.5e} electron {:11.5e} alpha {:11.5e}",
                                 heatingcoolingrates.dep_gamma, heatingcoolingrates.dep_positron,
                                 heatingcoolingrates.dep_electron, heatingcoolingrates.dep_alpha);
  if (any_fission) {
    estimators_file << std::format(" spfission {:11.5e}", heatingcoolingrates.dep_spfission);
  }
  estimators_file << '\n';

  estimators_file << std::format(
      "heating: ff {:11.5e} bf {:11.5e} coll {:11.5e}       dep {:11.5e} heating_dep/total_dep {:.3f}\n",
      heatingcoolingrates.heating_ff, heatingcoolingrates.heating_bf, heatingcoolingrates.heating_collisional,
      heatingcoolingrates.heating_dep, heatingcoolingrates.dep_frac_heating);
  estimators_file << std::format("cooling: ff {:11.5e} fb {:11.5e} coll {:11.5e} adiabatic {:11.5e}\n",
                                 heatingcoolingrates.cooling_ff, heatingcoolingrates.cooling_fb,
                                 heatingcoolingrates.cooling_collisional, heatingcoolingrates.cooling_adiabatic);

  estimators_file << '\n';

  estimators_file.flush();

  const auto write_estim_duration = std::time(nullptr) - sys_time_start_write_estimators;
  if (write_estim_duration >= 1) {
    printlnlog("writing estimators for timestep {} cell {} took {} seconds", timestep, mgi, write_estim_duration);
  }
}

void solve_Te_nltepops(const int nonemptymgi, const int nts, const int nts_prev,
                       HeatingCoolingRates& heatingcoolingrates) {
  const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  // bfheating coefficients are needed for the T_e solver, but they only depend on the radiation field, which is fixed
  // during the iterations below
  printlog("calculate_bfheatingcoeffs for timestep {} cell {}...", nts, mgi);
  const auto sys_time_start_calculate_bfheatingcoeffs = std::time(nullptr);
  THREADLOCALONHOST auto bfheatingcoeffs = std::vector<double>(get_includedlevels());

  calculate_bfheatingcoeffs(nonemptymgi, bfheatingcoeffs);
  printlnlog("took {} seconds", std::time(nullptr) - sys_time_start_calculate_bfheatingcoeffs);

  constexpr double nne_reltol = 0.04;
  constexpr double T_e_reltol = 0.04;
  for (int nlte_iter = 0; nlte_iter <= NLTEITER; nlte_iter++) {
    const auto sys_time_start_spencerfano = std::time(nullptr);
    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      // SF solution depends on the ionisation balance, and weakly on nne
      nonthermal::solve_spencerfano(nonemptymgi, nts, nlte_iter);
    }
    const auto duration_solve_spencerfano = std::time(nullptr) - sys_time_start_spencerfano;

    const auto sys_time_start_partfuncs_or_gamma = std::time(nullptr);
    for (int element = 0; element < get_nelements(); element++) {
      if (!elem_has_nlte_levels(element)) {
        calculate_cellpartfuncts(nonemptymgi, element);
      }
    }
    const auto duration_solve_partfuncs_or_gamma = std::time(nullptr) - sys_time_start_partfuncs_or_gamma;

    const double prev_T_e = grid::get_Te(nonemptymgi);
    const auto sys_time_start_Te = std::time(nullptr);

    // Find T_e as solution for thermal balance
    call_T_e_finder(nonemptymgi, globals::timesteps[nts_prev].mid, MINTEMP, MAXTEMP, heatingcoolingrates,
                    bfheatingcoeffs);

    const auto duration_solve_T_e = std::time(nullptr) - sys_time_start_Te;

    if (globals::total_nlte_levels == 0) {
      const auto sys_time_start_pops = std::time(nullptr);
      calculate_ion_balance_nne(nonemptymgi);
      const auto duration_solve_pops = std::time(nullptr) - sys_time_start_pops;

      printlnlog(
          "Grid solver cell {} timestep {}: time spent on: Spencer-Fano {}s, partfuncs/gamma {}s, T_e {}s, populations "
          "{}s",
          mgi, nts, duration_solve_spencerfano, duration_solve_partfuncs_or_gamma, duration_solve_T_e,
          duration_solve_pops);
      break;  // no iteration is needed without nlte pops
    }

    const double fracdiff_T_e = fabs((grid::get_Te(nonemptymgi) / prev_T_e) - 1);
    const auto sys_time_start_nltepops = std::time(nullptr);
    // fractional difference between previous and current iteration's (nne or max(ground state
    // population change))
    for (int element = 0; element < get_nelements(); element++) {
      if (get_nions(element) > 0 && elem_has_nlte_levels(element)) {
        solve_nlte_pops_element(element, nonemptymgi, nts, nlte_iter);
        calculate_cellpartfuncts(nonemptymgi, element);
      }
    }
    const auto duration_solve_nltepops = std::time(nullptr) - sys_time_start_nltepops;

    const double nne_prev = grid::get_nne(nonemptymgi);
    calculate_ion_balance_nne(nonemptymgi);  // sets nne
    const auto fracdiff_nne = fabs((grid::get_nne(nonemptymgi) / nne_prev) - 1);
    printlnlog(
        "NLTE solver mgi {} timestep {} iteration {}: time spent on: Spencer-Fano {}s, T_e {}s, NLTE populations {}s",
        mgi, nts, nlte_iter, duration_solve_spencerfano, duration_solve_T_e, duration_solve_nltepops);
    printlnlog(
        "NLTE (Spencer-Fano/Te/pops) solver mgi {} timestep {} iteration {}: prev_iter nne {:e}, new nne is {:e}, "
        "fracdiff {:g}, prev T_e {:g} new T_e {:g} fracdiff {:g}",
        mgi, nts, nlte_iter, nne_prev, grid::get_nne(nonemptymgi), fracdiff_nne, prev_T_e, grid::get_Te(nonemptymgi),
        fracdiff_T_e);
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
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
        const int groundcontindex = get_groundcontindex(element, ion);
        if (groundcontindex < 0) {
          continue;
        }
        const int ionestimindex = (nonemptymgi * globals::nbfcontinua_ground) + groundcontindex;

        globals::gammaestimator[ionestimindex] *= estimator_normfactor / H;
#ifdef DO_TITER
        if (globals::gammaestimator_save[ionestimindex] >= 0) {
          globals::gammaestimator[ionestimindex] =
              (globals::gammaestimator[ionestimindex] + globals::gammaestimator_save[ionestimindex]) / 2;
        }
        globals::gammaestimator_save[ionestimindex] = globals::gammaestimator[ionestimindex];
#endif

        globals::corrphotoionrenorm[ionestimindex] =
            globals::gammaestimator[ionestimindex] / get_corrphotoioncoeff_ana(element, ion, 0, 0, nonemptymgi);

        assert_always(std::isfinite(globals::corrphotoionrenorm[ionestimindex]));
      }

      // 2012-01-11. These loops should terminate here to precalculate *ALL* corrphotoionrenorm
      // values so that the values are known when required by the call to get_corrphotoioncoeff in
      // the following loops. Otherwise get_corrphotoioncoeff tries to renormalize by the closest
      // corrphotoionrenorm in frequency space which can lead to zero contributions to the total
      // photoionsation rate!
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
      const int ionestimindex = (nonemptymgi * globals::nbfcontinua_ground) + groundcontindex;

      if (!elem_has_nlte_levels(element)) {
        globals::gammaestimator[ionestimindex] = calculate_iongamma_per_gspop(nonemptymgi, element, ion);
      }

      if constexpr (USE_ION_BFHEATING_ESTIMATORS) {
        globals::bfheatingestimator[ionestimindex] *= estimator_normfactor;
#ifdef DO_TITER
        if (globals::bfheatingestimator_save[ionestimindex] >= 0) {
          globals::bfheatingestimator[ionestimindex] =
              (globals::bfheatingestimator[ionestimindex] + globals::bfheatingestimator_save[ionestimindex]) / 2;
        }
        globals::bfheatingestimator_save[ionestimindex] = globals::bfheatingestimator[ionestimindex];
#endif
        // Now convert bfheatingestimator into the bfheating renormalisation coefficient used in
        // for the remaining part of update_grid. At the start of the next update_packets, it will be reset

        const double bfheatingcoeff_ground = calculate_bfheatingcoeff(element, ion, 0, 0, nonemptymgi);
        globals::bfheatingestimator[ionestimindex] /= bfheatingcoeff_ground;

        assert_always(std::isfinite(globals::bfheatingestimator[ionestimindex]));
      }
    }
  }
}

#ifdef DO_TITER
static void titer_average_estimators(const int nonemptymgi) {
  if (globals::ffheatingestimator_save[nonemptymgi] >= 0) {
    globals::ffheatingestimator[nonemptymgi] =
        (globals::ffheatingestimator[nonemptymgi] + globals::ffheatingestimator_save[nonemptymgi]) / 2;
  }
  globals::ffheatingestimator_save[nonemptymgi] = globals::ffheatingestimator[nonemptymgi];
  if constexpr (!DIRECT_COL_HEAT) {
    if (globals::colheatingestimator_save[nonemptymgi] >= 0) {
      globals::colheatingestimator[nonemptymgi] =
          (globals::colheatingestimator[nonemptymgi] + globals::colheatingestimator_save[nonemptymgi]) / 2;
    }
    globals::colheatingestimator_save.at(nonemptymgi) = globals::colheatingestimator.at(nonemptymgi);
  }
}
#endif

void update_grid_cell(const int nonemptymgi, const int nts, const int nts_prev, const int titer, const double tratmid,
                      const double deltat, HeatingCoolingRates& heatingcoolingrates) {
  const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);

  const double deltaV =
      grid::get_modelcell_assocvolume_tmin(mgi) * pow(globals::timesteps[nts_prev].mid / globals::tmin, 3);
  const auto sys_time_start_update_cell = std::time(nullptr);

  printlnlog("update_grid_cell: working on mgi {} before timestep {} titeration {}...", mgi, nts, titer);

  // Update current mass density of cell
  const auto rho = static_cast<float>(grid::get_rho_tmin(mgi) / pow(tratmid, 3));
  grid::set_rho(nonemptymgi, rho);

  // Update elemental abundances with radioactive decays
  decay::update_abundances(nonemptymgi, globals::timesteps[nts].mid);
  nonthermal::calculate_deposition_rate_density(nonemptymgi, nts, heatingcoolingrates);

  const double estimator_normfactor = 1 / deltaV / deltat / globals::nprocs;
  const double estimator_normfactor_over4pi = ONEOVER4PI * estimator_normfactor;

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
    // last timestep. Therefore it has no valid Gamma estimators and must
    // be treated in LTE at restart.
    if (grid::thick_allcells[nonemptymgi] != 1 && grid::get_W(nonemptymgi) == 1) {
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
    // Normalise estimators and make sure that they are finite.
    // Then update T_R and W using the estimators.
    // (This could in principle also be done for empty cells)

    const auto sys_time_start_temperature_corrections = std::time(nullptr);

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
      grid::set_TR(nonemptymgi, T_J);
      grid::set_Te(nonemptymgi, T_J);
      grid::set_TJ(nonemptymgi, T_J);
      grid::set_W(nonemptymgi, 1);

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
    printlnlog("Temperature/NLTE solution for cell {} timestep {} took {} seconds", mgi, nts,
               std::time(nullptr) - sys_time_start_temperature_corrections);
  }

  const auto nne = grid::get_nne(nonemptymgi);
  const double compton_optical_depth_across_cell = SIGMA_T * nne * grid::propcell_width_tmin(mgi, 0) * tratmid;

  if constexpr (RPKT_GREY_TYPE == RpktGreyType::JUST2022_TEMP_LANTHANIDEFRAC) {
    grid::set_kappagrey(nonemptymgi, grid::calculate_cell_kappagrey(nonemptymgi));
  }
  const double radial_pos = grid::get_modelcell_mean_radial_pos(mgi, tratmid);
  const double grey_optical_depth_across_cell =
      grid::get_kappagrey(nonemptymgi) * grid::get_rho(nonemptymgi) * grid::propcell_width_tmin(mgi, 0) * tratmid;
  // cube corners will have radial pos > rmax, so clamp to 0.
  const double dist_to_obs = std::max(0., (globals::rmax * tratmid) - radial_pos);
  const auto grey_optical_depth =
      static_cast<float>(grid::get_kappagrey(nonemptymgi) * grid::get_rho(nonemptymgi) * dist_to_obs);
  printlnlog("modelgridindex {}, compton optical depth (/propgridcell) {:g}, grey optical depth (/propgridcell) {:g}",
             mgi, compton_optical_depth_across_cell, grey_optical_depth_across_cell);
  printlnlog("radial_pos {:g}, distance_to_obs {:g}, tau_dist {:g}", radial_pos, dist_to_obs, grey_optical_depth);

  grid::grey_depth_allcells[nonemptymgi] = grey_optical_depth;

  // grey_optical_depth = compton_optical_depth;

  if ((grey_optical_depth >= globals::cell_is_optically_thick) && (nts < globals::num_grey_timesteps)) {
    printlnlog("timestep {} cell {} is treated in grey approximation (chi_grey {:g} [cm2/g], tau {:g} >= {:g})", nts,
               mgi, grid::get_kappagrey(nonemptymgi), grey_optical_depth, globals::cell_is_optically_thick);
    grid::thick_allcells[nonemptymgi] = 1;
  } else if (VPKT_ON && (grey_optical_depth > vpkt::cell_is_optically_thick_vpkt)) {
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
    const auto sys_time_start_calc_kpkt_rates = std::time(nullptr);

    printlog("calculating cooling_rates for timestep {} cell {}...", nts, mgi);

    // don't pass pointer to heatingcoolingrates because current populations and rates weren't
    // used to determine T_e
    kpkt::calculate_cooling_rates(nonemptymgi, nullptr);

    printlnlog("took {} seconds", std::time(nullptr) - sys_time_start_calc_kpkt_rates);
  }

  if constexpr (EXPANSIONOPACITIES_ON || RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value()) {
    if (grid::thick_allcells[nonemptymgi] != 1) {
      calculate_expansion_opacities(nonemptymgi);
    }
  }

  const auto update_grid_cell_seconds = std::time(nullptr) - sys_time_start_update_cell;
  if (update_grid_cell_seconds > 0) {
    printlnlog("update_grid_cell for cell {} timestep {} took {} seconds", mgi, nts, update_grid_cell_seconds);
  }
}

}  // anonymous namespace

//  update the matter quantities in the grid cells at the start of the new timestep.
void update_grid(std::ostream& estimators_file, const int nts, const int nts_prev, const int titer,
                 const std::time_t real_time_start) {
  const auto my_rank = globals::my_rank;
  const auto sys_time_start_update_grid = std::time(nullptr);

  printlnlog("timestep {}: time before update grid (tstartup + {} seconds) simtime ts_mid {:g} days", nts,
             sys_time_start_update_grid - real_time_start, globals::timesteps[nts].mid / DAY);

  globals::lte_iteration = (globals::timestep < globals::num_lte_timesteps);
  printlnlog("lte_iteration {}", globals::lte_iteration ? 1 : 0);
  assert_always(globals::num_lte_timesteps > 0);  // The first time step must solve the ionisation balance in LTE

  const double tratmid = globals::timesteps[nts].mid / globals::tmin;

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

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int nonemptymgi = nstart_nonempty; nonemptymgi < (nstart_nonempty + ndo_nonempty); nonemptymgi++) {
    update_grid_cell(nonemptymgi, nts, nts_prev, titer, tratmid, deltat,
                     heatingcoolingrates_thisrankcells.at(nonemptymgi - nstart_nonempty));
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
                                 heatingcoolingrates_thisrankcells.at(nonemptymgi - nstart_nonempty));
      } else {
        // modelgrid cells that are not represented in the simulation grid
        estimators_file << std::format("timestep {} modelgridindex {} EMPTYCELL\n\n", nts, mgi);
        estimators_file.flush();
      }
    }
  }

  globals::max_path_step = std::min(1.e35, globals::rmax / 10.);
  printlnlog("max_path_step {:g}", globals::max_path_step);

  const auto time_update_grid_end_thisrank = std::time(nullptr);
  printlnlog("timestep {}: finished update grid for rank {} (took {} seconds)", nts, my_rank,
             time_update_grid_end_thisrank - sys_time_start_update_grid);

  MPI_Barrier(MPI_COMM_WORLD);
  printlnlog("timestep {}: time after update grid on all processes (rank {} took {}, waited {}, total {} seconds)", nts,
             my_rank, time_update_grid_end_thisrank - sys_time_start_update_grid,
             std::time(nullptr) - time_update_grid_end_thisrank, std::time(nullptr) - sys_time_start_update_grid);
}
