// Level populations and ionisation balance in LTE and approximate NLTE: partition functions,
// Boltzmann/Saha level and ion populations, and the solver for a self-consistent free electron
// density (nne).

#include "ltepop.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <tuple>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <boost/math/tools/toms748_solve.hpp>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "mpi_logging.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "ratecoeff.h"
#include "sn3d.h"

namespace {

// The conditions behind the warnings below can persist across the many evaluations that the nne and
// T_e root solvers make for a single cell, so each one is reported only on its first occurrence for
// a given cell and timestep. The check is made at the warning site so that the warnings are still
// reported when the NLTE solver calls into these functions itself.
struct CellWarningMarker {
  ptrdiff_t nonemptymgi{-1};
  int timestep{-1};

  // returns true the first time it is called for a particular cell and timestep
  auto is_first_occurrence(const ptrdiff_t nonemptymgi_in) -> bool {
    if (nonemptymgi_in == nonemptymgi && globals::timestep == timestep) {
      return false;
    }
    nonemptymgi = nonemptymgi_in;
    timestep = globals::timestep;
    return true;
  }
};

THREADLOCALONHOST CellWarningMarker neutralcell_warned;
THREADLOCALONHOST CellWarningMarker nne_maxit_warned;
THREADLOCALONHOST CellWarningMarker uppermost_ion_warned;
THREADLOCALONHOST CellWarningMarker ionfract_zeroed_warned;

// use Saha equation for LTE ionisation balance
[[gnu::pure]] [[nodiscard]] auto phi_saha(const int element, const int ion, const int nonemptymgi) -> double {
  const auto partfunc_ion = get_ion_partfunct(nonemptymgi, element, ion);
  const auto partfunc_upperion = get_ion_partfunct(nonemptymgi, element, ion + 1);

  const auto T_e = grid::Te_allcells[nonemptymgi];
  const double ionpot = epsilon(element, ion + 1, 0) - epsilon(element, ion, 0);
  const double partfunct_ratio = partfunc_ion / partfunc_upperion;
  return partfunct_ratio * SAHACONST * pow(T_e, -1.5) * exp(ionpot / KB / T_e);
}

// Calculate population ratio (a saha factor) of two consecutive ionisation stages in nebular approximation phi_j,k* =
// N_j,k*/(N_j+1,k* * nne)
// should have units of [cm^3] so that when multiplied by nne it gives the population ratio of two consecutive
// ionisation stages
[[gnu::pure]] [[nodiscard]] auto phi_rate_balance(const int element, const int ion, const int nonemptymgi) -> double {
  testmodeassert_valid_ion(element, ion);

  assert_testmodeonly(!globals::lte_iteration);
  assert_testmodeonly(grid::thick_allcells[nonemptymgi] != 1);  // should use use phi_lte instead

  assert_testmodeonly(!elem_has_nlte_levels(element));  // don't use this function if the NLTE solver is active

  const int uniqueionindex = get_uniqueionindex(element, ion);
  const auto partfunc_ion = get_ion_partfunct(nonemptymgi, element, ion);

  const auto T_e = grid::Te_allcells[nonemptymgi];

  // photoionisation plus collisional ionisation rate coefficient per ground level pop
  const auto groundcontindex = get_groundcontindex(element, ion);
  const double Gamma_groundlevel =
      groundcontindex >= 0
          ? globals::gammaestimator[(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) +
                                    groundcontindex]
          : 0.;

  // Convert Gamma to the photoionisation rate per ion pop
  const double Gamma_ion = Gamma_groundlevel * stat_weight(element, ion, 0) / partfunc_ion;

  const double Alpha_sp = get_ion_spontrecombcoeff(uniqueionindex, T_e);
  constexpr bool include_collisional_recombination = false;
  const double Col_rec = include_collisional_recombination
                             ? calculate_ionrecombcoeff(nonemptymgi, T_e, element, ion + 1, false, true, false, false)
                             : 0.;

  const double gamma_nt = NT_ON ? nonthermal::nt_ionisation_ratecoeff(nonemptymgi, element, ion) : 0.;

  // gamma_nt should generally be higher than the Gamma term for nebular epoch

  assert_always((Gamma_ion + gamma_nt) > 0);
  // numerator: recombination rate coefficient, i.e. rate per upper ion pop per nne [cm^3/s]
  // denominator: ionisation rate per lower ion pop [1/s]
  const double phi = (Alpha_sp + Col_rec) / (Gamma_ion + gamma_nt);
  assert_always(phi > 0.);

  return phi;
}

// calculate the free electron contribution from an element
[[gnu::pure]] [[nodiscard]] auto get_element_nne_contrib(const int nonemptymgi, const int element) -> double {
  if (grid::get_elem_numberdens(nonemptymgi, element) <= 0.) {
    return 0.;
  }

  double nne = 0.;
  const int nions = get_nions(element);
  for (int ion = 0; ion < nions; ion++) {
    const auto nnion = get_nnion(nonemptymgi, element, ion);
    const int ioncharge = get_ionstage(element, ion) - 1;
    nne += ioncharge * nnion;
  }
  return nne;
}

// assume a value for nne and then calculate the resulting nne
// the difference between the assumed and calculated nne is returned
auto nne_solution_f(const double nne_assumed, const int nonemptymgi, const bool force_saha) -> double {
  double nne_after = 0.;  // the resulting nne after setting the ion balance with nne_assumed
  for (int element = 0; element < get_nelements(); element++) {
    const double nnelement = grid::get_elem_numberdens(nonemptymgi, element);
    if (nnelement > 0 && get_nions(element) > 0) {
      if (!force_saha && elem_has_nlte_levels(element)) {
        // populations from the NLTE solver are fixed during the nne solver
        nne_after += get_element_nne_contrib(nonemptymgi, element);
      } else {
        const bool use_phi_saha = force_saha || FORCE_SAHA_ION_BALANCE(get_atomicnumber(element));
        const auto ionfractions = calculate_ionfractions(element, nonemptymgi, nne_assumed, use_phi_saha);
        const int uppermost_ion = static_cast<int>(std::ssize(ionfractions) - 1);
        for (int ion = 0; ion <= uppermost_ion; ion++) {
          const double nnion = nnelement * ionfractions[ion];
          const int ioncharge = get_ionstage(element, ion) - 1;
          nne_after += ioncharge * nnion;
        }
      }

      assert_always(std::isfinite(nne_after));
    }
  }
  nne_after = std::max(MINPOP, nne_after);

  return nne_after - nne_assumed;
}

// return population and whether the population came from the nlte solver
auto calculate_levelpop_nominpop(const int nonemptymgi, const int element, const int ion, const int level)
    -> std::tuple<double, bool> {
  testmodeassert_valid_level(element, ion, level);

  if (level == 0) {
    return {get_groundlevelpop(nonemptymgi, element, ion), false};
  }

  if (elem_has_nlte_levels(element)) {
    if (is_nlte(element, ion, level)) {
      const double nltepop_over_rho = get_nlte_levelpop_over_rho(nonemptymgi, element, ion, level);
      if (nltepop_over_rho >= 0.) {
        const double nn = nltepop_over_rho * grid::get_rho(nonemptymgi);
        return {nn, true};
      }
    } else if (ion_has_superlevel(element, ion)) {
      // level is in the superlevel (or is an autoionising level, whose population is attached
      // to the superlevel outside the NLTE solver). Ions without a superlevel fall through to
      // the Boltzmann fallback below rather than reading a non-existent superlevel slot
      assert_testmodeonly(level_isinsuperlevel(element, ion, level) || level_isautoionising(element, ion, level));

      const double superlevelpop_over_rho = get_nlte_superlevelpop_over_rho_over_slpartfunc(nonemptymgi, element, ion);
      if (superlevelpop_over_rho >= 0.) {
        const double nn = superlevelpop_over_rho * grid::get_rho(nonemptymgi) *
                          superlevel_boltzmann(nonemptymgi, element, ion, level);
        return {nn, true};
      }
    }
  }

  return {calculate_levelpop_boltzmann(nonemptymgi, element, ion, level), false};
}

// Calculate the partition function for ion=ion of element=element in a cell modelgridindex
auto calculate_partfunct(const int element, const int ion, const int nonemptymgi) -> float {
  testmodeassert_valid_ion(element, ion);
  double pop_store{NAN};

  bool initial = false;
  if (get_groundlevelpop(nonemptymgi, element, ion) < MINPOP) {
    // either there really is none of this ion or this is a first pass through
    // in either case, we won't have any real nlte_populations so the actual value of
    // of groundlevelpop for this calculation doesn't matter, so long as it's not zero!
    pop_store = get_groundlevelpop(nonemptymgi, element, ion);
    initial = true;
    set_groundlevelpop(nonemptymgi, element, ion, 1.);
  }

  double U = 1.;

  const int nlevels = get_nlevels(element, ion);
  const double groundpop = get_groundlevelpop(nonemptymgi, element, ion);
  for (int level = 1; level < nlevels; level++) {
    const auto nn = std::get<0>(calculate_levelpop_nominpop(nonemptymgi, element, ion, level));
    U += nn / groundpop;
  }
  U *= stat_weight(element, ion, 0);
  const auto U_float = static_cast<float>(U);
  assert_always(U_float > 0.);

  if (initial) {
    // put back the zero, just in case it matters for something
    set_groundlevelpop(nonemptymgi, element, ion, static_cast<float>(pop_store));
  }

  return U_float;
}

// Set the cell's free electron density nne to the sum of every element's electron contribution (floored at MINPOP).
void set_calculated_nne(const int nonemptymgi) {
  double nne = 0.;  // free electron density
  for (int element = 0; element < get_nelements(); element++) {
    nne += get_element_nne_contrib(nonemptymgi, element);
  }

  grid::set_nne(nonemptymgi, static_cast<float>(std::max(MINPOP, nne)));
}

// Special case of only neutral ions, set nne to some finite value so that packets are not lost in kpkts
void set_groundlevelpops_neutral(const ptrdiff_t nonemptymgi) {
  if (neutralcell_warned.is_first_occurrence(nonemptymgi)) {
    printlnlog("[warning] set_groundlevelpops_neutral: only neutral ions in cell {} timestep {} (repeats suppressed)",
               grid::get_mgi_of_nonemptymgi(nonemptymgi), globals::timestep);
  }
  for (int element = 0; element < get_nelements(); element++) {
    const auto nnelement = grid::get_elem_numberdens(nonemptymgi, element);
    const int nions = get_nions(element);
    // Assign the species population to the neutral ion and set higher ions to MINPOP
    for (int ion = 0; ion < nions; ion++) {
      double nnion{NAN};
      if (ion == 0) {
        nnion = nnelement;
      } else if (nnelement > 0.) {
        nnion = MINPOP;
      } else {
        nnion = 0.;
      }
      const auto groundpop =
          static_cast<float>(nnion * stat_weight(element, ion, 0) / get_ion_partfunct(nonemptymgi, element, ion));

      set_groundlevelpop(nonemptymgi, element, ion, groundpop);
    }
  }
}

// Solve for the free electron density nne in [0, nne_max] that makes the ion balance self-consistent, using the
// Boost TOMS 748 root-finder on the charge-conservation residual.
auto find_converged_nne(const int nonemptymgi, double nne_max, const bool force_lte) -> float {
  // search for nne solution in [0.,nne_max]

  const auto f_nne = [nonemptymgi, force_lte](const double nne) { return nne_solution_f(nne, nonemptymgi, force_lte); };

  constexpr double nne_min = 0.;
  const auto f_nne_min = f_nne(nne_min);
  const auto f_nne_max = f_nne(nne_max);
  assert_always(f_nne_min * f_nne_max <= 0.);  // there should be a root in the interval

  constexpr double fractional_accuracy = 1e-3;
  constexpr auto maxit = 50U;

  // use TOMS 748 solver from Boost
  uintmax_t iter = maxit;
  const auto result =
      boost::math::tools::toms748_solve(f_nne, nne_min, nne_max, f_nne_min, f_nne_max, ftol<fractional_accuracy>, iter);
  const double nne_solution = 0.5 * (result.first + result.second);
  if (iter >= maxit && nne_maxit_warned.is_first_occurrence(nonemptymgi)) {
    printlnlog(
        "[warning] find_converged_nne: cell {} timestep {}: nne did not converge within {} iterations "
        "(repeats suppressed)",
        grid::get_mgi_of_nonemptymgi(nonemptymgi), globals::timestep, iter);
  }

  return static_cast<float>(std::max(MINPOP, nne_solution));
}

}  // anonymous namespace

// Determine the highest ion stage of an element that retains a non-negligible population in a cell (higher ions
// are truncated), using either Saha or photoionisation/recombination rate balance.
[[nodiscard]] auto find_uppermost_ion(const int nonemptymgi, const int element, const double nne_hi,
                                      const bool force_saha) -> int {
  const int nions = get_nions(element);
  if (nions == 0) {
    return -1;
  }
  if (!force_saha && elem_has_nlte_levels(element)) {
    return nions - 1;
  }
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const bool use_phi_saha = force_saha || FORCE_SAHA_ION_BALANCE(get_atomicnumber(element));
  int uppermost_ion = nions - 1;
  if (!use_phi_saha) {
    for (int ion = 0; ion < nions - 1; ion++) {
      if (iongamma_is_zero(nonemptymgi, element, ion) &&
          (!NT_ON || ((globals::dep_estimator_gamma[nonemptymgi] == 0.) &&
                      (grid::get_modelinitnucmassfrac(modelgridindex, decay::get_nucindex(24, 48)) == 0.) &&
                      (grid::get_modelinitnucmassfrac(modelgridindex, decay::get_nucindex(28, 56)) == 0.)))) {
        uppermost_ion = ion;
        break;
      }
    }
  }

  // running product of nne * phi over the ions below, i.e. the population ratio of the ground ion to ion+1.
  // Once it overflows to infinity the higher ions are utterly negligible, so the list is truncated there.
  double pop_ratio_ground_to_upper = 1.;
  for (int ion = 0; ion < uppermost_ion; ion++) {
    const auto phifactor =
        use_phi_saha ? phi_saha(element, ion, nonemptymgi) : phi_rate_balance(element, ion, nonemptymgi);
    pop_ratio_ground_to_upper *= nne_hi * phifactor;

    if (!std::isfinite(pop_ratio_ground_to_upper)) {
      if (uppermost_ion_warned.is_first_occurrence(nonemptymgi)) {
        printlnlog(
            "[info] find_uppermost_ion: cell {} timestep {}: uppermost ion stage limited by phi factor overflow for "
            "Z={} ionstage {} (repeats suppressed)",
            modelgridindex, globals::timestep, get_atomicnumber(element), get_ionstage(element, ion));
      }
      return ion;
    }
  }
  return uppermost_ion;
}

// Calculate the fractions of an element's population in each ionisation stage based on Saha LTE or ionisation
// equilibrium
[[nodiscard]] auto calculate_ionfractions(const int element, const int nonemptymgi, const double nne,
                                          const bool use_phi_saha) -> std::vector<double> {
  assert_testmodeonly(element < get_nelements());
  const int uppermost_ion = grid::get_elements_uppermost_ion(nonemptymgi, element);

  if (uppermost_ion < 0) {
    return {};
  }

  std::vector<double> ionfractions(uppermost_ion + 1);
  ionfractions[uppermost_ion] = 1;

  double normfactor = 1.;

  for (int ion = uppermost_ion - 1; ion >= 0; ion--) {
    const auto phifactor =
        use_phi_saha ? phi_saha(element, ion, nonemptymgi) : phi_rate_balance(element, ion, nonemptymgi);
    ionfractions[ion] = ionfractions[ion + 1] * nne * phifactor;
    normfactor += ionfractions[ion];
  }

  for (int ion = 0; ion <= uppermost_ion; ion++) {
    ionfractions[ion] = ionfractions[ion] / normfactor;

    if (normfactor == 0. || !std::isfinite(ionfractions[ion])) {
      if (ionfract_zeroed_warned.is_first_occurrence(nonemptymgi)) {
        printlnlog(
            "[warning] calculate_ionfractions: cell {} timestep {}: non-finite ionfract set to zero for Z={} "
            "ionstage {} (T_e {:g} [K], T_R {:g} [K]; repeats suppressed)",
            grid::get_mgi_of_nonemptymgi(nonemptymgi), globals::timestep, get_atomicnumber(element),
            get_ionstage(element, ion), grid::Te_allcells[nonemptymgi], grid::TR_allcells[nonemptymgi]);
      }
      ionfractions[ion] = 0;
    }
  }
  return ionfractions;
}

// Calculate occupation population of a level assuming LTE excitation
[[gnu::pure]] [[nodiscard]] auto calculate_levelpop_boltzmann(const int nonemptymgi, const int element, const int ion,
                                                              const int level) -> double {
  testmodeassert_valid_level(element, ion, level);
  const auto nnground = get_groundlevelpop(nonemptymgi, element, ion);
  if (level == 0) {
    return nnground;
  }

  const auto T_exc = LTEPOP_EXCITATION_USE_TJ ? grid::TJ_allcells[nonemptymgi] : grid::Te_allcells[nonemptymgi];
  const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);

  const double E_aboveground = epsilon(ionuniquelevelindexstart + level) - epsilon(ionuniquelevelindexstart);

  return (nnground * stat_weight(ionuniquelevelindexstart + level) / stat_weight(ionuniquelevelindexstart) *
          exp(-E_aboveground / KB / T_exc));
}

// Return the population of a level, applying the MINPOP floor when the element is present (or 0 when it is absent),
// unless the underlying calculation flags that the floor should be skipped.
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto calculate_levelpop(const int nonemptymgi, const int element, const int ion,
                                                                const int level) -> double {
  const auto [nn, skipminpop] = calculate_levelpop_nominpop(nonemptymgi, element, ion, level);
  if (!skipminpop && nn < MINPOP) {
    if (grid::get_elem_massfrac(nonemptymgi, element) > 0) {
      return MINPOP;
    }
    return 0.;
  }

  return nn;
}

// Set the partition function for an element in a cell based on the current level populations
void calculate_cellpartfuncts(const int nonemptymgi, const int element) {
  const int nions = get_nions(element);
  for (int ion = 0; ion < nions; ion++) {
    set_ion_partfunct(nonemptymgi, element, ion, calculate_partfunct(element, ion, nonemptymgi));
  }
}

// If not already set by the NLTE solver, set the ground level populations from either Saha LTE or
// ionisation/recombination balance (Photoionisation Equilibrium)
void set_groundlevelpops(const int nonemptymgi, const int element, const float nne, const bool force_saha) {
  const int nions = get_nions(element);

  if (nions <= 0) {
    return;
  }

  // calculate number density of the current element (derived from mass fraction)
  const double nnelement = grid::get_elem_numberdens(nonemptymgi, element);

  const bool use_phi_saha = force_saha || FORCE_SAHA_ION_BALANCE(get_atomicnumber(element));

  const auto ionfractions =
      (nnelement > 0) ? calculate_ionfractions(element, nonemptymgi, nne, use_phi_saha) : std::vector<double>();

  // -1 when the element is absent (ionfractions is empty); cast to int before subtracting to avoid
  // unsigned wraparound to a huge positive value.
  const int uppermost_ion = static_cast<int>(ionfractions.size()) - 1;
  const ptrdiff_t nincludedions = get_includedions();

  // Use ion fractions to calculate the groundlevel populations
  for (int ion = 0; ion < nions; ion++) {
    const int uniqueionindex = get_uniqueionindex(element, ion);
    double nnion{NAN};
    if (nnelement <= 0) {
      // absent element: every ion has zero population (do not floor to MINPOP)
      nnion = 0.;
    } else if (ion <= uppermost_ion) {
      nnion = std::max(MINPOP, nnelement * ionfractions[ion]);
    } else {
      nnion = MINPOP;
    }

    const auto groundpop =
        static_cast<float>(nnion * stat_weight(element, ion, 0) /
                           grid::ion_partfuncts_allcells[(nonemptymgi * nincludedions) + uniqueionindex]);
    assert_always(groundpop >= 0.);

    grid::ion_groundlevelpops_allcells[(nonemptymgi * nincludedions) + uniqueionindex] = groundpop;
  }
}

// Determine the electron number density for a given cell using a root
// solver and calculate the dependent level populations.
auto calculate_ion_balance_nne(const int nonemptymgi) -> void {
  const bool force_saha = globals::lte_iteration || grid::thick_allcells[nonemptymgi] == 1;

  const double nne_max = grid::get_rho(nonemptymgi) / MH;

  bool only_lowest_ionstage = true;  // could be completely neutral, or just at each element's lowest ion stage
  for (int element = 0; element < get_nelements(); element++) {
    if (grid::get_elem_massfrac(nonemptymgi, element) > 0) {
      const auto uppermost_ion = find_uppermost_ion(nonemptymgi, element, nne_max, force_saha);
      grid::set_elements_uppermost_ion(nonemptymgi, element, uppermost_ion);

      only_lowest_ionstage = only_lowest_ionstage && (uppermost_ion <= 0);
    } else {
      grid::set_elements_uppermost_ion(nonemptymgi, element, get_nions(element) - 1);
    }
  }

  if (only_lowest_ionstage) {
    set_groundlevelpops_neutral(nonemptymgi);
  } else {
    const auto nne_solution = find_converged_nne(nonemptymgi, nne_max, force_saha);
    grid::set_nne(nonemptymgi, nne_solution);

    for (int element = 0; element < get_nelements(); element++) {
      // avoid overwriting the ground level populations set by the NLTE pop solver
      const bool already_set_by_nlte_solver = !force_saha && elem_has_nlte_levels(element);
      if (!already_set_by_nlte_solver) {
        set_groundlevelpops(nonemptymgi, element, nne_solution, force_saha);
      }
    }
  }

  set_calculated_nne(nonemptymgi);
}
