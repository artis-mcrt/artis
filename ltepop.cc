#include "ltepop.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include <cmath>
#include <memory>

#include "atomic.h"
#include "decay.h"
#include "grid.h"
#include "macroatom.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "ratecoeff.h"
#include "sn3d.h"
#include "update_grid.h"

static auto interpolate_ions_spontrecombcoeff(const int element, const int ion, const double T) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_always(T >= MINTEMP);
  int const lowerindex = floor(log(T / MINTEMP) / T_step_log);
  if (lowerindex < TABLESIZE - 1) {
    int const upperindex = lowerindex + 1;
    double const T_lower = MINTEMP * exp(lowerindex * T_step_log);
    double const T_upper = MINTEMP * exp(upperindex * T_step_log);

    double const f_upper = globals::elements[element].ions[ion].Alpha_sp[upperindex];
    double const f_lower = globals::elements[element].ions[ion].Alpha_sp[lowerindex];

    return f_lower + (f_upper - f_lower) / (T_upper - T_lower) * (T - T_lower);
  }
  return globals::elements[element].ions[ion].Alpha_sp[TABLESIZE - 1];
}

static auto phi(const int element, const int ion, const int modelgridindex) -> double
/// Calculates population ratio (a saha factor) of two consecutive ionisation stages
/// in nebular approximation phi_j,k* = N_j,k*/(N_j+1,k* * nne)
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));

  double phi = 0;
  auto partfunc_ion = grid::modelgrid[modelgridindex].composition[element].partfunct[ion];
  auto partfunc_upperion = grid::modelgrid[modelgridindex].composition[element].partfunct[ion + 1];

  const auto T_e = grid::get_Te(modelgridindex);
  /// Newest ionisation formula

  const bool use_lte_ratio = (globals::lte_iteration || grid::modelgrid[modelgridindex].thick == 1);

  if (use_lte_ratio) {
    const double ionpot = epsilon(element, ion + 1, 0) - epsilon(element, ion, 0);
    // printout("ionpot for element %d, ion %d is %g\n", element, ion, ionpot / EV);
    const double partfunct_ratio = partfunc_ion / partfunc_upperion;
    return partfunct_ratio * SAHACONST * pow(T_e, -1.5) * exp(ionpot / KB / T_e);
  }

  double Gamma = 0.;
  if constexpr (!USE_LUT_PHOTOION) {
    Gamma = calculate_iongamma_per_gspop(modelgridindex, element, ion);
  } else {
    Gamma = globals::gammaestimator[get_ionestimindex(modelgridindex, element, ion)];
  }

  // Gamma is the photoionization rate per ground level pop
  const double Gamma_ion = Gamma * stat_weight(element, ion, 0) / partfunc_ion;

  if (Gamma == 0. && (!NT_ON || (globals::rpkt_emiss[modelgridindex] == 0. &&
                                 grid::get_modelinitradioabund(modelgridindex, decay::get_nucindex(24, 48)) == 0. &&
                                 grid::get_modelinitradioabund(modelgridindex, decay::get_nucindex(28, 56)) == 0.))) {
    printout("Fatal: Gamma = 0 for element %d, ion %d in phi ... abort\n", element, ion);
    abort();
  }

  double Alpha_sp = 0.;
  if constexpr (NLTE_POPS_ON) {
    Alpha_sp =
        calculate_ionrecombcoeff(modelgridindex, T_e, element, ion + 1, false, false, false, false, false, false);
  } else {
    Alpha_sp = interpolate_ions_spontrecombcoeff(element, ion, T_e);
  }

  // const double Col_rec = calculate_ionrecombcoeff(modelgridindex, T_e, element, ion + 1, false, true, false, false,
  // false);
  const double Col_rec = 0.;

  const double gamma_nt = NT_ON ? nonthermal::nt_ionization_ratecoeff(modelgridindex, element, ion) : 0.;

  phi = (Alpha_sp + Col_rec) / (Gamma_ion + gamma_nt);

  // Y_nt should generally be higher than the Gamma term for nebular epoch

  if (!std::isfinite(phi) || phi == 0.) {
    printout(
        "[fatal] phi: phi %g exceeds numerically possible range for element %d, ion %d, T_e %g ... remove higher or "
        "lower ionisation stages\n",
        phi, element, ion, T_e);
    printout("[fatal] phi: Alpha_sp %g, Gamma %g, partfunct %g, stat_weight %g\n", Alpha_sp, Gamma, partfunc_ion,
             stat_weight(element, ion, 0));
    printout("[fatal] phi: upperionpartfunct %g, upperionstatweight %g\n", partfunc_upperion,
             stat_weight(element, ion + 1, 0));
    printout("[fatal] phi: gamma_nt %g Col_rec %g grid::get_nne(modelgridindex) %g\n", gamma_nt, Col_rec,
             grid::get_nne(modelgridindex));
    abort();
  }

  return phi;
}

static auto get_ionfractions(const int element, const int modelgridindex, const double nne) -> std::vector<double>
// Calculate the fractions of an element's population in each ionization stage
{
  const int uppermost_ion = grid::get_elements_uppermost_ion(modelgridindex, element);
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(uppermost_ion <= std::max(0, get_nions(element) - 1));

  std::vector<double> ionfractions(uppermost_ion + 1);
  std::vector<double> nnionfactor(uppermost_ion + 1);
  nnionfactor[uppermost_ion] = 1;

  double denominator = 1.;

  for (int ion = uppermost_ion - 1; ion >= 0; ion--) {
    nnionfactor[ion] = nnionfactor[ion + 1] * nne * phi(element, ion, modelgridindex);
    denominator += nnionfactor[ion];
  }

  for (int ion = 0; ion <= uppermost_ion; ion++) {
    ionfractions[ion] = nnionfactor[ion] / denominator;

    if (!std::isfinite(ionfractions[ion])) {
      if (modelgridindex != grid::get_npts_model()) {
        printout("[warning] ionfract set to zero for ionstage %d of Z=%d in cell %d with T_e %g, T_R %g\n",
                 get_ionstage(element, ion), get_atomicnumber(element), modelgridindex, grid::get_Te(modelgridindex),
                 grid::get_TR(modelgridindex));
        // abort();
        ionfractions[ion] = 0;
      }
    }
    // printout("ionfract(%d,%d,%d,%g) = %g\n", element, ion, modelgridindex, nne, ionfractions[ion]);
  }
  return ionfractions;
}

static auto nne_solution_f(double nne, void *paras) -> double
// assume a value for nne and then calculate the resulting nne
// the difference between the assumed and calculated nne is returned
{
  const int modelgridindex = (static_cast<struct nne_solution_paras *>(paras))->cellnumber;

  double outersum = 0.;
  for (int element = 0; element < get_nelements(); element++) {
    const double massfrac = grid::modelgrid[modelgridindex].composition[element].abundance;
    if (massfrac > 0 && get_nions(element) > 0) {
      const auto ionfractions = get_ionfractions(element, modelgridindex, nne);
      const int uppermost_ion = static_cast<int>(ionfractions.size() - 1);
      double elem_nne_contrib = 0.;
      for (int ion = 0; ion <= uppermost_ion; ion++) {
        elem_nne_contrib += (get_ionstage(element, ion) - 1) * ionfractions[ion];
      }

      const double elem_meanweight = grid::get_element_meanweight(modelgridindex, element);
      outersum += massfrac / elem_meanweight * elem_nne_contrib;

      if (!std::isfinite(outersum)) {
        printout("nne_solution_f: element %d uppermostion %d massfrac %g, mass %g\n", element,
                 grid::get_elements_uppermost_ion(modelgridindex, element), massfrac, elem_meanweight);
        printout("outersum %g\n", outersum);
        abort();
      }
    }
  }

  const double rho = grid::get_rho(modelgridindex);
  return rho * outersum - nne;
}

auto get_groundlevelpop(int modelgridindex, int element, int ion) -> double
/// Returns the given ions groundlevel population for modelgridindex which was precalculated
/// during update_grid and stored to the grid.
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));

  const double nn = grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion];
  if (nn < MINPOP) {
    if (grid::get_elem_abundance(modelgridindex, element) > 0) {
      return MINPOP;
    }
    return 0.;
  }
  return nn;
}

auto calculate_levelpop_lte(int modelgridindex, int element, int ion, int level) -> double
/// Calculates occupation population of a level assuming LTE excitation
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));

  const auto groundlevelpop = get_groundlevelpop(modelgridindex, element, ion);
  if (level == 0) {
    return groundlevelpop;
  }

  const double T_exc = LTEPOP_EXCITATION_USE_TJ ? grid::get_TJ(modelgridindex) : grid::get_Te(modelgridindex);

  const double E_aboveground = epsilon(element, ion, level) - epsilon(element, ion, 0);
  const double nnground = groundlevelpop;

  return (nnground * stat_weight(element, ion, level) / stat_weight(element, ion, 0) *
          exp(-E_aboveground / KB / T_exc));
}

static auto calculate_levelpop_nominpop(int modelgridindex, int element, int ion, int level, bool *skipminpop)
    -> double {
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));

  double nn = NAN;

  if (level == 0) {
    nn = get_groundlevelpop(modelgridindex, element, ion);
  } else if (elem_has_nlte_levels(element)) {
    if (is_nlte(element, ion, level)) {
      const double nltepop_over_rho =
          grid::modelgrid[modelgridindex].nlte_pops[globals::elements[element].ions[ion].first_nlte + level - 1];
      if (nltepop_over_rho < -0.9) {
        // Case for when no NLTE level information is available yet
        nn = calculate_levelpop_lte(modelgridindex, element, ion, level);
      } else {
        // printout("Using an nlte population!\n");
        nn = nltepop_over_rho * grid::get_rho(modelgridindex);
        if (!std::isfinite(nn)) {
          printout("[fatal] NLTE population failure.\n");
          printout("element %d ion %d level %d\n", element, ion, level);
          printout("nn %g nltepop_over_rho %g rho %g\n", nn, nltepop_over_rho, grid::get_rho(modelgridindex));
          printout("ground level %g\n", get_groundlevelpop(modelgridindex, element, ion));
          abort();
        }
        *skipminpop = true;
        return nn;
      }
    } else {
      // level is in the superlevel
      assert_testmodeonly(level_isinsuperlevel(element, ion, level));

      const int sl_nlte_index = globals::elements[element].ions[ion].first_nlte + get_nlevels_nlte(element, ion);

      const double superlevelpop_over_rho = grid::modelgrid[modelgridindex].nlte_pops[sl_nlte_index];
      if (superlevelpop_over_rho < -0.9)  // TODO: should change this to less than zero?
      {
        // Case for when no NLTE level information is available yet
        nn = calculate_levelpop_lte(modelgridindex, element, ion, level);
      } else {
        // printout("Using a superlevel population!\n");
        nn = superlevelpop_over_rho * grid::get_rho(modelgridindex) *
             superlevel_boltzmann(modelgridindex, element, ion, level);
        if (!std::isfinite(nn)) {
          printout("[fatal] NLTE population failure.\n");
          printout("element %d ion %d level %d\n", element, ion, level);
          printout("nn %g superlevelpop_over_rho %g rho %g\n", nn, superlevelpop_over_rho,
                   grid::get_rho(modelgridindex));
          printout("ground level %g\n", get_groundlevelpop(modelgridindex, element, ion));
          abort();
        }
        *skipminpop = true;
        return nn;
      }
    }
  } else {
    nn = calculate_levelpop_lte(modelgridindex, element, ion, level);
  }
  *skipminpop = false;
  return nn;
}

auto calculate_levelpop(int modelgridindex, int element, int ion, int level) -> double {
  bool skipminpop = false;
  double nn = calculate_levelpop_nominpop(modelgridindex, element, ion, level, &skipminpop);
  if (!skipminpop && nn < MINPOP) {
    if (grid::get_elem_abundance(modelgridindex, element) > 0) {
      nn = MINPOP;
      // nn = calculate_levelpop_lte(modelgridindex, element, ion, level);
    } else {
      nn = 0.;
    }
  }

  return nn;
}

auto get_levelpop(int modelgridindex, int element, int ion, int level) -> double
/// Calculates the population of a level from either LTE or NLTE information
{
  double nn = 0.;
  if (use_cellhist) {
    assert_testmodeonly(modelgridindex == globals::cellhistory[tid].cellnumber);
    nn = globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].population;
  } else {
    nn = calculate_levelpop(modelgridindex, element, ion, level);
  }

  assert_testmodeonly(nn >= 0.);
  assert_testmodeonly(std::isfinite(nn));

  return nn;
}

static auto calculate_partfunct(int element, int ion, int modelgridindex) -> double
/// Calculates the partition function for ion=ion of element=element in
/// cell modelgridindex
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  double pop_store = NAN;
  // double E_level, E_ground, test;

  int initial = 0;
  if (get_groundlevelpop(modelgridindex, element, ion) < MINPOP) {
    // either there reall is none of this ion or this is a first pass through
    // in either case, we won't have any real nlte_populations so the actual value of
    // of groundlevelpop for this calculation doesn't matter, so long as it's not zero!
    pop_store = get_groundlevelpop(modelgridindex, element, ion);
    initial = 1;
    grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = 1.0;
  }

  // printout("groundlevelpop %g\n", get_groundlevelpop(modelgridindex,element,ion));

  double U = 1.;

  const int nlevels = get_nlevels(element, ion);
  const double groundpop = get_groundlevelpop(modelgridindex, element, ion);
  for (int level = 1; level < nlevels; level++) {
    bool skipminpop = false;
    const double nn = calculate_levelpop_nominpop(modelgridindex, element, ion, level, &skipminpop) / groundpop;
    // const double nn = get_levelpop(modelgridindex, element, ion, level) / groundpop;
    U += nn;
  }
  U *= stat_weight(element, ion, 0);

  if (!std::isfinite(U)) {
    printout("element %d ion %d\n", element, ion);
    printout("modelgridindex %d\n", modelgridindex);
    printout("nlevels %d\n", nlevels);
    printout("sw %g\n", stat_weight(element, ion, 0));
    abort();
  }

  if (initial == 1) {
    // put back the zero, just in case it matters for something
    grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = pop_store;
  }

  return U;
}

void calculate_cellpartfuncts(const int modelgridindex, const int element)
/// The partition functions depend only on T_R and W. This means they don't
/// change during any iteration on T_e. Therefore their precalculation was
/// taken out of calculate_ion_balance_nne to save runtime.
// TODO: not true if LTEPOP_EXCITATION_USE_TJ is true unless LTE mode only (TJ=TR=Te)
{
  /// Precalculate partition functions for each ion in every cell
  /// this saves a factor 10 in calculation time of Saha-Boltzman populations
  const int nions = get_nions(element);
  for (int ion = 0; ion < nions; ion++) {
    grid::modelgrid[modelgridindex].composition[element].partfunct[ion] =
        calculate_partfunct(element, ion, modelgridindex);
  }
}

auto calculate_sahafact(int element, int ion, int level, int upperionlevel, double T, double E_threshold) -> double
/// calculates saha factor in LTE: Phi_level,ion,element = nn_level,ion,element/(nne*nn_upper,ion+1,element)
{
  const double g_lower = stat_weight(element, ion, level);
  const double g_upper = stat_weight(element, ion + 1, upperionlevel);
  const double sf = SAHACONST * g_lower / g_upper * pow(T, -1.5) * exp(E_threshold / KB / T);
  // printout("element %d, ion %d, level %d, T, %g, E %g has sf %g (g_l %g g_u %g)\n", element, ion, level, T,
  // E_threshold, sf,stat_weight(element,ion,level),stat_weight(element,ion+1,0) );
  if (sf < 0) {
    printout(
        "[fatal] calculate_sahafact: Negative Saha factor. sfac %g element %d ion %d level %d upperionlevel %d "
        "g_lower %g g_upper %g T %g E_threshold %g exppart %g\n",
        sf, element, ion, level, upperionlevel, g_lower, g_upper, T, E_threshold, exp(E_threshold / KB / T));
    abort();
  }
  return sf;
}

auto ionstagepop(int modelgridindex, int element, int ion) -> double
/// Calculates the given ionstage's population density using the ground level population and partition function
{
  return get_groundlevelpop(modelgridindex, element, ion) *
         grid::modelgrid[modelgridindex].composition[element].partfunct[ion] / stat_weight(element, ion, 0);
}

static auto find_uppermost_ion(const int modelgridindex, const int element, const double nne_hi) -> int {
  int uppermost_ion = 0;
  const int nions = get_nions(element);
  if (globals::lte_iteration || grid::modelgrid[modelgridindex].thick == 1) {
    uppermost_ion = nions - 1;
  } else {
    int ion = -1;
    for (ion = 0; ion < nions - 1; ion++) {
      double Gamma = 0.;
      if constexpr (!USE_LUT_PHOTOION) {
        Gamma = calculate_iongamma_per_gspop(modelgridindex, element, ion);
      } else {
        Gamma = globals::gammaestimator[get_ionestimindex(modelgridindex, element, ion)];
      }

      if ((Gamma == 0) &&
          (!NT_ON || ((globals::rpkt_emiss[modelgridindex] == 0.) &&
                      (grid::get_modelinitradioabund(modelgridindex, decay::get_nucindex(24, 48)) == 0.) &&
                      (grid::get_modelinitradioabund(modelgridindex, decay::get_nucindex(28, 56)) == 0.)))) {
        break;
      }
    }
    uppermost_ion = ion;
  }

  double factor = 1.;
  int ion = 0;
  for (ion = 0; ion < uppermost_ion; ion++) {
    factor *= nne_hi * phi(element, ion, modelgridindex);

    if (!std::isfinite(factor)) {
      printout(
          "[info] calculate_ion_balance_nne: uppermost_ion limited by phi factors for element "
          "Z=%d, ionstage %d in cell %d\n",
          get_atomicnumber(element), get_ionstage(element, ion), modelgridindex);
      return ion;
    }
  }
  uppermost_ion = ion;
  return uppermost_ion;
}

static auto calculate_nntot_nne(const int modelgridindex) -> double {
  double nne = 0.;  // free electron density
  double nntot = 0.;

  for (int element = 0; element < get_nelements(); element++) {
    // calculate number density of the current element (abundances are given by mass)
    const double nnelement = grid::get_elem_numberdens(modelgridindex, element);

    // Use ionization fractions to calculate the free electron contributions
    if (nnelement > 0) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        const auto nnion = ionstagepop(modelgridindex, element, ion);
        const int ioncharge = get_ionstage(element, ion) - 1;
        nne += ioncharge * nnion;
        nntot += nnion;
      }
    }
  }
  nntot += nne;

  nne = std::max(MINPOP, nne);
  grid::set_nne(modelgridindex, nne);
  return nntot;
}

auto calculate_ion_balance_nne(const int modelgridindex) -> double
/// Determines the electron number density for a given cell using one of
/// libgsl's root_solvers and calculates the depending level populations.
{
  const bool allow_nlte =
      (globals::total_nlte_levels > 0) && !globals::lte_iteration && (grid::modelgrid[modelgridindex].thick != 1);

  if (allow_nlte) {
    return calculate_nntot_nne(modelgridindex);
  }

  double nne_hi = grid::get_rho(modelgridindex) / MH;

  bool only_lowest_ionstage = true;  // could be completely neutral, or just at each element's lowest ion stage
  for (int element = 0; element < get_nelements(); element++) {
    if (grid::get_elem_abundance(modelgridindex, element) > 0) {
      const int uppermost_ion =
          allow_nlte ? get_nions(element) - 1 : find_uppermost_ion(modelgridindex, element, nne_hi);
      grid::set_elements_uppermost_ion(modelgridindex, element, uppermost_ion);

      only_lowest_ionstage = only_lowest_ionstage && (uppermost_ion <= 0);
    } else {
      grid::set_elements_uppermost_ion(modelgridindex, element, get_nions(element) - 1);
    }
  }

  if (only_lowest_ionstage) {
    /// Special case of only neutral ions, set nne to some finite value that
    /// packets are not lost in kpkts
    /// Introduce a flag variable which is sent to the T_e solver so that
    /// we get this info only once when T_e is converged and not for each
    /// iteration step.
    printout("[warning] calculate_ion_balance_nne: only neutral ions in cell modelgridindex %d\n", modelgridindex);
    /// Now calculate the ground level populations in nebular approximation and store them to the
    /// grid
    double nne = 0.;
    double nntot = 0.;
    for (int element = 0; element < get_nelements(); element++) {
      /// calculate number density of the current element (abundances are given by mass)
      const auto nnelement = grid::get_elem_numberdens(modelgridindex, element);
      const int nions = get_nions(element);
      /// Assign the species population to the neutral ion and set higher ions to MINPOP
      for (int ion = 0; ion < nions; ion++) {
        double nnion = NAN;
        if (ion == 0) {
          nnion = nnelement;
        } else if (nnelement > 0.) {
          nnion = MINPOP;
        } else {
          nnion = 0.;
        }
        nntot += nnion;
        nne += nnion * (get_ionstage(element, ion) - 1);
        grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion] =
            (nnion * stat_weight(element, ion, 0) /
             grid::modelgrid[modelgridindex].composition[element].partfunct[ion]);

        if (!std::isfinite(grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion])) {
          printout(
              "[warning] calculate_ion_balance_nne: groundlevelpop infinite in connection with "
              "MINPOP\n");
        }
      }
    }
    nntot += nne;
    nne = std::max(MINPOP, nne);
    grid::set_nne(modelgridindex, nne);
    return nntot;
  }

  /// Search solution for nne in [nne_lo,nne_hi]

  struct nne_solution_paras paras = {.cellnumber = modelgridindex};
  gsl_function f = {.function = &nne_solution_f, .params = &paras};

  double nne_lo = 0.;  // MINPOP;
  if (nne_solution_f(nne_lo, f.params) * nne_solution_f(nne_hi, f.params) > 0) {
    const auto T_R = grid::get_TR(modelgridindex);
    const auto T_e = grid::get_Te(modelgridindex);
    const auto W = grid::get_W(modelgridindex);
    printout("n, nne_lo, nne_hi, T_R, T_e, W, rho %d, %g, %g, %g, %g, %g, %g\n", modelgridindex, nne_lo, nne_hi, T_R,
             T_e, W, grid::get_rho(modelgridindex));
    printout("nne@x_lo %g\n", nne_solution_f(nne_lo, f.params));
    printout("nne@x_hi %g\n", nne_solution_f(nne_hi, f.params));

    for (int element = 0; element < get_nelements(); element++) {
      printout("cell %d, element %d, uppermost_ion is %d\n", modelgridindex, element,
               grid::get_elements_uppermost_ion(modelgridindex, element));

      if constexpr (USE_LUT_PHOTOION) {
        for (int ion = 0; ion <= grid::get_elements_uppermost_ion(modelgridindex, element); ion++) {
          printout("element %d, ion %d, gammaionest %g\n", element, ion,
                   globals::gammaestimator[get_ionestimindex(modelgridindex, element, ion)]);
        }
      }
    }
  }

  double nne = 0.;

  gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

  gsl_root_fsolver_set(solver, &f, nne_lo, nne_hi);
  constexpr int maxit = 50;
  constexpr double fractional_accuracy = 1e-3;
  int status = GSL_CONTINUE;
  for (int iter = 0; iter <= maxit; iter++) {
    gsl_root_fsolver_iterate(solver);
    nne = gsl_root_fsolver_root(solver);
    nne_lo = gsl_root_fsolver_x_lower(solver);
    nne_hi = gsl_root_fsolver_x_upper(solver);
    status = gsl_root_test_interval(nne_lo, nne_hi, 0, fractional_accuracy);
    if (status != GSL_CONTINUE) {
      break;
    }
  }
  if (status == GSL_CONTINUE) {
    printout("[warning] calculate_ion_balance_nne: nne did not converge within %d iterations\n", maxit);
  }

  gsl_root_fsolver_free(solver);

  nne = std::max(MINPOP, nne);

  grid::set_nne(modelgridindex, nne);

  /// Now calculate the ground level populations in nebular approximation and store them to the
  /// grid
  double nntot = nne;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    /// calculate number density of the current element (abundances are given by mass)
    const double nnelement = grid::get_elem_numberdens(modelgridindex, element);

    const auto ionfractions = (nnelement > 0) ? get_ionfractions(element, modelgridindex, nne) : std::vector<double>();

    const int uppermost_ion = static_cast<int>(ionfractions.size() - 1);

    /// Use ionizationfractions to calculate the groundlevel populations
    for (int ion = 0; ion < nions; ion++) {
      double nnion = NAN;
      if (ion <= uppermost_ion) {
        if (nnelement > 0) {
          nnion = std::max(MINPOP, nnelement * ionfractions[ion]);
        } else {
          nnion = 0.;
        }
      } else {
        nnion = MINPOP;
      }
      nntot += nnion;

      const double groundpop =
          (nnion * stat_weight(element, ion, 0) / grid::modelgrid[modelgridindex].composition[element].partfunct[ion]);

      if (!std::isfinite(groundpop)) {
        printout("[warning] calculate_ion_balance_nne: groundlevelpop infinite in connection with MINPOP\n");
      }

      grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = groundpop;
    }
  }

  return nntot;
}