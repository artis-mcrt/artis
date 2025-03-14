#include "ltepop.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "sn3d.h"

namespace {

struct nneSolutionParas {
  int nonemptymgi;
  bool force_saha;
};

auto interpolate_ions_spontrecombcoeff(const int uniqueionindex, const double T) -> double {
  const int lowerindex = floor(log(T / MINTEMP) / T_step_log);
  assert_testmodeonly(lowerindex >= 0);
  if (lowerindex < TABLESIZE - 1) {
    const int upperindex = lowerindex + 1;
    const double T_lower = MINTEMP * exp(lowerindex * T_step_log);
    const double T_upper = MINTEMP * exp(upperindex * T_step_log);

    const double f_upper = globals::ion_alpha_sp[(uniqueionindex * TABLESIZE) + upperindex];
    const double f_lower = globals::ion_alpha_sp[(uniqueionindex * TABLESIZE) + lowerindex];

    return f_lower + ((f_upper - f_lower) / (T_upper - T_lower) * (T - T_lower));
  }
  return globals::ion_alpha_sp[(uniqueionindex * TABLESIZE) + TABLESIZE - 1];
}

// use Saha equation for LTE ionization balance
auto phi_saha(const int element, const int ion, const int nonemptymgi) -> double {
  const int uniqueionindex = get_uniqueionindex(element, ion);
  const auto partfunc_ion =
      grid::ion_partfuncts_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) + uniqueionindex];
  const auto partfunc_upperion =
      grid::ion_partfuncts_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) + uniqueionindex + 1];

  const auto T_e = grid::get_Te(nonemptymgi);
  const double ionpot = epsilon(element, ion + 1, 0) - epsilon(element, ion, 0);
  const double partfunct_ratio = partfunc_ion / partfunc_upperion;
  return partfunct_ratio * SAHACONST * pow(T_e, -1.5) * exp(ionpot / KB / T_e);
}

// Calculate population ratio (a saha factor) of two consecutive ionisation stages in nebular approximation phi_j,k* =
// N_j,k*/(N_j+1,k* * nne)
auto phi_rate_balance(const int element, const int ion, const int nonemptymgi) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));

  assert_testmodeonly(!globals::lte_iteration);
  assert_testmodeonly(grid::modelgrid[nonemptymgi].thick != 1);  // should use use phi_lte instead

  assert_testmodeonly(!elem_has_nlte_levels(element));  // don't use this function if the NLTE solver is active

  const int uniqueionindex = get_uniqueionindex(element, ion);
  const auto partfunc_ion =
      grid::ion_partfuncts_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) + uniqueionindex];

  const auto T_e = grid::get_Te(nonemptymgi);

  // photoionisation plus collisional ionisation rate coefficient per ground level pop
  const double Gamma = globals::gammaestimator[get_ionestimindex_nonemptymgi(nonemptymgi, element, ion)];

  // Gamma is the photoionization rate per ground level pop
  const double Gamma_ion = Gamma * stat_weight(element, ion, 0) / partfunc_ion;

  const double Alpha_sp = interpolate_ions_spontrecombcoeff(uniqueionindex, T_e);

  // const double Col_rec = calculate_ionrecombcoeff(nonemptymgi, T_e, element, ion + 1, false, true, false, false,
  // false);
  const double Col_rec = 0.;

  const double gamma_nt = NT_ON ? nonthermal::nt_ionization_ratecoeff(nonemptymgi, element, ion) : 0.;

  if ((Gamma_ion + gamma_nt) == 0) {
    printout("Fatal: Gamma = 0 for element %d, ion %d in phi ... abort\n", element, ion);
    std::abort();
  }

  const double phi = (Alpha_sp + Col_rec) / (Gamma_ion + gamma_nt);

  // Y_nt should generally be higher than the Gamma term for nebular epoch

  if (!std::isfinite(phi) || phi == 0.) {
    const auto partfunc_upperion =
        grid::ion_partfuncts_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) + uniqueionindex + 1];
    printout(
        "[fatal] phi: phi %g exceeds numerically possible range for element %d, ion %d, T_e %g ... remove higher or "
        "lower ionisation stages\n",
        phi, element, ion, T_e);
    printout("[fatal] phi: Alpha_sp %g, Gamma %g, partfunct %g, stat_weight %g\n", Alpha_sp, Gamma_ion, partfunc_ion,
             stat_weight(element, ion, 0));
    printout("[fatal] phi: upperionpartfunct %g, upperionstatweight %g\n", partfunc_upperion,
             stat_weight(element, ion + 1, 0));
    printout("[fatal] phi: gamma_nt %g Col_rec %g grid::get_nne(nonemptymgi) %g\n", gamma_nt, Col_rec,
             grid::get_nne(nonemptymgi));
    std::abort();
  }

  return phi;
}

// calculate the free electron contribution from an element
auto get_element_nne_contrib(const int nonemptymgi, const int element) -> double {
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
auto nne_solution_f(const double nne_assumed, void *const voidparas) -> double {
  const auto *paras = static_cast<const nneSolutionParas *>(voidparas);
  const int nonemptymgi = paras->nonemptymgi;
  const bool force_saha = paras->force_saha;

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

auto calculate_levelpop_nominpop(const int nonemptymgi, const int element, const int ion, const int level,
                                 bool *const skipminpop) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));

  double nn{NAN};

  if (level == 0) {
    nn = get_groundlevelpop(nonemptymgi, element, ion);
  } else if (elem_has_nlte_levels(element)) {
    if (is_nlte(element, ion, level)) {
      // first_nlte refers to the first excited state (level=1)
      const double nltepop_over_rho = get_nlte_levelpop_over_rho(nonemptymgi, element, ion, level);
      if (nltepop_over_rho < -0.9) {
        // Case for when no NLTE level information is available yet
        nn = calculate_levelpop_boltzmann(nonemptymgi, element, ion, level);
      } else {
        nn = nltepop_over_rho * grid::get_rho(nonemptymgi);
        if (!std::isfinite(nn)) {
          printout("[fatal] NLTE population failure.\n");
          printout("element %d ion %d level %d\n", element, ion, level);
          printout("nn %g nltepop_over_rho %g rho %g\n", nn, nltepop_over_rho, grid::get_rho(nonemptymgi));
          printout("ground level %g\n", get_groundlevelpop(nonemptymgi, element, ion));
          std::abort();
        }
        *skipminpop = true;
        return nn;
      }
    } else {
      // level is in the superlevel
      assert_testmodeonly(level_isinsuperlevel(element, ion, level));

      const double superlevelpop_over_rho = get_nlte_superlevelpop_over_rho(nonemptymgi, element, ion);
      if (superlevelpop_over_rho < -0.9)  // TODO: should change this to less than zero?
      {
        // Case for when no NLTE level information is available yet
        nn = calculate_levelpop_boltzmann(nonemptymgi, element, ion, level);
      } else {
        nn = superlevelpop_over_rho * grid::get_rho(nonemptymgi) *
             superlevel_boltzmann(nonemptymgi, element, ion, level);
        if (!std::isfinite(nn)) {
          printout("[fatal] NLTE population failure.\n");
          printout("element %d ion %d level %d\n", element, ion, level);
          printout("nn %g superlevelpop_over_rho %g rho %g\n", nn, superlevelpop_over_rho, grid::get_rho(nonemptymgi));
          printout("ground level %g\n", get_groundlevelpop(nonemptymgi, element, ion));
          std::abort();
        }
        *skipminpop = true;
        return nn;
      }
    }
  } else {
    nn = calculate_levelpop_boltzmann(nonemptymgi, element, ion, level);
  }
  *skipminpop = false;
  return nn;
}

auto calculate_partfunct(const int element, const int ion, const int nonemptymgi) -> double
// Calculates the partition function for ion=ion of element=element in
// cell modelgridindex
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  double pop_store{NAN};

  const int uniqueionindex = get_uniqueionindex(element, ion);

  bool initial = false;
  if (get_groundlevelpop(nonemptymgi, element, ion) < MINPOP) {
    // either there really is none of this ion or this is a first pass through
    // in either case, we won't have any real nlte_populations so the actual value of
    // of groundlevelpop for this calculation doesn't matter, so long as it's not zero!
    pop_store = get_groundlevelpop(nonemptymgi, element, ion);
    initial = true;
    grid::ion_groundlevelpops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) + uniqueionindex] =
        1.;
  }

  double U = 1.;

  const int nlevels = get_nlevels(element, ion);
  const double groundpop = get_groundlevelpop(nonemptymgi, element, ion);
  for (int level = 1; level < nlevels; level++) {
    bool skipminpop = false;
    const double nn = calculate_levelpop_nominpop(nonemptymgi, element, ion, level, &skipminpop) / groundpop;
    U += nn;
  }
  U *= stat_weight(element, ion, 0);

  if (!std::isfinite(U)) {
    printout("element %d ion %d\n", element, ion);
    printout("modelgridindex %d\n", grid::get_mgi_of_nonemptymgi(nonemptymgi));
    printout("nlevels %d\n", nlevels);
    printout("sw %g\n", stat_weight(element, ion, 0));
    std::abort();
  }

  if (initial) {
    // put back the zero, just in case it matters for something
    grid::ion_groundlevelpops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) + uniqueionindex] =
        pop_store;
  }

  return U;
}

auto find_uppermost_ion(const int nonemptymgi, const int element, const double nne_hi, const bool force_saha) -> int {
  const int nions = get_nions(element);
  if (nions == 0) {
    return -1;
  }
  if (!force_saha && elem_has_nlte_levels(element)) {
    return nions - 1;
  }
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const bool use_phi_saha = force_saha || FORCE_SAHA_ION_BALANCE(get_atomicnumber(element));
  int uppermost_ion = 0;

  uppermost_ion = nions - 1;
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

  double factor = 1.;
  int ion = 0;
  for (ion = 0; ion < uppermost_ion; ion++) {
    const auto phifactor =
        use_phi_saha ? phi_saha(element, ion, nonemptymgi) : phi_rate_balance(element, ion, nonemptymgi);
    factor *= nne_hi * phifactor;

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

void set_calculated_nne(const int nonemptymgi) {
  double nne = 0.;  // free electron density
  for (int element = 0; element < get_nelements(); element++) {
    nne += get_element_nne_contrib(nonemptymgi, element);
  }

  grid::set_nne(nonemptymgi, std::max(MINPOP, nne));
}

// Special case of only neutral ions, set nne to some finite value so that packets are not lost in kpkts
void set_groundlevelpops_neutral(const int nonemptymgi) {
  printout("[warning] calculate_ion_balance_nne: only neutral ions in cell modelgridindex %d\n",
           grid::get_mgi_of_nonemptymgi(nonemptymgi));
  for (int element = 0; element < get_nelements(); element++) {
    const auto nnelement = grid::get_elem_numberdens(nonemptymgi, element);
    const int nions = get_nions(element);
    // Assign the species population to the neutral ion and set higher ions to MINPOP
    for (int ion = 0; ion < nions; ion++) {
      const int uniqueionindex = get_uniqueionindex(element, ion);
      double nnion{NAN};
      if (ion == 0) {
        nnion = nnelement;
      } else if (nnelement > 0.) {
        nnion = MINPOP;
      } else {
        nnion = 0.;
      }
      const double groundpop =
          (nnion * stat_weight(element, ion, 0) /
           grid::ion_partfuncts_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) + uniqueionindex]);

      grid::ion_groundlevelpops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) + uniqueionindex] =
          groundpop;
    }
  }
}

auto find_converged_nne(const int nonemptymgi, double nne_hi, const bool force_lte) -> float {
  // Search solution for nne in [nne_lo,nne_hi]

  nneSolutionParas paras = {.nonemptymgi = nonemptymgi, .force_saha = force_lte};
  gsl_function f = {.function = &nne_solution_f, .params = &paras};

  double nne_lo = 0.;  // MINPOP;
  if (nne_solution_f(nne_lo, f.params) * nne_solution_f(nne_hi, f.params) > 0) {
    const auto T_R = grid::get_TR(nonemptymgi);
    const auto T_e = grid::get_Te(nonemptymgi);
    const auto W = grid::get_W(nonemptymgi);
    const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    printout("n, nne_lo, nne_hi, T_R, T_e, W, rho %d, %g, %g, %g, %g, %g, %g\n", modelgridindex, nne_lo, nne_hi, T_R,
             T_e, W, grid::get_rho(nonemptymgi));
    printout("nne@x_lo %g\n", nne_solution_f(nne_lo, f.params));
    printout("nne@x_hi %g\n", nne_solution_f(nne_hi, f.params));

    for (int element = 0; element < get_nelements(); element++) {
      printout("modelgridindex %d, element %d, uppermost_ion is %d\n", modelgridindex, element,
               grid::get_elements_uppermost_ion(nonemptymgi, element));

      if constexpr (USE_LUT_PHOTOION) {
        for (int ion = 0; ion <= grid::get_elements_uppermost_ion(nonemptymgi, element); ion++) {
          printout("element %d, ion %d, gammaionest %g\n", element, ion,
                   globals::gammaestimator[get_ionestimindex_nonemptymgi(nonemptymgi, element, ion)]);
        }
      }
    }
  }

  double nne_solution = 0.;

  gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

  gsl_root_fsolver_set(solver, &f, nne_lo, nne_hi);
  constexpr int maxit = 50;
  constexpr double fractional_accuracy = 1e-3;
  int status = GSL_CONTINUE;
  int iter = 0;
  for (iter = 0; iter <= maxit; iter++) {
    gsl_root_fsolver_iterate(solver);
    nne_solution = gsl_root_fsolver_root(solver);
    nne_lo = gsl_root_fsolver_x_lower(solver);
    nne_hi = gsl_root_fsolver_x_upper(solver);
    status = gsl_root_test_interval(nne_lo, nne_hi, 0, fractional_accuracy);
    if (status != GSL_CONTINUE) {
      break;
    }
  }
  if (status == GSL_CONTINUE) {
    printout("[warning] calculate_ion_balance_nne: nne did not converge within %d iterations\n", iter + 1);
  }

  gsl_root_fsolver_free(solver);

  return std::max(MINPOP, nne_solution);
}

}  // anonymous namespace

// Calculate the fractions of an element's population in each ionization stage based on Saha LTE or ionisation
// equilibrium
[[nodiscard]] auto calculate_ionfractions(const int element, const int nonemptymgi, const double nne,
                                          const bool use_phi_saha) -> std::vector<double> {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
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
      printout("[warning] ionfract set to zero for ionstage %d of Z=%d in cell %d with T_e %g, T_R %g\n",
               get_ionstage(element, ion), get_atomicnumber(element), modelgridindex, grid::get_Te(nonemptymgi),
               grid::get_TR(nonemptymgi));
      ionfractions[ion] = 0;
    }
  }
  return ionfractions;
}

// Return the given ions groundlevel population for modelgridindex which was precalculated
// during update_grid and stored to the grid.
auto get_groundlevelpop(const int nonemptymgi, const int element, const int ion) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  const double nn = grid::ion_groundlevelpops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) +
                                                       get_uniqueionindex(element, ion)];
  if (nn < MINPOP) {
    if (grid::get_elem_abundance(nonemptymgi, element) > 0) {
      return MINPOP;
    }
    return 0.;
  }
  return nn;
}

// Calculate occupation population of a level assuming LTE excitation
auto calculate_levelpop_boltzmann(const int nonemptymgi, const int element, const int ion, const int level) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  const auto nnground = get_groundlevelpop(nonemptymgi, element, ion);
  if (level == 0) {
    return nnground;
  }

  const auto T_exc = LTEPOP_EXCITATION_USE_TJ ? grid::get_TJ(nonemptymgi) : grid::get_Te(nonemptymgi);

  const double E_aboveground = epsilon(element, ion, level) - epsilon(element, ion, 0);

  return (nnground * stat_weight(element, ion, level) / stat_weight(element, ion, 0) *
          exp(-E_aboveground / KB / T_exc));
}

auto calculate_levelpop(const int nonemptymgi, const int element, const int ion, const int level) -> double {
  bool skipminpop = false;
  double nn = calculate_levelpop_nominpop(nonemptymgi, element, ion, level, &skipminpop);
  if (!skipminpop && nn < MINPOP) {
    if (grid::get_elem_abundance(nonemptymgi, element) > 0) {
      nn = MINPOP;
    } else {
      nn = 0.;
    }
  }

  return nn;
}

// Calculate the population of a level from either LTE or NLTE information
__host__ __device__ auto get_levelpop(const int nonemptymgi, const int element, const int ion, const int level)
    -> double {
  double nn = 0.;
  if (use_cellcache) {
    assert_testmodeonly(globals::cellcache[cellcacheslotid].nonemptymgi == nonemptymgi);
    nn = globals::cellcache[cellcacheslotid].chelements[element].chions[ion].chlevels[level].population;
  } else {
    nn = calculate_levelpop(nonemptymgi, element, ion, level);
  }

  assert_testmodeonly(nn >= 0.);
  assert_testmodeonly(std::isfinite(nn));

  return nn;
}

// The partition functions depend only on T_R and W. This means they don't
// change during any iteration on T_e. Therefore their precalculation was
// taken out of calculate_ion_balance_nne to save runtime.
// TODO: not true if LTEPOP_EXCITATION_USE_TJ is true unless LTE mode only (TJ=TR=Te)
void calculate_cellpartfuncts(const int nonemptymgi, const int element) {
  const int nions = get_nions(element);
  for (int ion = 0; ion < nions; ion++) {
    grid::ion_partfuncts_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) +
                                  get_uniqueionindex(element, ion)] = calculate_partfunct(element, ion, nonemptymgi);
  }
}

// calculates saha factor in LTE: Phi_level,ion,element = nn_level,ion,element/(nne*nn_upper,ion+1,element)
__host__ __device__ auto calculate_sahafact(const int element, const int ion, const int level, const int upperionlevel,
                                            const double T, const double E_threshold) -> double {
  const double g_lower = stat_weight(element, ion, level);
  const double g_upper = stat_weight(element, ion + 1, upperionlevel);
  const double sf = SAHACONST * g_lower / g_upper * pow(T, -1.5) * exp(E_threshold / KB / T);
  if (sf < 0) {
    printout(
        "[fatal] calculate_sahafact: Negative Saha factor. sfac %g element %d ion %d level %d upperionlevel %d "
        "g_lower %g g_upper %g T %g E_threshold %g exppart %g\n",
        sf, element, ion, level, upperionlevel, g_lower, g_upper, T, E_threshold, exp(E_threshold / KB / T));
    std::abort();
  }
  return sf;
}

// Use the ground level population and partition function to get an ion population
[[nodiscard]] __host__ __device__ auto get_nnion(const int nonemptymgi, const int element, const int ion) -> double {
  const auto nnion = get_groundlevelpop(nonemptymgi, element, ion) *
                     grid::ion_partfuncts_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) +
                                                   get_uniqueionindex(element, ion)] /
                     stat_weight(element, ion, 0);
  assert_testmodeonly(nnion >= 0.);
  assert_testmodeonly(std::isfinite(nnion));
  return nnion;
}

// If not already set by the NLTE solver, set the ground level populations from either Saha LTE or
// ionization/recombination balance (Photoionization Equilibrium)
void set_groundlevelpops(const int nonemptymgi, const int element, const float nne, const bool force_saha) {
  const int nions = get_nions(element);

  if (nions <= 0) {
    return;
  }

  // calculate number density of the current element (abundances are given by mass)
  const double nnelement = grid::get_elem_numberdens(nonemptymgi, element);

  const bool use_phi_saha = force_saha || FORCE_SAHA_ION_BALANCE(get_atomicnumber(element));

  const auto ionfractions =
      (nnelement > 0) ? calculate_ionfractions(element, nonemptymgi, nne, use_phi_saha) : std::vector<double>();

  const int uppermost_ion = static_cast<int>(ionfractions.size() - 1);
  const ptrdiff_t nincludedions = get_includedions();

  // Use ion fractions to calculate the groundlevel populations
  for (int ion = 0; ion < nions; ion++) {
    const int uniqueionindex = get_uniqueionindex(element, ion);
    double nnion{NAN};
    if (ion <= uppermost_ion) {
      if (nnelement > 0) {
        nnion = std::max(MINPOP, nnelement * ionfractions[ion]);
      } else {
        nnion = 0.;
      }
    } else {
      nnion = MINPOP;
    }

    const double groundpop = nnion * stat_weight(element, ion, 0) /
                             grid::ion_partfuncts_allcells[(nonemptymgi * nincludedions) + uniqueionindex];

    if (!std::isfinite(groundpop)) {
      printout("[warning] calculate_ion_balance_nne: groundlevelpop infinite in connection with MINPOP\n");
    }

    grid::ion_groundlevelpops_allcells[(nonemptymgi * nincludedions) + uniqueionindex] = groundpop;
  }
}

// Determine the electron number density for a given cell using one of
// libgsl's root_solvers and calculates the depending level populations.
auto calculate_ion_balance_nne(const int nonemptymgi) -> void {
  const bool force_saha = globals::lte_iteration || grid::modelgrid[nonemptymgi].thick == 1;

  const double nne_hi = grid::get_rho(nonemptymgi) / MH;

  bool only_lowest_ionstage = true;  // could be completely neutral, or just at each element's lowest ion stage
  for (int element = 0; element < get_nelements(); element++) {
    if (grid::get_elem_abundance(nonemptymgi, element) > 0) {
      const int uppermost_ion = find_uppermost_ion(nonemptymgi, element, nne_hi, force_saha);
      grid::set_elements_uppermost_ion(nonemptymgi, element, uppermost_ion);

      only_lowest_ionstage = only_lowest_ionstage && (uppermost_ion <= 0);
    } else {
      grid::set_elements_uppermost_ion(nonemptymgi, element, get_nions(element) - 1);
    }
  }

  if (only_lowest_ionstage) {
    set_groundlevelpops_neutral(nonemptymgi);
  } else {
    const auto nne_solution = find_converged_nne(nonemptymgi, nne_hi, force_saha);
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
