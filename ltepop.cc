#include "ltepop.h"

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

auto nne_solution_f(double x, void *paras) -> double
/// For libgsl bracketing type solver
/// provides the equation which has to be solved to obtain the electron number
/// density (passed by x)
{
  const int modelgridindex = (static_cast<struct nne_solution_paras *>(paras))->cellnumber;
  const double rho = grid::get_rho(modelgridindex);

  double outersum = 0.;
  // printout("debug get_nelements() %d =========================\n",get_nelements());
  for (int element = 0; element < get_nelements(); element++) {
    const double abundance = grid::modelgrid[modelgridindex].composition[element].abundance;
    if (abundance > 0 && get_nions(element) > 0) {
      const double elem_meanweight = grid::get_element_meanweight(modelgridindex, element);
      double innersum = 0.;
      const int uppermost_ion = grid::get_elements_uppermost_ion(modelgridindex, element);

      auto ionfractions = get_ionfractions(element, modelgridindex, x, uppermost_ion);

      int ion = 0;
      for (ion = 0; ion <= uppermost_ion; ion++) {
        // printout("debug element %d, ion %d, ionfract(element,ion,T,x) %g\n",element,ion,ionfractions[ion]);
        innersum += (get_ionstage(element, ion) - 1) * ionfractions[ion];
      }
      assert_always(std::isfinite(innersum));
      outersum += abundance / elem_meanweight * innersum;
      if (!std::isfinite(outersum)) {
        printout("nne_solution_f: element %d ion %d uppermostion %d abundance %g, mass %g\n", element, ion,
                 grid::get_elements_uppermost_ion(modelgridindex, element), abundance, elem_meanweight);
        printout("outersum %g\n", outersum);
        abort();
      }
    }
  }

  return rho * outersum - x;
}

std::vector<double> get_ionfractions(const int element, const int modelgridindex, const double nne,
                                     const int uppermost_ion)
// Calculate the fractions of an element's population in each ionization stage
{
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

auto phi(const int element, const int ion, const int modelgridindex) -> double
/// Calculates population ratio (a saha factor) of two consecutive ionisation stages
/// in nebular approximation phi_j,k* = N_j,k*/(N_j+1,k* * nne)
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));

  double phi = 0;

  // double Y_nt, ionpot_in;
  // int element_in, ion_in, nions_in;
  // double rate_use;

  const auto T_e = grid::get_Te(modelgridindex);
  // double T_R = grid::get_TR(modelgridindex);

  // double W = globals::cell[cellnumber].W;

  /// Old ionisation formula
  // partfunct_ratio =
  // globals::cell[cellnumber].composition[element].partfunct[ion]/globals::cell[cellnumber].composition[element].partfunct[ion+1];
  // phi = 1./W * sqrt(T_R/T_e) * partfunct_ratio * SAHACONST * pow(T_R,-1.5) * exp(ionpot/KB/T_R);

  /// New ionisation formula with zeta
  // zeta = interpolate_zeta(element,ion,T_e);
  // phi = 1./W * 1./(zeta+W*(1-zeta)) * sqrt(T_R/T_e) * partfunct_ratio * SAHACONST * pow(T_R,-1.5) *
  // exp(ionpot/KB/T_R);

  /// Newest ionisation formula

  const bool use_lte_ratio = (globals::initial_iteration || grid::modelgrid[modelgridindex].thick == 1);

  if (use_lte_ratio) {
    const double ionpot = epsilon(element, ion + 1, 0) - epsilon(element, ion, 0);
    // printout("ionpot for element %d, ion %d is %g\n", element, ion, ionpot / EV);
    const double partfunct_ratio = grid::modelgrid[modelgridindex].composition[element].partfunct[ion] /
                                   grid::modelgrid[modelgridindex].composition[element].partfunct[ion + 1];
    phi = partfunct_ratio * SAHACONST * pow(T_e, -1.5) * exp(ionpot / KB / T_e);
  } else
  // elseif (NLTE_POPS_ALL_IONS_SIMULTANEOUS)
  //     {
  //       const float nne = grid::get_nne(modelgridindex);
  //       phi = ionstagepop(modelgridindex,element,ion) / ionstagepop(modelgridindex,element,ion+1) / nne;
  //     }
  // else
  {
    double Gamma = 0.;
    if constexpr (!USE_LUT_PHOTOION) {
      Gamma = calculate_iongamma_per_gspop(modelgridindex, element, ion);
    } else {
      Gamma =
          globals::gammaestimator[modelgridindex * get_nelements() * get_max_nions() + element * get_max_nions() + ion];
    }
    // printout("phicompare element %d ion %d T_e = %g gammaestimator %g calculate_iongamma_per_gspop %g\n",
    //          element, ion, T_e,
    //          globals::gammaestimator[modelgridindex * get_nelements() * get_max_nions() + element * get_max_nions() +
    //          ion], calculate_iongamma_per_gspop(modelgridindex, element, ion));

    // Gamma is the photoionization rate per ground level pop
    const double Gamma_ion =
        Gamma * stat_weight(element, ion, 0) / grid::modelgrid[modelgridindex].composition[element].partfunct[ion];

    if (Gamma == 0. && (!NT_ON || (globals::rpkt_emiss[modelgridindex] == 0. &&
                                   grid::get_modelinitradioabund(modelgridindex, decay::get_nucindex(24, 48)) == 0. &&
                                   grid::get_modelinitradioabund(modelgridindex, decay::get_nucindex(28, 56)) == 0.))) {
      printout("Fatal: Gamma = 0 for element %d, ion %d in phi ... abort\n", element, ion);
      abort();
    }

    // Alpha_st = stimrecombestimator[cellnumber*get_nelements()*get_max_nions()+element*get_max_nions()+ion];
    double const Alpha_st = 0.;  /// approximate treatment neglects stimulated recombination

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

    double Y_nt = 0.0;

    if constexpr (NT_ON) {
      Y_nt = nonthermal::nt_ionization_ratecoeff(modelgridindex, element, ion);
    }

    // || !std::isfinite(Gamma))
    // return phi_lte(element,ion,cellnumber);
    // gamma_lte = interpolate_photoioncoeff_below(element,ion,0,T_e) +
    // interpolate_photoioncoeff_above(element,ion,0,T_e); zeta = interpolate_zeta(element,ion,T_e); alpha_sp =
    // interpolate_spontrecombcoeff(element,ion,0,T_e); phi = gamma_lte*(Alpha_sp+Apha_st)/(Gamma*alpha_sp) *
    // partfunct_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e);

    phi = (Alpha_sp + Alpha_st + Col_rec) / (Gamma_ion + Y_nt);

    // Y_nt should generally be higher than the Gamma term for nebular epoch

    if (!std::isfinite(phi) || phi == 0.) {
      printout(
          "[fatal] phi: phi %g exceeds numerically possible range for element %d, ion %d, T_e %g ... remove higher or "
          "lower ionisation stages\n",
          phi, element, ion, T_e);
      printout("[fatal] phi: Alpha_sp %g, Alpha_st %g, Gamma %g, partfunct %g, stat_weight %g\n", Alpha_sp, Alpha_st,
               Gamma, grid::modelgrid[modelgridindex].composition[element].partfunct[ion],
               stat_weight(element, ion, 0));
      printout("[fatal] phi: upperionpartfunct %g, upperionstatweight %g\n",
               grid::modelgrid[modelgridindex].composition[element].partfunct[ion + 1],
               stat_weight(element, ion + 1, 0));
      printout("[fatal] phi: Y_nt %g Col_rec %g grid::get_nne(modelgridindex) %g\n", Y_nt, Col_rec,
               grid::get_nne(modelgridindex));
      // abort();
    }
  }

  return phi;
}

auto get_groundlevelpop(int modelgridindex, int element, int ion) -> double
/// Returns the given ions groundlevel population for modelgridindex which was precalculated
/// during update_grid and stored to the grid.
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  // double nn = grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion];
  // if (nn < MINPOP) nn = MINPOP;
  // return nn;

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
  if (level == 0) {
    return get_groundlevelpop(modelgridindex, element, ion);
  }

  const double T_exc = LTEPOP_EXCITATION_USE_TJ ? grid::get_TJ(modelgridindex) : grid::get_Te(modelgridindex);
  const double W = 1.;

  const double E_level = epsilon(element, ion, level);
  const double E_ground = epsilon(element, ion, 0);
  const double nnground = get_groundlevelpop(modelgridindex, element, ion);

  return (nnground * W * stat_weight(element, ion, level) / stat_weight(element, ion, 0) *
          exp(-(E_level - E_ground) / KB / T_exc));
}

static auto calculate_levelpop_nominpop(int modelgridindex, int element, int ion, int level, bool *skipminpop)
    -> double {
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));

  double nn = NAN;

  //  T_exc = MINTEMP;

  if (level == 0) {
    nn = get_groundlevelpop(modelgridindex, element, ion);
  } else if constexpr (NLTE_POPS_ON) {
    if (is_nlte(element, ion, level)) {
      // printout("Using an nlte population!\n");
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
    } else  // level is in the superlevel
    {
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

auto calculate_partfunct(int element, int ion, int modelgridindex) -> double
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
        "[fatal] calculate_sahafact: Negative Saha factor. sfac %g element %d ion %d level %d upperionlevel %d g_lower "
        "%g g_upper %g T %g E_threshold %g exppart %g\n",
        sf, element, ion, level, upperionlevel, g_lower, g_upper, T, E_threshold, exp(E_threshold / KB / T));
    abort();
  }
  return sf;
}

auto ionstagepop(int modelgridindex, int element, int ion) -> double
/// Calculates the given ionstages total population in nebular approximation for modelgridindex
/// The precalculated ground level population and partition function are used.
{
  return get_groundlevelpop(modelgridindex, element, ion) *
         grid::modelgrid[modelgridindex].composition[element].partfunct[ion] / stat_weight(element, ion, 0);
}
