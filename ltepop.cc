#include "sn3d.h"
#include "atomic.h"
#include "grid.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "nltepop.h"
#include "ratecoeff.h"
#include "update_grid.h"

// default option if not specified
#ifndef LTEPOP_EXCITATIONTEMPERATURE
  #define LTEPOP_EXCITATIONTEMPERATURE grid::get_Te(modelgridindex)
#endif

extern __host__ __device__ inline double calculate_sahafact(int element, int ion, int level, int upperionlevel, double T, double E_threshold);
extern __host__ __device__ inline double ionstagepop(int modelgridindex, int element, int ion);


double nne_solution_f(double x, void *paras)
/// For libgsl bracketing type solver
/// provides the equation which has to be solved to obtain the electron number
/// density (passed by x)
{
  const int modelgridindex = ((nne_solution_paras *) paras)->cellnumber;
  const double rho = grid::get_rho(modelgridindex);

  double outersum = 0.;
  //printout("debug get_nelements() %d =========================\n",get_nelements());
  for (int element = 0; element < get_nelements(); element++)
  {
    const double abundance = grid::modelgrid[modelgridindex].composition[element].abundance;
    if (abundance > 0)
    {
      double innersum = 0.;
      //printout("debug get_nions (element %d) %d =========================\n",element,get_nions(element));
      //uppermost_ion = globals::elements[element].uppermost_ion;
      const int uppermost_ion = grid::get_elements_uppermost_ion(modelgridindex, element);

      double ionfractions[uppermost_ion + 1];
      get_ionfractions(element, modelgridindex, x, ionfractions, uppermost_ion);

      int ion;
      for (ion = 0; ion <= uppermost_ion; ion++)
      {
        //printout("debug element %d, ion %d, ionfract(element,ion,T,x) %g\n",element,ion,ionfractions[ion]);
        innersum += (get_ionstage(element, ion) - 1) * ionfractions[ion];
      }
      assert_always(std::isfinite(innersum));
      outersum += abundance / globals::elements[element].mass * innersum;
      if (!std::isfinite(outersum))
      {
        printout("nne_solution_f: element %d ion %d uppermostion %d abundance %g, mass %g\n",
                 element, ion, grid::get_elements_uppermost_ion(modelgridindex, element),
                 abundance, globals::elements[element].mass);
        printout("outersum %g\n",outersum);
        abort();
      }
    }
  }

  return rho * outersum - x;
}


void get_ionfractions(int element, int modelgridindex, double nne, double ionfractions[], int uppermost_ion)
// Calculate the fractions of an element's population in each ionization stage
// size of ionfractions array must be >= uppermostion + 1
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(uppermost_ion < get_nions(element));

  double nnionfactor[uppermost_ion + 1];
  nnionfactor[uppermost_ion] = 1;

  double denominator = 1.;

  for (int ion = uppermost_ion - 1; ion >= 0; ion--)
  {
    nnionfactor[ion] = nnionfactor[ion + 1] * nne * phi(element, ion, modelgridindex);
    denominator += nnionfactor[ion];
  }

  for (int ion = 0; ion <= uppermost_ion; ion++)
  {
    const double numerator = nnionfactor[ion];
    ionfractions[ion] = numerator / denominator;

    if (!std::isfinite(ionfractions[ion]))
    {
      if (modelgridindex != grid::get_npts_model())
      {
        printout("[warning] ionfract set to zero for ionstage %d of Z=%d in cell %d with T_e %g, T_R %g\n",
                 get_ionstage(element,ion), get_element(element), modelgridindex, grid::get_Te(modelgridindex), grid::get_TR(modelgridindex));
        //abort();
        ionfractions[ion] = 0;
      }
    }
    // printout("ionfract(%d,%d,%d,%g) = %g\n", element, ion, modelgridindex, nne, ionfractions[ion]);
  }
}


__host__ __device__
static double interpolate_ions_spontrecombcoeff(const int element, const int ion, const double T)
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_always(T >= MINTEMP);
  int lowerindex = floor(log(T / MINTEMP) / T_step_log);
  if (lowerindex < TABLESIZE - 1)
  {
    int upperindex = lowerindex + 1;
    double T_lower =  MINTEMP * exp(lowerindex * T_step_log);
    double T_upper =  MINTEMP * exp(upperindex * T_step_log);

    double f_upper = globals::elements[element].ions[ion].Alpha_sp[upperindex];
    double f_lower = globals::elements[element].ions[ion].Alpha_sp[lowerindex];

    return f_lower + (f_upper - f_lower) / (T_upper - T_lower) * (T - T_lower);
  }
  else
    return globals::elements[element].ions[ion].Alpha_sp[TABLESIZE-1];
}


double phi(const int element, const int ion, const int modelgridindex)
/// Calculates population ratio (a saha factor) of two consecutive ionisation stages
/// in nebular approximation phi_j,k* = N_j,k*/(N_j+1,k* * nne)
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));

  double phi = 0;

  //double Y_nt, ionpot_in;
  //int element_in, ion_in, nions_in;
  //double rate_use;

  const float T_e = grid::get_Te(modelgridindex);
  //double T_R = grid::get_TR(modelgridindex);

  //double W = globals::cell[cellnumber].W;

  /// Old ionisation formula
  //partfunct_ratio = globals::cell[cellnumber].composition[element].partfunct[ion]/globals::cell[cellnumber].composition[element].partfunct[ion+1];
  //phi = 1./W * sqrt(T_R/T_e) * partfunct_ratio * SAHACONST * pow(T_R,-1.5) * exp(ionpot/KB/T_R);

  /// New ionisation formula with zeta
  //zeta = interpolate_zeta(element,ion,T_e);
  //phi = 1./W * 1./(zeta+W*(1-zeta)) * sqrt(T_R/T_e) * partfunct_ratio * SAHACONST * pow(T_R,-1.5) * exp(ionpot/KB/T_R);

  /// Newest ionisation formula

#ifdef FORCE_LTE
  const bool use_lte_ratio = true;
#else
  const bool use_lte_ratio = (globals::initial_iteration || grid::modelgrid[modelgridindex].thick == 1);
#endif

  if (use_lte_ratio)
  {
    const double ionpot = epsilon(element, ion + 1, 0) - epsilon(element, ion, 0);
    //printout("ionpot for element %d, ion %d is %g\n", element, ion, ionpot / EV);
    const double partfunct_ratio = grid::modelgrid[modelgridindex].composition[element].partfunct[ion] / grid::modelgrid[modelgridindex].composition[element].partfunct[ion + 1];
    phi = partfunct_ratio * SAHACONST * pow(T_e, -1.5) * exp(ionpot / KB / T_e);
  }
  else
// elseif (NLTE_POPS_ALL_IONS_SIMULTANEOUS)
//     {
//       const float nne = grid::get_nne(modelgridindex);
//       phi = ionstagepop(modelgridindex,element,ion) / ionstagepop(modelgridindex,element,ion+1) / nne;
//     }
// else
  {
    //Gamma = photoionestimator[cellnumber*get_nelements()*get_max_nions()+element*get_max_nions()+ion];
    #if NO_LUT_PHOTOION
      const double Gamma = calculate_iongamma_per_gspop(modelgridindex, element, ion);
    #else
      const double Gamma = globals::gammaestimator[modelgridindex * get_nelements() * get_max_nions() + element * get_max_nions() + ion];
    #endif
    // printout("phicompare element %d ion %d T_e = %g gammaestimator %g calculate_iongamma_per_gspop %g\n",
    //          element, ion, T_e,
    //          globals::gammaestimator[modelgridindex * get_nelements() * get_max_nions() + element * get_max_nions() + ion],
    //          calculate_iongamma_per_gspop(modelgridindex, element, ion));

    // Gamma is the photoionization rate per ground level pop
    const double Gamma_ion = Gamma * stat_weight(element, ion, 0) / grid::modelgrid[modelgridindex].composition[element].partfunct[ion];

    if (Gamma == 0. && (!NT_ON || (globals::rpkt_emiss[modelgridindex] == 0. && grid::get_modelinitradioabund(modelgridindex, 24, 48) == 0. && grid::get_modelinitradioabund(modelgridindex, 28, 56) == 0.)))
    {
      printout("Fatal: Gamma = 0 for element %d, ion %d in phi ... abort\n",element,ion);
      abort();
    }

    //Alpha_st = stimrecombestimator[cellnumber*get_nelements()*get_max_nions()+element*get_max_nions()+ion];
    double Alpha_st = 0.; ///approximate treatment neglects stimulated recombination

    double Alpha_sp = 0.;
    if (NLTE_POPS_ON)
    {
      Alpha_sp = calculate_ionrecombcoeff(modelgridindex, T_e, element, ion + 1, false, false, false, false, false, false);
    }
    else
    {
      Alpha_sp = interpolate_ions_spontrecombcoeff(element, ion, T_e);
    }

    // const double Col_rec = calculate_ionrecombcoeff(modelgridindex, T_e, element, ion + 1, false, true, false, false, false);
    const double Col_rec = 0.;

    double Y_nt = 0.0;

    if (NT_ON)
    {
      Y_nt = nonthermal::nt_ionization_ratecoeff(modelgridindex, element, ion);
    }

    // || !std::isfinite(Gamma))
    //return phi_lte(element,ion,cellnumber);
    //gamma_lte = interpolate_photoioncoeff_below(element,ion,0,T_e) + interpolate_photoioncoeff_above(element,ion,0,T_e);
    //zeta = interpolate_zeta(element,ion,T_e);
    //alpha_sp = interpolate_spontrecombcoeff(element,ion,0,T_e);
    //phi = gamma_lte*(Alpha_sp+Apha_st)/(Gamma*alpha_sp) * partfunct_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e);

    phi = (Alpha_sp + Alpha_st + Col_rec) / (Gamma_ion + Y_nt);

    // Y_nt should generally be higher than the Gamma term for nebular epoch

    if (!std::isfinite(phi) || phi == 0.)
    {
      printout("[fatal] phi: phi %g exceeds numerically possible range for element %d, ion %d, T_e %g ... remove higher or lower ionisation stages\n", phi, element, ion, T_e);
      printout("[fatal] phi: Alpha_sp %g, Alpha_st %g, Gamma %g, partfunct %g, stat_weight %g\n", Alpha_sp, Alpha_st, Gamma, grid::modelgrid[modelgridindex].composition[element].partfunct[ion], stat_weight(element, ion, 0));
      printout("[fatal] phi: upperionpartfunct %g, upperionstatweight %g\n", grid::modelgrid[modelgridindex].composition[element].partfunct[ion + 1], stat_weight(element, ion + 1, 0));
      printout("[fatal] phi: Y_nt %g Col_rec %g grid::get_nne(modelgridindex) %g\n", Y_nt, Col_rec, grid::get_nne(modelgridindex));
      //abort();
    }
  }

  return phi;
}


//double phi_lte(int element, int ion, int cellnumber)
/// Calculates population ratio (a saha factor) of two consecutive ionisation stages
/// in nebular approximation phi_j,k* = N_j,k*/N_j+1,k* * nne
/*{
  double partfunct_ratio;
  double phi;

  double ionpot = epsilon(element,ion+1,0) - epsilon(element,ion,0);
  float T_e = globals::cell[cellnumber].T_e;

  partfunct_ratio = globals::cell[cellnumber].composition[element].partfunct[ion]/globals::cell[cellnumber].composition[element].partfunct[ion+1];
  phi = partfunct_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e);

  if (!std::isfinite(phi))
  {
    printout("phi_lte: phi %g exceeds numerically possible range for element %d, ion %d, T_e %g, ... remove higher or lower ionisation stages\n",phi,element,ion,T_e);
    abort();
  }
  return phi;
}*/


/*
double calculate_ltepartfunct(int element, int ion, double T)
/// Calculates the LTE partition function for ion=ion of element=element at
/// temperature T
{
  double U;
  double epsilon_groundlevel;
  double oneoverkbt;
  int level;
  int nlevels;

  epsilon_groundlevel = epsilon(element,ion,0);
  oneoverkbt = 1/KB/T;
  U = 0.;
  nlevels = get_nlevels(element,ion);
  for (level = 0; level < nlevels; level++)
  {
    U += stat_weight(element,ion,level) * exp(-(epsilon(element,ion,level)-epsilon_groundlevel)*oneoverkbt);
  }

  if (!std::isfinite(U)) abort();
  return U;
}
*/


double calculate_partfunct(int element, int ion, int modelgridindex)
/// Calculates the partition function for ion=ion of element=element in
/// cell modelgridindex
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  double pop_store;
  //double E_level, E_ground, test;

  //double T_exc = grid::get_TJ(modelgridindex);
  //double W = 1.;

  //  T_exc = MINTEMP;

/*  if (T_exc <= MINTEMP)
  {
    T_exc = grid::get_TR(modelgridindex);
    W = grid::get_W(modelgridindex);
  }*/
  //double T_R = grid::modelgrid[modelgridindex].T_R;
  //double W = grid::modelgrid[modelgridindex].W;
  //double T_J = pow(W,1./4.)*T_R;
  //double oneoverkbtexc = 1/KB/T_exc;
  //double epsilon_groundlevel = epsilon(element,ion,0);

  int initial = 0;
  if (get_groundlevelpop(modelgridindex, element, ion) < MINPOP)
  {
    //either there reall is none of this ion or this is a first pass through
    //in either case, we won't have any real nlte_populations so the actual value of
    //of groundlevelpop for this calculation doesn't matter, so long as it's not zero!
    pop_store = get_groundlevelpop(modelgridindex, element, ion);
    initial = 1;
    grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = 1.0;
  }

  //printout("groundlevelpop %g\n", get_groundlevelpop(modelgridindex,element,ion));

  double U = 1.0;//stat_weight(element,ion,0);
  const int nlevels = get_nlevels(element, ion);
/*  if (T_exc <= MINTEMP)
  {
    for (level = 1; level < nlevels; level++)
    {
      if (globals::elements[element].ions[ion].levels[level].metastable)
      {
        T_exc = grid::get_TJ(modelgridindex);
        W = 1.;
      }
      else
      {
        T_exc = grid::get_TR(modelgridindex);
        W = grid::get_W(modelgridindex);
      }
      U += W*stat_weight(element,ion,level) * exp(-(epsilon(element,ion,level)-epsilon_groundlevel)*oneoverkbtexc);
    }
  }
  else
  {*/
    const double groundpop = get_groundlevelpop(modelgridindex, element, ion);
    for (int level = 1; level < nlevels; level++)
    {
      const double nn = calculate_exclevelpop(modelgridindex, element, ion, level) / groundpop;//*stat_weight(element,ion,0);

      //if (NLTE_POPS_ON)
      //{
      //if ((is_nlte(element,ion,level) != 1) || (test = grid::modelgrid[modelgridindex].nlte_pops[globals::elements[element].ions[ion].first_nlte+level-1]) < -0.9)
      //{
	  /* Case for when no NLTE level information is available yet */
      //  E_level = epsilon(element,ion,level);
      //  E_ground = epsilon(element,ion,0);
      //  nn = W * stat_weight(element,ion,level) * exp(-(E_level-E_ground)/KB/T_exc);
      //	}
      //else
      //	{
	  //printout("Using an nlte population!\n");
      //	  nn = test * grid::modelgrid[modelgridindex].rho / get_groundlevelpop(modelgridindex,element,ion)*stat_weight(element,ion,0);
      //	  if (!std::isfinite(nn))
      //	    {
      //
      //      printout("[fatal] NLTE population failure.\n");
      //	      printout("element %d ion %d level %d\n", element, ion, level);
      //	      printout("nn %g test %g rho %g\n", nn, test, grid::modelgrid[modelgridindex].rho);
      //	      printout("ground level %g\n", get_groundlevelpop(modelgridindex,element,ion));
      //	      abort();
      //	    }

      //	}
      //}
	     U += nn;//W*stat_weight(element,ion,level) * exp(-(epsilon(element,ion,level)-epsilon_groundlevel)*oneoverkbtexc);
    }
//   }

  U *= stat_weight(element,ion,0);

  if (!std::isfinite(U))
  {
    printout("element %d ion %d\n",element,ion);
    printout("modelgridindex %d\n",modelgridindex);
    printout("nlevels %d\n",nlevels);
    printout("sw %g\n",stat_weight(element,ion,0));
    //printout("T_exc %g \n",T_exc);
    abort();
  }

  if (initial == 1)
  {
    //put back the zero, just in case it matters for something
    grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = pop_store;
  }

  return U;
}


/*
double calculate_groundlevelpop(int element, int ion, double T, int cellnumber, double nne, double nnnextion)
///calculates ground level population for ion=ion of element=element at
///temperature T and electron number density nne
///further the total population number nnnextion of the next higher ionisation stage is needed
{
  double partfunct(int element, int ion, double T);

  double deltaE = epsilon(element,ion+1,0) - epsilon(element,ion,0);
  double n0;

  //n0 = nnnextion * nne * stat_weight(element,ion,0)/partfunct(element,ion+1,T) * C * pow(T,-1.5) * exp(deltaE/KB/T);
  n0 = nnnextion * nne * stat_weight(element,ion,0)/globals::cell[cellnumber].composition[element].partfunct[ion+1] * SAHACONST * pow(T,-1.5) * exp(deltaE/KB/T);

  return n0;
}
*/


__host__ __device__
double get_groundlevelpop(int modelgridindex, int element, int ion)
/// Returns the given ions groundlevel population for modelgridindex which was precalculated
/// during update_grid and stored to the grid.
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  //double nn = grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion];
  //if (nn < MINPOP) nn = MINPOP;
  //return nn;

  const double nn = grid::modelgrid[modelgridindex].composition[element].groundlevelpop[ion];
  if (nn < MINPOP)
  {
    if (grid::get_elem_abundance(modelgridindex,element) > 0)
      return MINPOP;
    else
      return 0.;
  }
  else
  {
    return nn;
  }
}


__host__ __device__
double get_groundmultiplet_pop(int modelgridindex, int element, int ion)
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  const int nlevels_gm = get_nlevels_groundterm(element, ion);

  double gmpop = 0.;
  for (int level = 0; level < nlevels_gm; level++)
  {
    gmpop += calculate_exclevelpop(modelgridindex, element, ion, level);
  }

  return gmpop;
}


/*static double calculate_exclevelpop_old(int modelgridindex, int element, int ion, int level)
/// Calculates occupation number of level relative to the ions ground level population
/// using a modified version of the Boltzmann formula, which fulfills the diluted BB
/// approximation (or nebular approximation).
{
  double T_exc = grid::get_TJ(modelgridindex);
  double W = 1.;

  //  T_exc = MINTEMP;


  // if (T_exc <= MINTEMP)
  // {
  //   if (globals::elements[element].ions[ion].levels[level].metastable)
  //   {
  //     T_exc = grid::get_TJ(modelgridindex);
  //   }
  //   else
  //   {
  //     T_exc = grid::get_TR(modelgridindex);
  //     W = grid::get_W(modelgridindex);
  //   }
  // }

  double nn;
  if (level == 0)
    nn = get_groundlevelpop(modelgridindex,element,ion);
  else
  {
    double E_level = epsilon(element,ion,level);
    double E_ground = epsilon(element,ion,0);
    nn = get_groundlevelpop(modelgridindex,element,ion) * W * stat_weight(element,ion,level) / stat_weight(element,ion,0) * exp(-(E_level-E_ground)/KB/T_exc);
  }

  if (nn < MINPOP)
  {
    if (grid::get_elem_abundance(modelgridindex,element) > 0)
      nn = MINPOP;
    else
      nn = 0.;
  }

  if (!std::isfinite(nn))
  {
    printout("[fatal] calculate_exclevelpop: level %d of ion %d of element %d has infinite level population %g\n",level,ion,element,nn);
    printout("[fatal] calculate_exclevelpop: associated ground level has pop %g\n",get_groundlevelpop(modelgridindex,element,ion));
    printout("[fatal] calculate_exclevelpop: associated ion has pop %g\n",ionstagepop(modelgridindex,element,ion));
    printout("[fatal] calculate_exclevelpop: associated partition function %g\n",grid::modelgrid[modelgridindex].composition[element].partfunct[ion]);
  }
  return nn;
}*/


__host__ __device__
double calculate_levelpop_lte(int modelgridindex, int element, int ion, int level)
/// Calculates occupation population of a level assuming LTE excitation
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  if (level == 0)
    return get_groundlevelpop(modelgridindex, element, ion);

  const double T_exc = LTEPOP_EXCITATIONTEMPERATURE;
  const double W = 1.;

  const double E_level = epsilon(element, ion, level);
  const double E_ground = epsilon(element, ion, 0);
  const double nnground = get_groundlevelpop(modelgridindex, element, ion);

  return (
    nnground * W *
    stat_weight(element, ion, level) / stat_weight(element, ion, 0) *
    exp(- (E_level - E_ground) / KB / T_exc));
}


__host__ __device__
double calculate_exclevelpop(int modelgridindex, int element, int ion, int level)
/// Calculates the population of a level from either LTE or NLTE information
{
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  double nn;

  //  T_exc = MINTEMP;


/*  if (T_exc <= MINTEMP)
  {
    if (globals::elements[element].ions[ion].levels[level].metastable)
    {
      T_exc = grid::get_TJ(modelgridindex);
    }
    else
    {
      T_exc = grid::get_TR(modelgridindex);
      W = grid::get_W(modelgridindex);
    }
  }*/

  if (level == 0)
  {
    nn = get_groundlevelpop(modelgridindex,element,ion);
  }
  else if (NLTE_POPS_ON)
  {
    if (is_nlte(element,ion,level))
    {
      //printout("Using an nlte population!\n");
      const double nltepop_over_rho = grid::modelgrid[modelgridindex].nlte_pops[globals::elements[element].ions[ion].first_nlte + level - 1];
      if (nltepop_over_rho < -0.9)
      {
        // Case for when no NLTE level information is available yet
        nn = calculate_levelpop_lte(modelgridindex, element, ion, level);
      }
      else
      {
        //printout("Using an nlte population!\n");
        nn = nltepop_over_rho * grid::get_rho(modelgridindex);
        if (!std::isfinite(nn))
        {
          printout("[fatal] NLTE population failure.\n");
          printout("element %d ion %d level %d\n", element, ion, level);
          printout("nn %g nltepop_over_rho %g rho %g\n", nn, nltepop_over_rho, grid::get_rho(modelgridindex));
          printout("ground level %g\n", get_groundlevelpop(modelgridindex, element, ion));
          abort();
        }
        return nn;
      }
    }
    else // level is in the superlevel
    {
      assert_always(level_isinsuperlevel(element, ion, level));

      const int sl_nlte_index = globals::elements[element].ions[ion].first_nlte + get_nlevels_nlte(element, ion);

      const double superlevelpop_over_rho = grid::modelgrid[modelgridindex].nlte_pops[sl_nlte_index];
      if (superlevelpop_over_rho < -0.9) //TODO: should change this to less than zero?
      {
        // Case for when no NLTE level information is available yet
        nn = calculate_levelpop_lte(modelgridindex, element, ion, level);
      }
      else
      {
        //printout("Using a superlevel population!\n");
        nn = superlevelpop_over_rho * grid::get_rho(modelgridindex) * superlevel_boltzmann(modelgridindex, element, ion, level);
        if (!std::isfinite(nn))
        {
          printout("[fatal] NLTE population failure.\n");
          printout("element %d ion %d level %d\n", element, ion, level);
          printout("nn %g superlevelpop_over_rho %g rho %g\n", nn, superlevelpop_over_rho, grid::get_rho(modelgridindex));
          printout("ground level %g\n", get_groundlevelpop(modelgridindex, element, ion));
          abort();
        }
        return nn;
      }
    }
  }
  else
  {
    nn = calculate_levelpop_lte(modelgridindex, element, ion, level);
  }

  if (nn < MINPOP)
  {
    if (grid::get_elem_abundance(modelgridindex,element) > 0)
      nn = MINPOP;
      // nn = calculate_levelpop_lte(modelgridindex, element, ion, level);
    else
      nn = 0.;
  }

  if (!std::isfinite(nn))
  {
    printout("[fatal] calculate_exclevelpop: level %d of ion %d of element %d has infinite level population %g\n",level,ion,element,nn);
    printout("[fatal] calculate_exclevelpop: associated ground level has pop %g\n", get_groundlevelpop(modelgridindex, element, ion));
    printout("[fatal] calculate_exclevelpop: associated ion has pop %g\n", ionstagepop(modelgridindex, element, ion));
    printout("[fatal] calculate_exclevelpop: associated partition function %g\n",grid::modelgrid[modelgridindex].composition[element].partfunct[ion]);
  }

  return nn;
}

/*
#ifdef FORCE_LTE
{
  double nn;

  double T_R = grid::get_TR(modelgridindex);
  //double W = globals::cell[cellnumber].W;

  if (level == 0) nn = get_groundlevelpop(modelgridindex,element,ion);
  else nn = get_groundlevelpop(modelgridindex,element,ion) * stat_weight(element,ion,level)/stat_weight(element,ion,0) * exp(-(epsilon(element,ion,level)-epsilon(element,ion,0))/KB/T_R);

  return nn;
}
#else
{
  double E_level,E_ground;
  double nn;

  double T_J = grid::get_TJ(modelgridindex);
  double T_R = grid::get_TR(modelgridindex);
  double W = grid::get_W(modelgridindex);
  //double T_J = pow(W,1./4.)*T_R;

  if (level == 0) nn = get_groundlevelpop(modelgridindex,element,ion);
  else
  {
    E_level = epsilon(element,ion,level);
    E_ground = epsilon(element,ion,0);
    nn = get_groundlevelpop(modelgridindex,element,ion) * stat_weight(element,ion,level)/stat_weight(element,ion,0) * exp(-(E_level-E_ground)/KB/T_J);
  }

  if (nn < MINPOP) nn = MINPOP;
  return nn;
}
#endif
*/


/*void calculate_levelpops(int modelgridindex)
/// Calculates the full level populations for a given grid cell
/// and stores them to the active entry of the cellhistory.
{
  for (int element = 0; element < get_nelements(); element++)
  {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      globals::cellhistory[tid].chelements[element].chions[ion].chlevels[0].population = get_groundlevelpop(modelgridindex, element, ion);
      //printout("element %d, ion %d, level 0: population %g\n",element,ion,groundlevelpop(cellnumber, element, ion));
      const int nlevels = get_nlevels(element,ion);
      for (int level = 1; level < nlevels; level++)
      {
        globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level].population = calculate_exclevelpop(modelgridindex,element,ion,level);
        //printout("element %d, ion %d, level %d: population %g\n",element,ion,level,exclevelpop(cellnumber,element,ion,level,T));
      }
    }
  }
}*/


