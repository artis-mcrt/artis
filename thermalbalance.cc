#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "thermalbalance.h"
#include "update_grid.h"


typedef struct Te_solution_paras
{
  double t_current;
  int modelgridindex;
  heatingcoolingrates_t *heatingcoolingrates;
} Te_solution_paras;

typedef struct
{
  double nu_edge;
  int modelgridindex;
  float T_R;
  float *photoion_xs;
} gsl_integral_paras_bfheating;


#if (!NO_LUT_BFHEATING)
  double get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, float T, double W)
  {
    /// The correction factor for stimulated emission in gammacorr is set to its
    /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
    /// correction may be evaluated at T_R!
    double bfheatingcoeff = 0.;

    /*double nnlevel = calculate_exclevelpop(cellnumber,element,ion,level);
    bfheating = nnlevel * W * interpolate_bfheatingcoeff_below(element,ion,level,T_R);*/
    const int lowerindex = floor(log(T/MINTEMP)/T_step_log);
    if (lowerindex < TABLESIZE - 1)
    {
      const int upperindex = lowerindex + 1;
      const double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
      const double T_upper =  MINTEMP * exp(upperindex*T_step_log);

      const double f_upper = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfheating_coeff[upperindex];
      const double f_lower = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfheating_coeff[lowerindex];

      bfheatingcoeff = (f_lower + (f_upper - f_lower)/(T_upper - T_lower) * (T - T_lower));
    }
    else
      bfheatingcoeff = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfheating_coeff[TABLESIZE-1];

    return W * bfheatingcoeff;
  }
#endif


__host__ __device__
static double integrand_bfheatingcoeff_custom_radfield(double nu, void *voidparas)
/// Integrand to calculate the rate coefficient for bfheating using gsl integrators.
{
  const gsl_integral_paras_bfheating *const params = (gsl_integral_paras_bfheating *) voidparas;

  const int modelgridindex = params->modelgridindex;
  const double nu_edge = params->nu_edge;
  const double T_R = params->T_R;
  // const double Te_TR_factor = params->Te_TR_factor; // = sqrt(T_e/T_R) * sahafac(Te) / sahafac(TR)

  const float sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge, nu);

  // const float T_e = get_Te(modelgridindex);
  // return sigma_bf * (1 - nu_edge/nu) * radfield(nu,modelgridindex) * (1 - Te_TR_factor * exp(-HOVERKB * nu / T_e));

  return sigma_bf * (1 - nu_edge/nu) * radfield(nu,modelgridindex) * (1 - exp(-HOVERKB*nu/T_R));
}


static double calculate_bfheatingcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex)
{
  // const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
  // const double E_threshold = epsilon(element, ion + 1, upperionlevel) - epsilon(element, ion, level);
  const double E_threshold = get_phixs_threshold(element, ion, level, phixstargetindex);

  const double nu_threshold = ONEOVERH * E_threshold;
  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge; // nu of the uppermost point in the phixs table

  // const float T_e = get_Te(modelgridindex);
  // const double T_R = get_TR(modelgridindex);
  // const double sf_Te = calculate_sahafact(element,ion,level,upperionlevel,T_e,E_threshold);
  // const double sf_TR = calculate_sahafact(element,ion,level,upperionlevel,T_R,E_threshold);

  gsl_integral_paras_bfheating intparas;
  intparas.nu_edge = nu_threshold;
  intparas.modelgridindex = modelgridindex;
  intparas.T_R = get_TR(modelgridindex);
  intparas.photoion_xs = elements[element].ions[ion].levels[level].photoion_xs;
  // intparas.Te_TR_factor = sqrt(T_e/T_R) * sf_Te / sf_TR;

  double bfheating = 0.0;

#if (!CUDA_ENABLED || !USECUDA_BFHEATING || CUDA_VERIFY_CPUCONSISTENCY)
  double error = 0.0;
  const double epsrel = 1e-3;
  const double epsrelwarning = 1e-1;
  const double epsabs = 0.;

  gsl_function F_bfheating;
  F_bfheating.function = &integrand_bfheatingcoeff_custom_radfield;
  F_bfheating.params = &intparas;

  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);

  const int status = gsl_integration_qag(
    &F_bfheating, nu_threshold, nu_max_phixs, epsabs, epsrel,
     GSLWSIZE, GSL_INTEG_GAUSS61, gslworkspace, &bfheating, &error);

  if (status != 0 && (status != 18 || (error / bfheating) > epsrelwarning))
  {
    printout("bf_heating integrator gsl warning %d. modelgridindex %d Z=%d ionstage %d lower %d phixstargetindex %d integral %g error %g\n",
             status, modelgridindex, get_element(element), get_ionstage(element, ion), level, phixstargetindex, bfheating, error);
  }
  gsl_set_error_handler(previous_handler);
#endif

#if (CUDA_ENABLED && USECUDA_BFHEATING)

  void *dev_intparas;
  checkCudaErrors(cudaMalloc(&dev_intparas, sizeof(gsl_integral_paras_bfheating)));
  checkCudaErrors(cudaMemcpy((void**) dev_intparas, (void *) &intparas, sizeof(gsl_integral_paras_bfheating), cudaMemcpyHostToDevice));

  const double bfheating_gpu = calculate_integral_gpu<integrand_bfheatingcoeff_custom_radfield>(dev_intparas, nu_threshold, nu_max_phixs);

  cudaFree(dev_intparas);

  #if CUDA_VERIFY_CPUCONSISTENCY
  if (bfheating > 0 && abs(bfheating_gpu / bfheating - 1.) > 0.3)
  {
    printout("WARNING bfheating CUDA test: element %d ion %d level %d phixstargetindex %d modelgridindex %d gsl %.2e gpu %.2e\n", element, ion, level, phixstargetindex, modelgridindex, bfheating, bfheating_gpu);
  }
  #endif

  bfheating = bfheating_gpu;

#endif

  bfheating *= FOURPI * get_phixsprobability(element, ion, level, phixstargetindex);

  return bfheating;
}


static double get_bfheatingcoeff(int element, int ion, int level, int phixstargetindex)
// depends only the radiation field
// no dependence on T_e or populations
{
  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].bfheatingcoeff;
}


void calculate_bfheatingcoeffs(int modelgridindex)
{
  for (int element = 0; element < nelements; element++)
  {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      const int nlevels = get_nlevels(element,ion);
      for (int level = 0; level < nlevels; level++)
      {
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
        {
        #if NO_LUT_BFHEATING
          const double bfheatingcoeff = calculate_bfheatingcoeff(element, ion, level, phixstargetindex, modelgridindex);
        #else
          /// The correction factor for stimulated emission in gammacorr is set to its
          /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
          /// correction may be evaluated at T_R!
          const double T_R = get_TR(modelgridindex);
          const double W = get_W(modelgridindex);
          double bfheatingcoeff = W * get_bfheatingcoeff_ana(element,ion,level,phixstargetindex, T_R, W);
          const int index_in_groundlevelcontestimator = elements[element].ions[ion].levels[level].closestgroundlevelcont;
          if (index_in_groundlevelcontestimator >= 0)
            bfheatingcoeff *= bfheatingestimator[modelgridindex*nelements*maxion + index_in_groundlevelcontestimator];

          if (!isfinite(bfheatingcoeff))
          {
            printout("[fatal] get_bfheatingcoeff returns a NaN! W %g get_bfheatingcoeff(element,ion,level,phixstargetindex) %g index_in_groundlevelcontestimator %d bfheatingestimator[modelgridindex*nelements*maxion+index_in_groundlevelcontestimator] %g",
                     W,get_bfheatingcoeff(element,ion,level,phixstargetindex),index_in_groundlevelcontestimator,bfheatingestimator[modelgridindex*nelements*maxion+index_in_groundlevelcontestimator]);
            abort();
          }
        #endif

          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].bfheatingcoeff = bfheatingcoeff;
        }
      }
    }
  }
  cellhistory[tid].bfheating_mgi = modelgridindex;
}


static double get_heating_ion_coll_deexc(const int modelgridindex, const int element, const int ion, const double T_e, const double nne)
{
  double C_deexc = 0.;
  const int nlevels = get_nlevels(element, ion);

  for (int level = 0; level < nlevels; level++)
  {
    const double nnlevel = calculate_exclevelpop(modelgridindex, element, ion, level);
    const double epsilon_level = epsilon(element, ion, level);
    // Collisional heating: deexcitation to same ionization stage
    // ----------------------------------------------------------
    const int ndowntrans = get_ndowntrans(element, ion, level);
    for (int ii = 0; ii < ndowntrans; ii++)
    {
      const int lineindex = elements[element].ions[ion].levels[level].downtrans_lineindicies[ii];
      const int lower = linelist[lineindex].lowerlevelindex;
      const double epsilon_trans = epsilon_level - epsilon(element, ion, lower);
      const double C = nnlevel * col_deexcitation_ratecoeff(T_e, nne, epsilon_trans, lineindex) * epsilon_trans;
      C_deexc += C;
    }
  }
  // const double nnion = ionstagepop(modelgridindex, element, ion);
  // printout("ion_col_deexc_heating: T_e %g nne %g Z=%d ionstage %d nnion %g heating_contrib %g contrib/nnion %g\n", T_e, nne, get_element(element), get_ionstage(element, ion), nnion, C_deexc, C_deexc / nnion);
  return C_deexc;
}


static void calculate_heating_rates(
  const int modelgridindex, const double T_e, const double nne, heatingcoolingrates_t *heatingcoolingrates)
/// Calculate the heating rates for a given cell. Results are returned
/// via the elements of the global heatingrates data structure.
{
  /*
  gsl_function F_bfheating;
  F_bfheating.function = &bfheating_integrand_gsl;
  gslintegration_bfheatingparas bfheatingparas;
  bfheatingparas.cellnumber = cellnumber;
  double bfhelper;
  */

/*  PKT dummypkt;
  dummypkt.where = cellnumber;
  PKT *pkt_ptr;
  pkt_ptr = &dummypkt;*/

  double C_deexc = 0.;
  //double C_recomb = 0.;
  double bfheating = 0.;
  double ffheating = 0.;

  assert(cellhistory[tid].bfheating_mgi == modelgridindex);

  //int nlevels_lowerion = 0;
  for (int element = 0; element < nelements; element++)
  {
    // mastate[tid].element = element;
    #ifdef DIRECT_COL_HEAT
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      C_deexc += get_heating_ion_coll_deexc(modelgridindex, element, ion, T_e, nne);
    }
    #endif
//
//         /// Collisional heating: recombination to lower ionization stage
//         /// ------------------------------------------------------------
//         /// For the moment we deal only with ionisations to the next ions groundlevel.
//         /// For speed issues (reduced number of calls to epsilon) this is checked here
//         /// instead of the more general way to check in col_recomb!
//         if (ion > 0 && level == 0) /// Check whether lower ionisation stage available
//         {
//           for (lower = 0; lower < nlevels_lowerion; lower++)
//           {
//             epsilon_trans = epsilon_current - epsilon(element,ion-1,lower);
//             C = col_recombination(pkt_ptr,lower,epsilon_trans)*epsilon_trans;
//             C_recomb += C;
//           }
//         }
//
//
//         /*
//         /// Bound-free heating (analytical calculation)
//         /// -------------------------------------------
//         /// We allow bound free-transitions only if there is a higher ionisation stage
//         /// left in the model atom to match the bound-free absorption in the rpkt routine.
//         /// There this condition is needed as we can only ionise to existing ionisation
//         /// stage even if there would be further ionisation stages in nature which
//         /// are not included in the model atom.
//         if (ion < nions-1)
//         {
//           epsilon_trans = epsilon(element,ion+1,0) - epsilon_current;
//           /// NB: W comes from the fact, that the W coming from the radiation field was factored
//           /// out in the precalculation of the bf-heating coefficient (this is justified by the
//           /// linear dependence on W).
//           /// The rate coefficient is calculated under the assumption T_e=T_R because its direct
//           /// T_e dependence is very weak. This means we have to pass T_R as the temperature
//           /// even if we are iterating here on T_e. (Otherwise we would allow a large temperature
//           /// range for T_R which changes the coefficient strongly).
//           //C = interpolate_bfheatingcoeff(element,ion,level,T_R)*W*nnlevel;
//           C = nnlevel * (W*interpolate_bfheatingcoeff_below(element,ion,level,T_R));// + W_D*interpolate_bfheatingcoeff_above(element,ion,level,T_D));
//           if (element == 6)
//           {
//             if (ion == 1)
//             {
//               printout("T_e %g, nnlevel %g, W %g, heatingcoeff %g\n",T_e,nnlevel,W,interpolate_bfheatingcoeff_below(element,ion,level,T_R));
//             }
//           }
//
//           /// Exact calculation of the bf-heating coefficients using integrators.
//           /// This makes things SLOW!!!
//           //bfheatingparas.nu_edge = epsilon_trans/H;
//           //F_bfheating.params = &bfheatingparas;
//           /// Discuss about the upper frequency limit (here 1e16 Hz) which we should choose
//           //gsl_integration_qag(&F_bfheating, bfheatingparas.nu_edge, 10*bfheatingparas.nu_edge, 0, intaccuracy, 1024, 6, wspace, &bfhelper, &error);
//           //C = bfhelper * FOURPI * nnlevel;
//           bfheating += C;
//         }
//         */
//       }

        /// Bound-free heating (from estimators)
        /// ------------------------------------
        //if (ion < nions-1) bfheating += bfheatingestimator[cellnumber*nelements*maxion+element*maxion+ion];

        /// Bound-free heating (renormalised analytical calculation)
        /// --------------------------------------------------------
        /// We allow bound free-transitions only if there is a higher ionisation stage
        /// left in the model atom to match the bound-free absorption in the rpkt routine.
        /// There this condition is needed as we can only ionise to existing ionisation
        /// stage even if there would be further ionisation stages in nature which
        /// are not included in the model atom.
        //nlevels_currention = get_nlevels(element,ion);
        //nlevels_currention = get_ionisinglevels(element,ion);
  //      if (ion > 0) nlevels_lowerion = get_nlevels(element,ion-1);

    for (int ion = 0; ion < nions - 1; ion++)
    {
      const int nbflevels = get_ionisinglevels(element, ion);
      for (int level = 0; level < nbflevels; level++)
      {
        const double nnlevel = calculate_exclevelpop(modelgridindex,element,ion,level);
        const int nphixstargets = get_nphixstargets(element,ion,level);
        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++)
        {
          //int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
          //double epsilon_upper = epsilon(element,ion+1,upper);
          //double epsilon_trans = epsilon_upper - epsilon_current;
          bfheating += nnlevel * get_bfheatingcoeff(element, ion, level, phixstargetindex);
        }
      }
    }
  }


  /// Free-free heating
  /// -----------------
  /// From estimators
  ffheating = ffheatingestimator[modelgridindex];

  /// Analytical calculation using T, and populations
  /// This is always taken as an additional process to the other importantheatingterms
  /// because its not done per ion but overall ions.
  /*
  gsl_function F_ffheating;
  F_ffheating.function = &ffheating_integrand_gsl;

  gslintegration_ffheatingparas intparas;
  intparas.T_e = T_e;
  intparas.cellnumber = cellnumber;
  F_ffheating.params = &intparas;

  /// Discuss about the upper frequency limit (here 1e16 Hz) which we should choose
  gsl_integration_qag(&F_ffheating, 0, 1e16, 0, intaccuracy, 1024, GSL_INTEG_GAUSS61, wspace, &ffheating, &error);
  /// or this integrator
  //gsl_integration_qng(&F_ffheating, 0, 1e16, 0, intaccuracy, &ffheating, &error, &neval);
  ffheating *= FOURPI;
  */

#ifdef DIRECT_COL_HEAT
  heatingcoolingrates->heating_collisional = C_deexc;
#else
  /// Collisional heating (from estimators)
  heatingcoolingrates->heating_collisional = colheatingestimator[modelgridindex];//C_deexc + C_recomb;
#endif

  heatingcoolingrates->heating_bf = bfheating;
  heatingcoolingrates->heating_ff = ffheating;
}


static double T_e_eqn_heating_minus_cooling(const double T_e, void *paras)
/// Thermal balance equation on which we have to iterate to get T_e
{
  const int modelgridindex = ((Te_solution_paras *) paras)->modelgridindex;
  const double t_current = ((Te_solution_paras *) paras)->t_current;
  heatingcoolingrates_t *heatingcoolingrates = ((Te_solution_paras *) paras)->heatingcoolingrates;

  /// Set new T_e guess for the current cell and update populations
  //cell[cellnumber].T_e = T_e;
  set_Te(modelgridindex, T_e);
  double nntot;
  if (NLTE_POPS_ON && NLTE_POPS_ALL_IONS_SIMULTANEOUS)
    nntot = calculate_electron_densities(modelgridindex);
  else
    nntot = calculate_populations(modelgridindex);

  /// Then calculate heating and cooling rates
  const float nne = get_nne(modelgridindex);
  calculate_cooling_rates(modelgridindex, heatingcoolingrates);
  calculate_heating_rates(modelgridindex, T_e, nne, heatingcoolingrates);

  /// If selected take direct gamma heating into account
  if (do_rlc_est == 3)
  {
    const double nt_frac_heating = get_nt_frac_heating(modelgridindex);
    heatingcoolingrates->heating_gamma = get_deposition_rate_density(modelgridindex) * nt_frac_heating;
    heatingcoolingrates->nt_frac_heating = nt_frac_heating;
  }
  else
  {
    heatingcoolingrates->heating_gamma = 0.;
  }

  /// Adiabatic cooling term
  const double p = nntot * KB * T_e;
  const double volumetmin = vol_init_modelcell(modelgridindex);
  const double dV = 3 * volumetmin / pow(tmin, 3) * pow(t_current, 2);
  const double V = volumetmin * pow(t_current / tmin, 3);
  //printout("nntot %g, p %g, dV %g, V %g\n",nntot,p,dV,V);
  heatingcoolingrates->cooling_adiabatic = p * dV / V;

  const double total_heating_rate = heatingcoolingrates->heating_ff + heatingcoolingrates->heating_bf + heatingcoolingrates->heating_collisional + heatingcoolingrates->heating_gamma;
  const double total_coolingrate = heatingcoolingrates->cooling_ff + heatingcoolingrates->cooling_fb + heatingcoolingrates->cooling_collisional + heatingcoolingrates->cooling_adiabatic;

  return total_heating_rate - total_coolingrate; // - 0.01*(heatingrates_thisthread->bf+coolingrates[tid].fb)/2;
}


void call_T_e_finder(const int modelgridindex, const int timestep, const double t_current, const double T_min, const double T_max, heatingcoolingrates_t *heatingcoolingrates)
{
  const double T_e_old = get_Te(modelgridindex);
  printout("Finding T_e in cell %d at timestep %d...", modelgridindex, timestep);

  //double deltat = (T_max - T_min) / 100;

  /// Force tb_info switch to 1 for selected cells in serial calculations
  /*
  if (modelgridindex % 14 == 0)
  {
    tb_info = 1;
  }
  else tb_info = 0;
  */
  /// Force tb_info for all cells
  //tb_info = 1;

  /// Check whether the thermal balance equation has a root in [T_min,T_max]
  //gsl_root_fsolver *solver;
  //solver = gsl_root_fsolver_alloc(solvertype);
  //mintemp_f.function = &mintemp_solution_f;
  //maxtemp_f.function = &maxtemp_solution_f;

  Te_solution_paras paras;
  paras.modelgridindex = modelgridindex;
  paras.t_current = t_current;
  paras.heatingcoolingrates = heatingcoolingrates;

  gsl_function find_T_e_f;
  find_T_e_f.function = &T_e_eqn_heating_minus_cooling;
  find_T_e_f.params = &paras;

  double thermalmin = T_e_eqn_heating_minus_cooling(T_min, find_T_e_f.params);
  double thermalmax = T_e_eqn_heating_minus_cooling(T_max, find_T_e_f.params);
  // printout("(heating - cooling) at T_min: %g, at T_max: %g\n",thermalmin,thermalmax);
  if (!isfinite(thermalmin) || !isfinite(thermalmax))
  {
    printout("[abort request] call_T_e_finder: non-finite results in modelcell %d (T_R=%g, W=%g). T_e forced to be MINTEMP\n",
             modelgridindex, get_TR(modelgridindex), get_W(modelgridindex));
    thermalmax = thermalmin = -1;
  }

  //double thermalmax = find_T_e(T_max,find_T_e_f.params);
  double T_e;
  if (thermalmin * thermalmax < 0)
  {
    /// If it has, then solve for the root T_e
    /// but first of all printout heating and cooling rates if the tb_info switch is set
    /*if (tb_info == 1)
    {
      for (i = 0; i <= 100; i++)
      {
        T_e = T_min + i*deltat;
        thermalmin = find_T_e(T_e,find_T_e_f.params);
        //fprintf(tb_file,"%d %g %g %g %g %g %g %g %g %g %g %g\n",modelgridindex,T_e,heatingrates[tid].ff,heatingrates[tid].bf,heatingrates[tid].collbb,heatingrates[tid].collbf,heatingrates[tid].gamma,coolingrates[tid].ff,coolingrates[tid].fb,coolingrates[tid].collbb,coolingrates[tid].collbf,coolingrates[tid].adiabatic);
        fprintf(tb_file,"%d %g %g %g %g %g %g %g %g %g\n",modelgridindex,T_e,heatingrates[tid].ff,heatingrates[tid].bf,heatingrates[tid].collisional, heatingrates[tid].gamma,coolingrates[tid].ff,coolingrates[tid].fb,coolingrates[tid].collisional,coolingrates[tid].adiabatic);
      }
    }*/
    /// now start so iterate on T_e solution
    ///onedimensional gsl root solver, derivative type
    /*const gsl_root_fdfsolver_type *solvertype;
    solvertype = gsl_root_fdfsolver_newton;
    gsl_root_fdfsolver *solver;
    solver = gsl_root_fdfsolver_alloc(solvertype);

    gsl_function_fdf fdf;
    fdf.f = &nne_solution_f;
    fdf.df = &nne_solution_deriv;
    fdf.fdf = &nne_solution_fdf;*/

    // one-dimensional gsl root solver, bracketing type
    const gsl_root_fsolver_type *solvertype = gsl_root_fsolver_brent;
    gsl_root_fsolver *T_e_solver = gsl_root_fsolver_alloc(solvertype);

    gsl_root_fsolver_set(T_e_solver, &find_T_e_f, T_min, T_max);
    const double fractional_accuracy = 1e-3;
    const int maxit = 100;
    int status;
    for (int iternum = 0; iternum < maxit; iternum++)
    {
      gsl_root_fsolver_iterate(T_e_solver);
      T_e = gsl_root_fsolver_root(T_e_solver);
      const double T_e_min = gsl_root_fsolver_x_lower(T_e_solver);
      const double T_e_max = gsl_root_fsolver_x_upper(T_e_solver);
      status = gsl_root_test_interval(T_e_min, T_e_max, 0, fractional_accuracy);
      // printout("iter %d, T_e interval [%g, %g], guess %g, status %d\n", iternum, T_e_min, T_e_max, T_e, status);
      if (status != GSL_CONTINUE)
      {
        printout("after %d iterations, T_e = %g K, interval [%g, %g]\n", iternum + 1, T_e, T_e_min, T_e_max);
        break;
      }
    }

    if (status == GSL_CONTINUE)
      printout("[warning] call_T_e_finder: T_e did not converge within %d iterations\n", maxit);

    gsl_root_fsolver_free(T_e_solver);
  }
  /// Quick solver style: works if we can assume that there is either one or no
  /// solution on [MINTEM.MAXTEMP] (check that by doing a plot of heating-cooling
  /// vs. T_e using the tb_info switch)
  else if (thermalmax < 0)
  {
    /// Thermal balance equation always negative ===> T_e = T_min
    /// Calculate the rates again at this T_e to print them to file
    T_e = MINTEMP;
    //if (nts_global >= 15) T_e = MAXTEMP; ///DEBUG ONLY!!!!!!!!!!!!!!!!!!!!!!!!! invert boundaries!
    printout("[warning] call_T_e_finder: cooling bigger than heating at lower T_e boundary %g in modelcell %d (T_R=%g,W=%g). T_e forced to be MINTEMP\n",MINTEMP,modelgridindex,get_TR(modelgridindex),get_W(modelgridindex));
  }
  else
  {
    /// Thermal balance equation always negative ===> T_e = T_max
    /// Calculate the rates again at this T_e to print them to file
    T_e = MAXTEMP;
    //if (nts_global >= 15) T_e = MINTEMP; ///DEBUG ONLY!!!!!!!!!!!!!!!!!!!!!!!!! invert boundaries!
    printout("[warning] call_T_e_finder: heating bigger than cooling over the whole T_e range [%g,%g] in modelcell %d (T_R=%g,W=%g). T_e forced to be MAXTEMP\n",MINTEMP,MAXTEMP,modelgridindex,get_TR(modelgridindex),get_W(modelgridindex));
    //printout("hot cell %d, with T_R %g, T_e %g, W %g, nne %g\n",cellnumber,cell[cellnumber].T_R,T_e,cell[cellnumber].W,cell[cellnumber].nne);
  }

  /// Otherwise the more expensive solver style _might_ work
  /*
  else
  {
      /// Search for the first temperature point where the cooling rates dominate over
      /// the heating rates. This determines the interval in which the solution must be
      double helper = -99.;
      for (i = 0; i <= 100; i++)
      {
        T_e = T_min + i*deltat;
        thermalmin = find_T_e(T_e,find_T_e_f.params);
        if (tb_info == 1)
        {
          fprintf(tb_file,"%d %g %g %g %g %g %g %g %g %g\n",modelgridindex,T_e,heatingrates[tid].ff,heatingrates[tid].bf,heatingrates[tid].collisional, heatingrates[tid].gamma,coolingrates[tid].ff,coolingrates[tid].fb,coolingrates[tid].collisional,coolingrates[tid].adiabatic);
        }
        if (thermalmin <= 0 && helper < 0)
        {
          helper = T_e;
          if (tb_info != 1) break;
        }
        //if (find_T_e(T_e,find_T_e_f.params) <= 0) break;
      }
      if (helper < 0) helper = MAXTEMP;
      T_min = helper-deltat;
      T_max = helper;
      //T_min = T_e-deltat;
      //T_max = T_e;
      //if (tb_info == 1) printout("T_min %g, T_max %g for cell %d\n",T_min,T_max,modelgridindex);

      if (T_max == MAXTEMP)
      {
        printout("[warning] call_T_e_finder: heating bigger than cooling over the whole T_e range [%g,%g] in cell %d (T_R=%g,W=%g). T_e forced to be MAXTEMP\n",MINTEMP,MAXTEMP,modelgridindex,get_TR(modelgridindex),get_W(modelgridindex));
        T_e = MAXTEMP;
      }
      else if (T_min < MINTEMP)
      {
        printout("[warning] call_T_e_finder: cooling bigger than heating at lower T_e boundary %g in cell %d (T_R=%g,W=%g). T_e forced to be MINTEMP\n",MINTEMP,modelgridindex,modelgridindex,get_TR(modelgridindex),get_W(modelgridindex));
        T_e =  MINTEMP;
      }
      else
      {
        /// We found an interval which has a solution. Solve for the root
        T_e_solver = gsl_root_fsolver_alloc(solvertype);
        gsl_root_fsolver_set(T_e_solver, &find_T_e_f, T_min, T_max);
        iter2 = 0;
        do
        {
          iter2++;
          status = gsl_root_fsolver_iterate(T_e_solver);
          T_e = gsl_root_fsolver_root(T_e_solver);
          //cell[cellnumber].T_e = T_e;
          set_Te(modelgridindex,T_e);
          T_e_min = gsl_root_fsolver_x_lower(T_e_solver);
          T_e_max = gsl_root_fsolver_x_upper(T_e_solver);
          status = gsl_root_test_interval(T_e_min,T_e_max,0,fractional_accuracy);
          //printout("[debug] find T_e:   %d [%g, %g] %g %g\n",iter2,x_lo,x_hi,x_0,x_hi-x_lo);
        }
        while (status == GSL_CONTINUE && iter2 < maxit);
        if (status == GSL_CONTINUE) printout("[warning] call_T_e_finder: T_e did not converge within %d iterations\n",maxit);
        gsl_root_fsolver_free(T_e_solver);
        //printout("%d %g %g %g %g %g %g %g %g %g %g %g %g\n",cellnumber,T_e,ffheatingestimator[cellnumber*nelements*maxion+0*maxion+0],bfheatingestimator[cellnumber*nelements*maxion+0*maxion+0],heatingrates[tid].collbb,heatingrates[tid].collbf,heatingrates[tid].gamma,coolingrates[tid].ff,coolingrates[tid].fb,coolingrates[tid].collbb,coolingrates[tid].collbf,coolingrates[tid].adiabatic,ffheatingestimator[cellnumber*nelements*maxion+0*maxion+0]+bfheatingestimator[cellnumber*nelements*maxion+0*maxion+0]+heatingrates[tid].collbb+heatingrates[tid].collbf+heatingrates[tid].gamma-coolingrates[tid].ff-coolingrates[tid].fb-coolingrates[tid].collbb-coolingrates[tid].collbf-coolingrates[tid].adiabatic);
      }
  }
  */

  //fprintf(heating_file,"%d %g %g %g %g %g %g %g %g\n",modelgridindex,heatingrates[tid].ff,heatingrates[tid].bf,heatingrates[tid].collisional, heatingrates[tid].gamma,coolingrates[tid].ff,coolingrates[tid].fb,coolingrates[tid].collisional,coolingrates[tid].adiabatic);

  if (neutral_flag)
    printout("[info] call_T_e_finder: cell %d contains only neutral ions\n", modelgridindex);

  if (T_e > 2 * T_e_old)
  {
    T_e = 2 * T_e_old;
    printout("use T_e damping in cell %d\n", modelgridindex);
    if (T_e > MAXTEMP)
      T_e = MAXTEMP;
  }
  else if (T_e < 0.5 * T_e_old)
  {
    T_e = 0.5 * T_e_old;
    printout("use T_e damping in cell %d\n", modelgridindex);
    if (T_e < MINTEMP)
      T_e = MINTEMP;
  }

  set_Te(modelgridindex, T_e);
}


// void calculate_oldheating_rates(int cellnumber)
// /// Calculate the heating rates for a given cell. Results are returned
// /// via the elements of the global heatingrates data structure.
// {
//   double interpolate_bfheatingcoeff_below(int element, int ion, int level, double T_R);
//   double interpolate_bfheatingcoeff_above(int element, int ion, int level, double T_R);
//   double bfheating_integrand_gsl(double nu, void *paras);
//   double ffheating_integrand_gsl(double nu, void *paras);
//   double col_deexcitation(PKT *pkt_ptr, int lower, double epsilon_trans, double statweight_target, int lineindex);
//   double col_recombination(PKT *pkt_ptr, int lower, double epsilon_trans);
//   double calculate_exclevelpop(int cellnumber, int element, int ion, int level);
//
//   double epsilon_current,epsilon_target,epsilon_trans,nnlevel;
//   double C;
//   double statweight_target;
//   double T_e,T_R,T_D,W,W_D;
//
//   int nions,nlevels_currention,nlevels_lowerion,ndowntrans;
//   int element,ion,level,lower;
//   int ii,lineindex;
//
//   double intaccuracy = 1e-2;
//   double error;
//   //size_t neval; ///for qng integrator
//
//   /*
//   gsl_function F_bfheating;
//   F_bfheating.function = &bfheating_integrand_gsl;
//   gslintegration_bfheatingparas bfheatingparas;
//   bfheatingparas.cellnumber = cellnumber;
//   double bfhelper;
//   */
//
//   PKT dummypkt;
//   dummypkt.where = cellnumber;
//   PKT *pkt_ptr;
//   pkt_ptr = &dummypkt;
//
//   T_e = cell[cellnumber].T_e; ///WHY EQUALED THAT T_R ??? TYPO OR DEEPER SENSE ????????????????????????????????????????????????????????????????????
//   T_R = cell[cellnumber].T_R;
// //  T_D = cell[cellnumber].T_D;
//   W = cell[cellnumber].W;
// //  W_D = cell[cellnumber].W_D;
//
//   double C_deexc = 0.;
//   double C_recomb = 0.;
//   double bfheating = 0.;
//   double ffheating = 0.;
//
//   nlevels_lowerion = 0;
//   for (element = 0; element < nelements; element++)
//   {
//     mastate[tid].element = element;
//     nions = get_nions(element);
//     for (ion = 0; ion < nions; ion++)
//     {
//       mastate[tid].ion = ion;
//       nlevels_currention = get_nlevels(element,ion);
//       if (ion > 0) nlevels_lowerion = get_nlevels(element,ion-1);
//
// //       for (level = 0; level < nlevels_currention; level++)
// //       {
// //         epsilon_current = epsilon(element,ion,level);
// //         mastate[tid].level = level;
// //         nnlevel = calculate_exclevelpop(cellnumber,element,ion,level);
// //
// //
// //         /// Collisional heating: deexcitation to same ionization stage
// //         /// ----------------------------------------------------------
// //         ndowntrans = get_ndowntrans(element, ion, level);
// //         for (ii = 1; ii <= ndowntrans; ii++)
// //         {
// //           lower = elements[element].ions[ion].levels[level].downtrans[ii].targetlevel;
// //           epsilon_target = elements[element].ions[ion].levels[level].downtrans[ii].epsilon;
// //           lineindex = elements[element].ions[ion].levels[level].downtrans_lineindicies[ii];
// //           epsilon_trans = epsilon_current - epsilon_target;
// //           C = col_deexcitation(pkt_ptr,lower,epsilon_trans,statweight_target,lineindex)*epsilon_trans;
// //           C_deexc += C;
// //         }
// //
// //         /// Collisional heating: recombination to lower ionization stage
// //         /// ------------------------------------------------------------
// //         /// For the moment we deal only with ionisations to the next ions groundlevel.
// //         /// For speed issues (reduced number of calls to epsilon) this is checked here
// //         /// instead of the more general way to check in col_recomb!
// //         if (ion > 0 && level == 0) /// Check whether lower ionisation stage available
// //         {
// //           for (lower = 0; lower < nlevels_lowerion; lower++)
// //           {
// //             epsilon_trans = epsilon_current - epsilon(element,ion-1,lower);
// //             C = col_recombination(pkt_ptr,lower,epsilon_trans)*epsilon_trans;
// //             C_recomb += C;
// //           }
// //         }
// //
// //
// //         /*
// //         /// Bound-free heating (analytical calculation)
// //         /// -------------------------------------------
// //         /// We allow bound free-transitions only if there is a higher ionisation stage
// //         /// left in the model atom to match the bound-free absorption in the rpkt routine.
// //         /// There this condition is needed as we can only ionise to existing ionisation
// //         /// stage even if there would be further ionisation stages in nature which
// //         /// are not included in the model atom.
// //         if (ion < nions-1)
// //         {
// //           epsilon_trans = epsilon(element,ion+1,0) - epsilon_current;
// //           /// NB: W comes from the fact, that the W coming from the radiation field was factored
// //           /// out in the precalculation of the bf-heating coefficient (this is justified by the
// //           /// linear dependence on W).
// //           /// The rate coefficient is calculated under the assumption T_e=T_R because its direct
// //           /// T_e dependence is very weak. This means we have to pass T_R as the temperature
// //           /// even if we are iterating here on T_e. (Otherwise we would allow a large temperature
// //           /// range for T_R which changes the coefficient strongly).
// //           //C = interpolate_bfheatingcoeff(element,ion,level,T_R)*W*nnlevel;
// //           C = nnlevel * (W*interpolate_bfheatingcoeff_below(element,ion,level,T_R));// + W_D*interpolate_bfheatingcoeff_above(element,ion,level,T_D));
// //           if (element == 6)
// //           {
// //             if (ion == 1)
// //             {
// //               printout("T_e %g, nnlevel %g, W %g, heatingcoeff %g\n",T_e,nnlevel,W,interpolate_bfheatingcoeff_below(element,ion,level,T_R));
// //             }
// //           }
// //
// //           /// Exact calculation of the bf-heating coefficients using integrators.
// //           /// This makes things SLOW!!!
// //           //bfheatingparas.nu_edge = epsilon_trans/H;
// //           //F_bfheating.params = &bfheatingparas;
// //           /// Discuss about the upper frequency limit (here 1e16 Hz) which we should choose
// //           //gsl_integration_qag(&F_bfheating, bfheatingparas.nu_edge, 10*bfheatingparas.nu_edge, 0, intaccuracy, 1024, 6, wspace, &bfhelper, &error);
// //           //C = bfhelper * FOURPI * nnlevel;
// //           bfheating += C;
// //         }
// //         */
// //       }
//
//       /// Bound-free heating (from estimators)
//       /// ------------------------------------
//       if (ion < nions-1) bfheating += bfheatingestimator_save[cellnumber*nelements*maxion+element*maxion+ion];
//     }
//   }
//
//
//
//
//
//
//   /// Free-free heating
//   /// -----------------
//   /// From estimators
//   ffheating = ffheatingestimator_save[cellnumber];
//
//   /*
//   /// Analytical calculation using T, and populations
//   /// This is always taken as an additional process to the other importantheatingterms
//   /// because its not done per ion but overall ions.
//
//   gsl_function F_ffheating;
//   F_ffheating.function = &ffheating_integrand_gsl;
//
//   gslintegration_ffheatingparas intparas;
//   intparas.T_e = T_e;
//   intparas.cellnumber = cellnumber;
//   F_ffheating.params = &intparas;
//
//   /// Discuss about the upper frequency limit (here 1e16 Hz) which we should choose
//   gsl_integration_qag(&F_ffheating, 0, 1e16, 0, intaccuracy, 1024, 6, wspace, &ffheating, &error);
//   /// or this integrator
//   //gsl_integration_qng(&F_ffheating, 0, 1e16, 0, intaccuracy, &ffheating, &error, &neval);
//   ffheating *= FOURPI;
//   */
//
//   heatingrates[tid].collisional = colheatingestimator_save[cellnumber]; //C_deexc + C_recomb;
//   //heatingrates[tid].collbb = C_deexc;
//   //heatingrates[tid].collbf = C_recomb;
//   heatingrates[tid].bf = bfheating;
//   heatingrates[tid].ff = ffheating;
//   //printout("ffheating %g, bfheating %g, colheating %g\n",ffheating,bfheating,C_deexc+C_recomb);
// }
