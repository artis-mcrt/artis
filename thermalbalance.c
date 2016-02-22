#include "sn3d.h"
#include <gsl/gsl_integration.h>

#ifndef FORCE_LTE

///****************************************************************************
double call_T_e_finder(int modelgridindex, double t_current, int tb_info, double T_min, double T_max)
{
  double find_T_e(double T_e, void *paras);

  double fractional_accuracy = 1e-2;
  int maxit = 100;
  int iter2,status;
  double T_e,T_e_min,T_e_max;
  double thermalmin,thermalmax;
  double deltat;
  int i;

  //gsl_root_fsolver *solver;
  //solver = gsl_root_fsolver_alloc(solvertype);
  //mintemp_f.function = &mintemp_solution_f;
  //maxtemp_f.function = &maxtemp_solution_f;

  ///onedimensional gsl root solver, derivative type
  /*const gsl_root_fdfsolver_type *solvertype;
  solvertype = gsl_root_fdfsolver_newton;
  gsl_root_fdfsolver *solver;
  solver = gsl_root_fdfsolver_alloc(solvertype);

  gsl_function_fdf fdf;
  fdf.f = &nne_solution_f;
  fdf.df = &nne_solution_deriv;
  fdf.fdf = &nne_solution_fdf;*/


  ///onedimensional gsl root solver, bracketing type
  const gsl_root_fsolver_type *solvertype;
  gsl_root_fsolver *T_e_solver;
  solvertype = gsl_root_fsolver_brent;

  Te_solution_paras paras;
  paras.cellnumber = modelgridindex;
  paras.t_current = t_current;

  gsl_function find_T_e_f;
  find_T_e_f.function = &find_T_e;
  find_T_e_f.params = &paras;

  deltat = (T_max-T_min)/100;

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
  thermalmin = find_T_e(T_min,find_T_e_f.params);
  thermalmax = find_T_e(T_max,find_T_e_f.params);
  if (!isfinite(thermalmin) || !isfinite(thermalmax))
  {
    printout("thermalmin %g, thermalmax %g\n",thermalmin,thermalmax);
    printout("[abort request] call_T_e_finder: non-finte results in modelcell %d (T_R=%g,W=%g). T_e forced to be MINTEMP\n",modelgridindex,get_TR(modelgridindex),get_W(modelgridindex));
    thermalmax = thermalmin = -1;
  }

  printout("thermalmin %g, thermalmax %g\n",thermalmin,thermalmax);
  //double thermalmax = find_T_e(T_max,find_T_e_f.params);
  if (thermalmin*thermalmax < 0)
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
      printout("[debug] find T_e:   iter %d, interval [%g, %g], guess %g, status %d\n",iter2,T_e_min,T_e_max,T_e,status);
    }
    while (status == GSL_CONTINUE && iter2 < maxit);

    if (status == GSL_CONTINUE)
      printout("[warning] call_T_e_finder: T_e did not converge within %d iterations\n",maxit);

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
    //cell[n].T_e = T_e;
  }
  else
  {
    /// Thermal balance equation always negative ===> T_e = T_max
    /// Calculate the rates again at this T_e to print them to file
    T_e = MAXTEMP;
    //if (nts_global >= 15) T_e = MINTEMP; ///DEBUG ONLY!!!!!!!!!!!!!!!!!!!!!!!!! invert boundaries!
    printout("[warning] call_T_e_finder: heating bigger than cooling over the whole T_e range [%g,%g] in modelcell %d (T_R=%g,W=%g). T_e forced to be MAXTEMP\n",MINTEMP,MAXTEMP,modelgridindex,get_TR(modelgridindex),get_W(modelgridindex));
    //printout("hot cell %d, with T_R %g, T_e %g, W %g, nne %g\n",cellnumber,cell[cellnumber].T_R,T_e,cell[cellnumber].W,cell[cellnumber].nne);
    //cell[n].T_e = T_e;
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

  if (neutral_flag == 1) printout("[info] call_T_e_finder: cell %d contains only neutral ions\n",modelgridindex);
  return T_e;
}



///***************************************************************************/
double find_T_e(double T_e, void *paras)
/// Thermal balance equation on which we have to iterate to get T_e
{
  double calculate_populations(int modelgridindex, int first_nonempty_cell);
  void calculate_cooling_rates(int modelgridindex);
  void calculate_heating_rates(int modelgridindex);
  double nntot;
  int element,ion,nions;

  int modelgridindex = ((Te_solution_paras *) paras)->cellnumber;
  double t_current = ((Te_solution_paras *) paras)->t_current;

  /// Set new T_e guess for the current cell and update populations
  //cell[cellnumber].T_e = T_e;
  set_Te(modelgridindex,T_e);
  nntot = calculate_populations(modelgridindex,0);

  /// Then calculate heating and cooling rates
  calculate_cooling_rates(modelgridindex);
  calculate_heating_rates(modelgridindex);
  /// These heating rates using estimators work only for hydrogen!!!
  //double ionstagepop(int cellnumber, int element, int ion);
  //double heating_ff,heating_bf;
  //double nne = cell[cellnumber].nne;
  //double W = cell[cellnumber].W;
  //heating_ff = cell[cellnumber].heating_ff * ionstagepop(cellnumber,0,1)/cell[cellnumber].composition[0].partfunct[1]*nne / sqrt(T_e);
  //heating_bf = cell[cellnumber].heating_bf * W*ionstagepop(cellnumber,0,0)/cell[cellnumber].composition[0].partfunct[0];

  /// If selected take direct gamma heating into account
  if (do_rlc_est == 3)
  {
    heatingrates[tid].gamma = rpkt_emiss[modelgridindex] * 1.e20; ///1.e20 since the emissivities are all scaled by this
    heatingrates[tid].gamma *= 4*PI; /// This was missing here! rpkt_emiss is normalised to give a real emissivity
                                     /// for the formal integral calculation!!!
    // Above is the gamma-ray bit. Below is *supposed* to be the kinetic energy of positrons created by 56Co and 48V. These formulae should be checked, however.
    heatingrates[tid].gamma += (0.610*0.19*MEV)*(exp(-1.*time_step[nts_global].mid/TCOBALT) - exp(-1.*time_step[nts_global].mid/TNICKEL))/(TCOBALT-TNICKEL)*modelgrid[modelgridindex].fni*get_rho(modelgridindex)/MNI56;
    heatingrates[tid].gamma += (0.290*0.499*MEV)*(exp(-1.*time_step[nts_global].mid/T48V) - exp(-1.*time_step[nts_global].mid/T48CR))/(T48V-T48CR)*modelgrid[modelgridindex].f48cr*get_rho(modelgridindex)/MCR48;
  }
  else
  {
    heatingrates[tid].gamma = 0.;
  }

  //heatingrates[tid].gamma = cell[cellnumber].f_ni*cell[cellnumber].rho_init/MNI56 * pow(tmin/t_current,3) * (ENICKEL/TNICKEL*exp(-t_current/TNICKEL) + ECOBALT/(TCOBALT-TNICKEL)*(exp(-t_current/TCOBALT)-exp(-t_current/TNICKEL)));
  //heatingrates[tid].gamma *= 0.01;
  //double factor = -1./(t_current*(-TCOBALT+TNICKEL));
  //factor *= (-ENICKEL*exp(-t_current/TNICKEL)*t_current*TCOBALT - ENICKEL*exp(-t_current/TNICKEL)*TNICKEL*TCOBALT + ENICKEL*exp(-t_current/TNICKEL)*t_current*TNICKEL + pow(TNICKEL,2)*ENICKEL*exp(-t_current/TNICKEL) - TCOBALT*t_current*ECOBALT*exp(-t_current/TCOBALT) - pow(TCOBALT,2)*ECOBALT*exp(-t_current/TCOBALT) + ECOBALT*t_current*TNICKEL*exp(-t_current/TNICKEL) + pow(TNICKEL,2)*ECOBALT*exp(-t_current/TNICKEL) + ENICKEL*TCOBALT*TNICKEL - ENICKEL*pow(TNICKEL,2) - pow(TNICKEL,2)*ECOBALT + ECOBALT*pow(TCOBALT,2));
  //heatingrates[tid].gamma = cell[cellnumber].f_ni*cell[cellnumber].rho_init/MNI56 * pow(tmin/t_current,3) * ((ENICKEL/TNICKEL*exp(-t_current/TNICKEL) + ECOBALT/(TCOBALT-TNICKEL)*(exp(-t_current/TCOBALT)-exp(-t_current/TNICKEL))) - factor/t_current);

  /// Adiabatic cooling term
  double p,dV,V;
  p = nntot*KB*T_e;
  dV = 3*pow(wid_init/tmin,3)*pow(t_current,2);
  V = pow(wid_init*t_current/tmin,3);
  //printout("nntot %g, p %g, dV %g, V %g\n",nntot,p,dV,V);
  coolingrates[tid].adiabatic = p*dV/V;

  return heatingrates[tid].ff + heatingrates[tid].bf + heatingrates[tid].collisional - coolingrates[tid].ff - coolingrates[tid].fb - coolingrates[tid].collisional + heatingrates[tid].gamma - coolingrates[tid].adiabatic; // - 0.01*(heatingrates[tid].bf+coolingrates[tid].fb)/2;
  //return heatingrates[tid].gamma - coolingrates[tid].fb; // - 0.01*(heatingrates[tid].bf+coolingrates[tid].fb)/2;
  //return heatingrates[tid].ff + heatingrates[tid].bf + heatingrates[tid].collbf - coolingrates[tid].ff - coolingrates[tid].fb - coolingrates[tid].collbf + heatingrates[tid].gamma - coolingrates[tid].adiabatic - 0.01*(heatingrates[tid].bf+coolingrates[tid].fb)/2;
}



///****************************************************************************
void calculate_heating_rates(int modelgridindex)
/// Calculate the heating rates for a given cell. Results are returned
/// via the elements of the global heatingrates data structure.
{
  double get_bfheatingcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
  double interpolate_bfheatingcoeff_below(int element, int ion, int level, double T_R);
  double interpolate_bfheatingcoeff_above(int element, int ion, int level, double T_R);
  double bfheating_integrand_gsl(double nu, void *paras);
  double ffheating_integrand_gsl(double nu, void *paras);
  double col_deexcitation(int modelgridindex, int lower, double epsilon_trans, double statweight_target, int lineindex);
//  double col_deexcitation(PKT *pkt_ptr, int lower, double epsilon_trans, double statweight_target, int lineindex);
//  double col_recombination(PKT *pkt_ptr, int lower, double epsilon_trans);
  double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);

  double epsilon_current,epsilon_target,epsilon_trans,nnlevel;
  double C;
  double statweight_target;
  double T_e,T_R,T_D,W,W_D;

  int nions,nlevels_currention,nlevels_lowerion,ndowntrans,phixstargetindex;
  int element,ion,level,lower;
  int ii,lineindex;

  gsl_integration_workspace *wspace = gsl_integration_workspace_alloc(1000);
  double intaccuracy = 1e-2;
  double error;
  //size_t neval; ///for qng integrator

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

  T_e = get_Te(modelgridindex);
  T_R = get_TR(modelgridindex);
  W = get_W(modelgridindex);

  double C_deexc = 0.;
  double C_recomb = 0.;
  double bfheating = 0.;
  double ffheating = 0.;

  nlevels_lowerion = 0;
  for (element = 0; element < nelements; element++)
  {
    mastate[tid].element = element;
    nions = get_nions(element);
    for (ion = 0; ion < nions; ion++)
    {
      #ifdef DIRECT_COL_HEAT
      {
        mastate[tid].ion = ion;
        nlevels_currention = get_nlevels(element,ion);
//      if (ion > 0) nlevels_lowerion = get_nlevels(element,ion-1);

       for (level = 0; level < nlevels_currention; level++)
       {
         epsilon_current = epsilon(element,ion,level);
         mastate[tid].level = level;
         nnlevel = calculate_exclevelpop(modelgridindex,element,ion,level);
         mastate[tid].nnlevel = nnlevel;
         mastate[tid].statweight = stat_weight(element,ion,level);
//
//
//         /// Collisional heating: deexcitation to same ionization stage
//         /// ----------------------------------------------------------
         ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
         for (ii = 1; ii <= ndowntrans; ii++)
         {
           lower = elements[element].ions[ion].levels[level].downtrans[ii].targetlevel;
           epsilon_target = elements[element].ions[ion].levels[level].downtrans[ii].epsilon;
           statweight_target = elements[element].ions[ion].levels[level].downtrans[ii].stat_weight;
           lineindex = elements[element].ions[ion].levels[level].downtrans[ii].lineindex;
           epsilon_trans = epsilon_current - epsilon_target;
           C = col_deexcitation(modelgridindex,lower,epsilon_trans,statweight_target,lineindex) * epsilon_trans;
           C_deexc += C;
         }
       }
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
//           //gsl_integration_qag(&F_bfheating, bfheatingparas.nu_edge, 10*bfheatingparas.nu_edge, 0, intaccuracy, 1000, 6, wspace, &bfhelper, &error);
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
      if (ion < nions-1)
      {
        //nlevels_currention = get_nlevels(element,ion);
        //nlevels_currention = get_ionisinglevels(element,ion);
        nlevels_currention = get_bfcontinua(element,ion);
        for (level = 0; level < nlevels_currention; level++)
        {
          nnlevel = calculate_exclevelpop(modelgridindex,element,ion,level);
          for (phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            bfheating += nnlevel * get_bfheatingcoeff(element,ion,level,phixstargetindex,modelgridindex);
          }
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
  gsl_integration_qag(&F_ffheating, 0, 1e16, 0, intaccuracy, 1000, 6, wspace, &ffheating, &error);
  /// or this integrator
  //gsl_integration_qng(&F_ffheating, 0, 1e16, 0, intaccuracy, &ffheating, &error, &neval);
  ffheating *= FOURPI;
  */


  /// Collisional heating (from estimators)
  /// -------------------------------------
  heatingrates[tid].collisional = colheatingestimator[modelgridindex];//C_deexc + C_recomb;
#ifdef DIRECT_COL_HEAT
  heatingrates[tid].collisional = C_deexc;
#endif

//  heatingrates[tid].collbb = C_deexc;
//  heatingrates[tid].collbf = C_recomb;
  heatingrates[tid].bf = bfheating;
  heatingrates[tid].ff = ffheating;
  //printout("ffheating %g, bfheating %g, colheating %g\n",ffheating,bfheating,C_deexc+C_recomb);

  gsl_integration_workspace_free(wspace);
}



///****************************************************************************
void calculate_cooling_rates(int modelgridindex)
/// Calculate the cooling rates for a given cell. Results are returned
/// via the elements of the global coolingrates data structure.
{
  double get_bfcooling(int element, int ion, int level, int phixstargetindex, int modelgridindex);
  double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
  int get_ionstage(int element, int ion);
  double ionstagepop(int modelgridindex, int element, int ion);
  double get_groundlevelpop(int modelgridindex, int element, int ion);

  double C,C_ff,C_fb,C_exc,C_ion;
  double nncurrention,nnnextionlevel;
  double epsilon_current,epsilon_upper,epsilon_trans;
  int element,ion,level,upper,phixstargetindex;
  int nlevels_currention,ionisinglevels;
  int ii,lineindex;
  int ioncharge;
  int nions,nuptrans;

  double nne = get_nne(modelgridindex);
  double T_e = get_Te(modelgridindex);

/*  PKT dummypkt;
  dummypkt.where = cellnumber;
  PKT *pkt_ptr;
  pkt_ptr = &dummypkt;*/

  /// calculate rates for
  C = 0.;
  C_ff = 0.;   /// free-free creation of rpkts
  C_fb = 0.;   /// free-bound creation of rpkt
  C_exc = 0.;  /// collisional excitation of macroatoms
  C_ion = 0.;  /// collisional ionisation of macroatoms
  for (element=0; element < nelements; element++)
  {
    //printout("[debug] do_kpkt: element %d\n",element);
    mastate[tid].element = element;
    nions = get_nions(element);
    for (ion=0; ion < nions; ion++)
    {
      //printout("[debug] do_kpkt: ion %d\n",ion);
      mastate[tid].ion = ion;
      nlevels_currention = get_nlevels(element,ion);
      //ionisinglevels = get_ionisinglevels(element,ion);
      ionisinglevels = get_bfcontinua(element,ion);
      nnnextionlevel = get_groundlevelpop(modelgridindex,element,ion+1);
      nncurrention = ionstagepop(modelgridindex,element,ion);

      /// ff creation of rpkt
      /// -------------------
      ioncharge = get_ionstage(element,ion) - 1;
      //printout("[debug] ioncharge %d, nncurrention %g, nne %g\n",ion,nncurrention,nne);
      if (ioncharge > 0)
      {
        C = 1.426e-27 * sqrt(T_e) * pow(ioncharge,2) * nncurrention * nne;
        C_ff += C;
      }

      for (level = 0; level < nlevels_currention; level++)
      {
        //printout("[debug] do_kpkt: element %d, ion %d, level %d\n",element,ion,level);
        epsilon_current = epsilon(element,ion,level);
        mastate[tid].level = level;
        mastate[tid].nnlevel = calculate_exclevelpop(modelgridindex,element,ion,level);

        /// excitation to same ionization stage
        /// -----------------------------------
        nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
        for (ii = 1; ii <= nuptrans; ii++)
        {
          upper = elements[element].ions[ion].levels[level].uptrans[ii].targetlevel;
          lineindex = elements[element].ions[ion].levels[level].uptrans[ii].lineindex;
          //printout("    excitation to level %d possible\n",upper);
          epsilon_trans = epsilon(element,ion,upper) - epsilon_current;
          C = col_excitation(modelgridindex,upper,lineindex,epsilon_trans) * epsilon_trans;
          C_exc += C;
        }

        if (ion < nions-1 && level < ionisinglevels) ///check whether further ionisation stage available
        {
          //printout("    ionisation possible\n");
          /// ionization to higher ionization stage
          /// -------------------------------------
          C = 0.0;
          for (phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            epsilon_upper = epsilon(element,ion+1,upper);
            epsilon_trans = epsilon_upper - epsilon_current;
            //printout("cooling list: col_ionization\n");
            C += col_ionization(modelgridindex,phixstargetindex,epsilon_trans) * epsilon_trans;
          }
          C_ion += C;

          /// fb creation of r-pkt
          /// free bound rates are calculated from the lower ion, but associated to the higher ion
          /// --------------------
          //upper = 0;
          //epsilon_upper = epsilon(element,ion+1,0);
          //E_threshold = epsilon_upper - epsilon_current;
          //E_threshold = epsilon_trans;
          //printout("get_bfcooling(%d,%d,%d,%d) for histindex %d\n",element,ion,level,cellnumber,histindex);
          C = 0.0;
          for (phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            C += get_bfcooling(element,ion,level,phixstargetindex,modelgridindex);
          }
          C_fb += C;
        }
      }
    }
  }

  coolingrates[tid].collisional = C_exc + C_ion;
  coolingrates[tid].collbb = C_exc;
  coolingrates[tid].collbf = C_ion;
  coolingrates[tid].fb = C_fb;
  coolingrates[tid].ff = C_ff;
}





///****************************************************************************
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
//   double stat_weight(int element, int ion, int level);
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
//   gsl_integration_workspace *wspace = gsl_integration_workspace_alloc (1000);
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
// //         mastate[tid].nnlevel = nnlevel;
// //         mastate[tid].statweight = stat_weight(element,ion,level);
// //
// //
// //         /// Collisional heating: deexcitation to same ionization stage
// //         /// ----------------------------------------------------------
// //         ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
// //         for (ii = 1; ii <= ndowntrans; ii++)
// //         {
// //           lower = elements[element].ions[ion].levels[level].downtrans[ii].targetlevel;
// //           epsilon_target = elements[element].ions[ion].levels[level].downtrans[ii].epsilon;
// //           statweight_target = elements[element].ions[ion].levels[level].downtrans[ii].stat_weight;
// //           lineindex = elements[element].ions[ion].levels[level].downtrans[ii].lineindex;
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
// //           //gsl_integration_qag(&F_bfheating, bfheatingparas.nu_edge, 10*bfheatingparas.nu_edge, 0, intaccuracy, 1000, 6, wspace, &bfhelper, &error);
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
//   gsl_integration_qag(&F_ffheating, 0, 1e16, 0, intaccuracy, 1000, 6, wspace, &ffheating, &error);
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
//
//   gsl_integration_workspace_free(wspace);
// }

#endif
