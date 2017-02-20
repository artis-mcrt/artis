#include "sn3d.h"
#include <gsl/gsl_integration.h>


///****************************************************************************
void tabulate_ratecoefficients_gsl()
/// Precalculates the rate coefficients for stimulated and spontaneous
/// recombination and photoionisation on a given temperature grid using
/// libgsl integrators.
/// NB: with the nebular approximation they only depend on T_e, T_R and W. 
/// W is easily factored out. For stimulated recombination we must assume
/// T_e = T_R for this precalculation.
{
  double calculate_sahafact(int element, int ion, int level, double T, double E_threshold);
  double interpolate_spontrecombcoeff(int element, int ion, int level, double T);
  double epsilon(int element, int ion, int level);
  int get_ionstage(int element, int ion);
  int get_element(int element);

  
  double alpha_sp_integrand_gsl(double nu, void *paras);
  double alpha_sp_E_integrand_gsl(double nu, void *paras);
  double gamma_integrand_gsl(double nu, void *paras);
  double gammacorr_integrand_gsl(double nu, void *paras);
  double approx_bfheating_integrand_gsl(double nu, void *paras);
  double bfcooling_integrand_gsl(double nu, void *paras);
  //double stimulated_bfcooling_integrand_gsl(double nu, void *paras);
  //double stimulated_recomb_integrand_gsl(double nu, void *paras);
  
  gsl_function F_alpha_sp, F_alpha_sp_E, F_gamma, F_gammacorr;
  gsl_function F_bfheating,F_bfcooling;//,F_stimulated_bfcooling,F_stimulated_recomb;
  F_alpha_sp.function = &alpha_sp_integrand_gsl;
  F_alpha_sp_E.function = &alpha_sp_E_integrand_gsl;
  F_gamma.function = &gamma_integrand_gsl;
  F_gammacorr.function = &gammacorr_integrand_gsl;
  F_bfheating.function = &approx_bfheating_integrand_gsl;
  F_bfcooling.function = &bfcooling_integrand_gsl;
  //F_stimulated_bfcooling.function = &stimulated_bfcooling_integrand_gsl;
  //F_stimulated_recomb.function = &stimulated_recomb_integrand_gsl;
  
  gsl_integration_workspace *w;
  gslintegration_paras intparas;
  FILE *ratecoeff_file;
  double intaccuracy = 1e-2;        /// Fractional accuracy of the integrator
  double E_threshold,nu_threshold;
  double sf,T_e;//,T_R;
  double alpha_sp,gammacorr,bfheating_coeff,bfcooling_coeff;
  double alpha_sp_E,gamma_below,gamma_above,gammacorr_above;
  double bfheating_coeff_above;//,stimulated_bfcooling_coeff,stimulated_recomb_coeff;
  double error;
  double zeta;
  //float ddum;
  float T_min,T_max;
  //size_t neval;   /// for qng integrator
  int nions,nlevels,ionisinglevels,nbfcont;
  int element,ion,level,iter;
  int check,calculate;
  int dum1,dum2,dum3,dum4;
  
  /// Determine the temperture grids gridsize
  T_step = (1.*MAXTEMP-MINTEMP) / (TABLESIZE-1.);               /// global variables
  T_step_log = (log(MAXTEMP)-log(MINTEMP)) / (TABLESIZE-1.);
  
  
  /// Then readin the precalculated rate coefficients from file
  if ((ratecoeff_file = fopen("ratecoeff.dat", "r")) != NULL)
  {
    /// Check whether current temperature range and atomic data match
    /// the precalculated rate coefficients
    fscanf(ratecoeff_file,"%g %g %d",&T_min,&T_max,&dum1);
    printout("Tmin %g Tmax %g TABLESIZE %d \n",T_min,T_max,dum1);
    if (T_min == MINTEMP && T_max == MAXTEMP && dum1 == TABLESIZE)
    {
      check = 1;
      for (element = 0; element < nelements; element++)
      {
        nions = get_nions(element);
        for (ion = 0; ion < nions; ion++)
        {
          fscanf(ratecoeff_file,"%d %d %d %d",&dum1,&dum2,&dum3,&dum4);
          nlevels = get_nlevels(element,ion);
          ionisinglevels = get_ionisinglevels(element,ion);
          if (get_element(element) == dum1 && get_ionstage(element,ion) == dum2 && nlevels == dum3 && ionisinglevels == dum4) 
            check = 1;
          else
          {
              check += 1;
          }
        }
      }
    }
    else check = 100;
  
    if (check == 1)
    {
      printout("[info] tabulate_ratecoefficients_gsl:  Matching ratecoeff.dat file found. Readin this file ...\n");
      for (element = 0; element < nelements; element++)
      {
        nions = get_nions(element) - 1;
        for (ion = 0; ion < nions; ion++)
        {
          //nlevels = get_nlevels(element,ion);
          nlevels = get_ionisinglevels(element,ion); /// number of ionising levels associated with current ion
          nbfcont = get_bfcontinua(element,ion);     /// number of ionising levels of the current ion which are used in the simulation
          for (level = 0; level < nlevels; level++)
          {
            /// Loop over the temperature grid
            for (iter = 0; iter < TABLESIZE; iter++)
            {
              fscanf(ratecoeff_file,"%lg %lg %lg %lg\n", &alpha_sp,&bfcooling_coeff,&gammacorr,&bfheating_coeff);
    
              if (level < nbfcont)
              {
                elements[element].ions[ion].levels[level].spontrecombcoeff[iter] = alpha_sp;
                elements[element].ions[ion].levels[level].bfcooling_coeff[iter] = bfcooling_coeff;
                elements[element].ions[ion].levels[level].corrphotoioncoeff[iter] = gammacorr;
                elements[element].ions[ion].levels[level].bfheating_coeff[iter] = bfheating_coeff;
              }
            }
          }
        }
      }
      fclose(ratecoeff_file);
      calculate = 0;
    }
    else
    {
      printout("[info] tabulate_ratecoefficients_gsl:  No matching ratecoeff.dat file found. Calculate ...\n");
      fclose(ratecoeff_file);
      calculate = 1;
    }
  }
  else
  {
    printout("[info] tabulate_ratecoefficients_gsl:  No ratecoeff.dat file available. Create new one ...\n");
    calculate = 1;
  }
  
  /// Check if we need to calculate the ratecoefficients or if we were able to read them from file
  if (calculate == 1)
  {
    /// Calculate the rate coefficients for each level of each ion of each element
    for (element = 0; element < nelements; element++)
    {
      nions = get_nions(element) - 1;
      for (ion = 0; ion < nions; ion++)
      {
        //nlevels = get_nlevels(element,ion);
        nlevels = get_ionisinglevels(element,ion);
        /// That's only an option for pure LTE
        //if (TAKE_N_BFCONTINUA < nlevels) nlevels = TAKE_N_BFCONTINUA;
        w = gsl_integration_workspace_alloc(1000);
        mastate[tid].element = element;   /// Global variable which passes the current element to all subfunctions of macroatom.c
        mastate[tid].ion = ion;   /// Global variable which passes the current ion to all subfunctions of macroatom.c
        for (level = 0; level < nlevels; level++)
        {
          //printout("element %d, ion %d, level %d, epsilon %g, continuum %g, nlevels %d\n",element,ion,level,epsilon(element,ion,level),epsilon(element,ion+1,0),nlevels);
          alpha_sp = 0.;
          gammacorr = 0.;
          bfheating_coeff = 0.;
          bfcooling_coeff = 0.;

          mastate[tid].level = level;   /// Global variable which passes the current level to all subfunctions of macroatom.c
          E_threshold = epsilon(element,ion+1,0) - epsilon(element,ion,level);
          nu_threshold = E_threshold/H;
          intparas.nu_edge = nu_threshold;   /// Global variable which passes the threshold to the integrator
          
          /// Loop over the temperature grid
          for (iter = 0; iter < TABLESIZE; iter++)
          {
            T_e = MINTEMP * exp(iter*T_step_log);
            //T_e = MINTEMP + iter*T_step;
            sf = calculate_sahafact(element,ion,level,T_e,E_threshold);
            //printout("%d %g\n",iter,T_e);
            
            intparas.T = T_e;
            F_alpha_sp.params = &intparas;
            //F_alpha_sp_E.params = &intparas;
            //F_gamma.params = &intparas;
            F_gammacorr.params = &intparas;
            F_bfcooling.params = &intparas;
            //F_stimulated_bfcooling.params = &intparas;
            //F_stimulated_recomb.params = &intparas;
            F_bfheating.params = &intparas;
            
            /// Spontaneous recombination and bf-cooling coefficient don't depend on the cutted radiation field
            gsl_integration_qag(&F_alpha_sp, nu_threshold, 10*nu_threshold, 0, intaccuracy, 1000, 6, w, &alpha_sp, &error);
            //gsl_integration_qng(&F_alpha_sp, nu_threshold, 10*nu_threshold, 0, intaccuracy, &alpha_sp, &error, &neval);
            alpha_sp *= FOURPI * sf;
            gsl_integration_qag(&F_bfcooling, nu_threshold, 10*nu_threshold, 0, intaccuracy, 1000, 6, w, &bfcooling_coeff, &error);
            bfcooling_coeff *= FOURPI * sf; 

            gsl_integration_qag(&F_gammacorr, nu_threshold, 10*nu_threshold, 0, intaccuracy, 1000, 6, w, &gammacorr, &error);
            gammacorr *= FOURPI; 
            
            gsl_integration_qag(&F_bfheating, nu_threshold, 10*nu_threshold, 0, intaccuracy, 1000, 6, w, &bfheating_coeff, &error);
            bfheating_coeff *= FOURPI; 
            
            /// Save the previously calculated coefficients to memory
            elements[element].ions[ion].levels[level].spontrecombcoeff[iter] = alpha_sp;
            elements[element].ions[ion].levels[level].bfcooling_coeff[iter] = bfcooling_coeff;
            elements[element].ions[ion].levels[level].corrphotoioncoeff[iter] = gammacorr;
            elements[element].ions[ion].levels[level].bfheating_coeff[iter] = bfheating_coeff;
          }
        }
        gsl_integration_workspace_free(w);
      }
    }
  
    /// And the master process writes them to file in a serial operation
    if (rank_global == 0)
    {
      if ((ratecoeff_file = fopen("ratecoeff.dat", "w")) == NULL)
      {
        printout("Cannot open ratecoeff.dat\n");
        exit(0);
      }
      fprintf(ratecoeff_file,"%g %g %d\n",MINTEMP,MAXTEMP,TABLESIZE);
      for (element = 0; element < nelements; element++)
      {
        nions = get_nions(element);
        for (ion = 0; ion < nions; ion++)
        {
          fprintf(ratecoeff_file,"%d %d %d %d\n",get_element(element), get_ionstage(element,ion), get_nlevels(element,ion),  get_ionisinglevels(element,ion));
        }
      }
      
      for (element = 0; element < nelements; element++)
      {
        nions = get_nions(element) - 1;
        for (ion = 0; ion < nions; ion++)
        {
          //nlevels = get_nlevels(element,ion);
          nlevels = get_ionisinglevels(element,ion);
          for (level = 0; level < nlevels; level++)
          {
            /// Loop over the temperature grid
            for (iter = 0; iter < TABLESIZE; iter++)
            {
              alpha_sp = elements[element].ions[ion].levels[level].spontrecombcoeff[iter];
              bfcooling_coeff = elements[element].ions[ion].levels[level].bfcooling_coeff[iter];
              gammacorr = elements[element].ions[ion].levels[level].corrphotoioncoeff[iter];
              bfheating_coeff = elements[element].ions[ion].levels[level].bfheating_coeff[iter];
              //fprintf(ratecoeff_file,"%g %g %g %g %g %g %g %g %g\n", alpha_sp,alpha_sp_E,bfcooling_coeff,gamma_below,gamma_above,gammacorr_below,gammacorr_above,bfheating_coeff_below,bfheating_coeff_above);
              fprintf(ratecoeff_file,"%g %g %g %g\n", alpha_sp,bfcooling_coeff,gammacorr,bfheating_coeff);
            }
          }
        }
      }
      fclose(ratecoeff_file);
    }
  }
  
  
  for (iter = 0; iter < TABLESIZE; iter++)
  {
    T_e = MINTEMP * exp(iter*T_step_log);
    //T_e = MINTEMP + iter*T_step;
    for (element = 0; element < nelements; element++)
    {
      nions = get_nions(element) - 1;
      for (ion = 0; ion < nions; ion++)
      {
        //nlevels = get_nlevels(element,ion); 
        //nlevels = get_ionisinglevels(element,ion); ///number of levels of the current ion which have an associated photoion cross section
        nlevels = get_bfcontinua(element,ion); /// number of ionising levels used in the simulation
        zeta = 0.;
        for (level = 0; level < nlevels; level++)
        {
          zeta += interpolate_spontrecombcoeff(element,ion,level,T_e);
        }
        elements[element].ions[ion].Alpha_sp[iter] = zeta;
//         zeta = interpolate_spontrecombcoeff(element,ion,0,T_e) / zeta;
//         elements[element].ions[ion].zeta[iter] = zeta;
      }
    }
  }
}



///****************************************************************************
/// The following functions define the integrands for these rate coefficients
/// for use with libgsl integrators.
double alpha_sp_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the rate coefficient for spontaneous recombination
/// using gsl integrators.
{
  double photoionization_crosssection(double nu_edge, double nu);
  //double epsilon(int element, int ion, int level);
  double x,sigma_bf;
  
  //int element = mastate[tid].element;
  //int ion = mastate[tid].ion;
  //int level = mastate[tid].level;
  
  double T = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;
  
  //double nu_edge = (epsilon(element,ion+1,0)-epsilon(element,ion,level))/H;
  //printout("[debug] alpha_sp_integrand: element, ion, level: %d, %d, %d\n",exchangepkt_ptr->MA_element,exchangepkt_ptr->MA_ion,exchangepkt_ptr->MA_level);
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  sigma_bf = photoionization_crosssection(nu_edge,nu);
  x = TWOOVERCLIGHTSQUARED * sigma_bf * pow(nu,2) * exp(-HOVERKB*nu/T);
  ///in formula this looks like
  ///x = sigma_bf/H/nu * 2*H*pow(nu,3)/pow(CLIGHT,2) * exp(-H*nu/KB/T);
  
  ///set contributions from Lyman continuum artificially to zero to overcome it's large opacity
  //if (exchangepkt_ptr->MA_level == 0) x = 0;
  return x;
}

double alpha_sp_E_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the rate coefficient for spontaneous recombination
/// using gsl integrators.
{
  double photoionization_crosssection(double nu_edge, double nu);
  //double epsilon(int element, int ion, int level);
  double x,sigma_bf;
  
  //int element = mastate[tid].element;
  //int ion = mastate[tid].ion;
  //int level = mastate[tid].level;
  
  double T = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;
  
  //double nu_edge = (epsilon(element,ion+1,0)-epsilon(element,ion,level))/H;
  //printout("[debug] alpha_sp_integrand: element, ion, level: %d, %d, %d\n",exchangepkt_ptr->MA_element,exchangepkt_ptr->MA_ion,exchangepkt_ptr->MA_level);
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  sigma_bf = photoionization_crosssection(nu_edge,nu);
  x = TWOOVERCLIGHTSQUARED * sigma_bf * pow(nu,3)/nu_edge * exp(-HOVERKB*nu/T);
  ///in formula this looks like
  ///x = sigma_bf/H/nu * 2*H*pow(nu,3)/pow(CLIGHT,2) * exp(-H*nu/KB/T);
  
  ///set contributions from Lyman continuum artificially to zero to overcome it's large opacity
  //if (exchangepkt_ptr->MA_level == 0) x = 0;
  return x;
}


double gamma_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the rate coefficient for photoionization 
/// using gsl integrators.
{
  double photoionization_crosssection(double nu_edge, double nu);
  double radfield2(double nu, double T, double W);
  double T = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;
  double x,sigma_bf;
  
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  sigma_bf = photoionization_crosssection(nu_edge,nu);
  
  /// Dependence on dilution factor W is linear. This allows to set it here to
  /// 1. and scale to its actual value later on.
  x = sigma_bf/H/nu * radfield2(nu,T,1.);
  //x = sigma_bf/H/nu * radfield2(nu,T,1.);
  //if (HOVERKB*nu/T < 1e-2) x = sigma_bf * pow(nu,2)/(HOVERKB*nu/T);
  //else if (HOVERKB*nu/T >= 1e2) x = sigma_bf * pow(nu,2)*exp(-HOVERKB*nu/T);
  //else x = sigma_bf * pow(nu,2)/(exp(HOVERKB*nu/T)-1);
  

  return x;
}

double gammacorr_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the rate coefficient for photoionization 
/// using gsl integrators.
{
  double photoionization_crosssection(double nu_edge, double nu);
  double radfield2(double nu, double T, double W);
  double T = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;
  double x,sigma_bf;
  
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  sigma_bf = photoionization_crosssection(nu_edge,nu);
  
  /// Dependence on dilution factor W is linear. This allows to set it here to
  /// 1. and scale to its actual value later on.
  /// Assumption T_e = T_R makes n_kappa/n_i * (n_i/n_kappa)* = 1
  x = sigma_bf/H/nu * radfield2(nu,T,1.) * (1-exp(-H*nu/KB/T));
  return x;
}

double approx_bfheating_integrand_gsl(double nu, void *paras)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximative way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  double calculate_sahafact(int element, int ion, int level, double T, double E_threshold);
  double photoionization_crosssection(double nu_edge, double nu);
  double radfield2(double nu, double T, double W);
  double x;
  
  
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;
  double sigma_bf = photoionization_crosssection(nu_edge,nu);
  
  /// Precalculation for T_e=T_R and W=1
  float T  = ((gslintegration_paras *) paras)->T;
  x = sigma_bf*(1-nu_edge/nu)*radfield2(nu,T,1) * (1-exp(-H*nu/KB/T));
  
  /// Precalculation for a (T_R,T_e)-grid, but still W is assumed to be 1. 
  /// The radfield part can be corrected later because of its linear dependence. 
  /// But not the W in the stimulated correction term!
  /*float T_e  = ((gslintegration_paras *) paras)->T;
  float T_R  = ((gslintegration_paras *) paras)->T2;
  int element = mastate[tid].element;
  int ion  = mastate[tid].ion;
  int level = mastate[tid].level;
  double E_threshold = nu_edge*H;
  double sf_Te = calculate_sahafact(element,ion,level,T_e,E_threshold);
  double sf_TR = calculate_sahafact(element,ion,level,T_R,E_threshold);
  x = sigma_bf*(1-nu_edge/nu)*radfield2(nu,T_R,1) * (1 - sqrt(T_e/T_R) * sf_Te/sf_TR * exp(-H*nu/KB/T_e));*/
  
  return x;
}

/*double bfheating_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the modified rate coefficient for photoionization 
/// using gsl integrators.
{
  double calculate_sahafact(int element, int ion, int level, double T, double E_threshold);
  double get_groundlevelpop(int cellnumber, int element, int ion);
  double photoionization_crosssection(double nu_edge, double nu);
  double radfield2(double nu, double T, double W);
  double x;
  
  int cellnumber = ((gslintegration_bfheatingparas *) paras)->cellnumber;
  double nu_edge = ((gslintegration_bfheatingparas *) paras)->nu_edge;
  
  float T_e = cell[cellnumber].T_e;
  float T_R = cell[cellnumber].T_R;
  float W = cell[cellnumber].W;
  float nne = cell[cellnumber].nne;
  
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  int element = mastate[tid].element;
  int ion  = mastate[tid].ion;
  int level = mastate[tid].level;
  double nnlevel = mastate[tid].nnlevel;
  double sigma_bf = photoionization_crosssection(nu_edge,nu);
  double E_threshold = nu_edge*H;
  double sf = calculate_sahafact(element,ion,level,T_e,E_threshold);
  double nnionlevel = get_groundlevelpop(cellnumber,element,ion+1);
  
  x = sigma_bf*(1-nu_edge/nu)*radfield2(nu,T_R,W) * (1-nnionlevel*nne/nnlevel*sf*exp(-H*nu/KB/T_e));
  return x;
}*/

/*double ffheating_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the free-free heating rate using gsl integrators.
{
  double ionstagepop(int cellnumber, int element, int ion);
  double radfield2(double nu, double T, double W);
  int get_ionstage(int element, int ion);
  int get_element(int element);
  
  double nne;//,nnion;//,nnlevel;
  double g_ff,kappa_ff;
  double T_R,T_D,W,W_D;
  double x;
  
  int element,ion;
  int nions,Z;
  
  double T_e = ((gslintegration_ffheatingparas *) paras)->T_e;
  int cellnumber = ((gslintegration_ffheatingparas *) paras)->cellnumber;
  
  nne = cell[cellnumber].nne;
  T_R = cell[cellnumber].T_R;
//  T_D = cell[cellnumber].T_D;
  W = cell[cellnumber].W;
//  W_D = cell[cellnumber].W_D;
  
  g_ff = 1;
  kappa_ff = 0.;
  for (element = 0; element < nelements; element++)
  {
    //Z = get_element(element); ///atomic number
    nions = get_nions(element);
    for (ion = 0; ion < nions; ion++)
    {
      /// Z is ionic charge in the following formula
      Z = get_ionstage(element,ion)-1;
      if (get_ionstage(element,ion) > 1)
      {
        //nnion = ionstagepop(cellnumber,element,ion);
        //kappa_ff += 3.69255e8 * pow(Z,2) / sqrt(T_e) * pow(nu,-3) * g_ff * nne * nnion * (1-exp(-HOVERKB*nu/T_e));
        
        kappa_ff += pow(Z,2) * ionstagepop(cellnumber,element,ion);
        
        //kappa_ff += ionstagepop(cellnumber,element,ion)*(1-exp(-HOVERKB*nu/T_e))*pow(Z,2)*nne* pow(nu,-3) ;
      }
    }
  }
  kappa_ff *= 3.69255e8 / sqrt(T_e) * pow(nu,-3) * g_ff * nne * (1-exp(-HOVERKB*nu/T_e));
//  if (nu <= nu_rfcut)
    x = kappa_ff * radfield2(nu,T_R,W);
//  else
//    x = kappa_ff * radfield2(nu,T_D,W_D);
  
  return x;
}*/

double bfcooling_integrand_gsl(double nu, void *paras)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximative way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  double photoionization_crosssection(double nu_edge, double nu);
  double x;
  
  float T  = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;
  
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection(nu_edge,nu);
  
  x = sigma_bf*(1-nu_edge/nu) * TWOHOVERCLIGHTSQUARED*pow(nu,3) * exp(-HOVERKB*nu/T);
  return x;
}

double bfcooling_integrand_gsl_2(double nu, void *paras)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximative way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  double photoionization_crosssection(double nu_edge, double nu);
  double x;
  
  float T  = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;
  
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection(nu_edge,nu);
  
  x = sigma_bf*(1/nu_edge-1/nu) * TWOOVERCLIGHTSQUARED*pow(nu,3) * exp(-HOVERKB*nu/T);
  return x;
}


double stimulated_bfcooling_integrand_gsl(double nu, void *paras)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximative way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  double photoionization_crosssection(double nu_edge, double nu);
  double radfield2(double nu, double T, double W);
  double x;
  
  float T  = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;
  
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection(nu_edge,nu);
  
  x = sigma_bf*(1-nu_edge/nu) * radfield2(nu, T, 1) * exp(-HOVERKB*nu/T);
  return x;
}


double stimulated_recomb_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the rate coefficient for spontaneous recombination
/// using gsl integrators.
{
  double photoionization_crosssection(double nu_edge, double nu);
  //double epsilon(int element, int ion, int level);
  double radfield2(double nu, double T, double W);
  double x,sigma_bf;
  
  //int element = mastate[tid].element;
  //int ion = mastate[tid].ion;
  //int level = mastate[tid].level;
  
  double T = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;
  
  //double nu_edge = (epsilon(element,ion+1,0)-epsilon(element,ion,level))/H;
  //printout("[debug] alpha_sp_integrand: element, ion, level: %d, %d, %d\n",exchangepkt_ptr->MA_element,exchangepkt_ptr->MA_ion,exchangepkt_ptr->MA_level);
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  sigma_bf = photoionization_crosssection(nu_edge,nu);
  x = sigma_bf/H/nu * radfield2(nu,T,1) * exp(-HOVERKB*nu/T);
  ///in formula this looks like
  ///x = sigma_bf/H/nu * 2*H*pow(nu,3)/pow(CLIGHT,2) * exp(-H*nu/KB/T);
  
  ///set contributions from Lyman continuum artificially to zero to overcome it's large opacity
  //if (exchangepkt_ptr->MA_level == 0) x = 0;
  return x;
}




///***************************************************************************/
/// The following functions are used to interpolate the rate coefficients
/// for a given temperature out of the precalculated values.
double interpolate_spontrecombcoeff(int element, int ion, int level, double T)
{
  double result;
  /*int lowerindex = floor((T-MINTEMP)/T_step);
  int upperindex = lowerindex + 1;
  double T_upper =  MINTEMP + upperindex*T_step;
  double T_lower =  MINTEMP + lowerindex*T_step;*/
  int lowerindex = floor(log(T/MINTEMP)/T_step_log);
  if (lowerindex < TABLESIZE-1)
  {
    int upperindex = lowerindex + 1;
    double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
    double T_upper =  MINTEMP * exp(upperindex*T_step_log);
    
    double f_upper = elements[element].ions[ion].levels[level].spontrecombcoeff[upperindex];
    double f_lower = elements[element].ions[ion].levels[level].spontrecombcoeff[lowerindex];
    //printout("interpolate_spontrecombcoeff element %d, ion %d, level %d, upper %g, lower %g\n",element,ion,level,f_upper,f_lower);
  
    result=f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
  }
  else result = elements[element].ions[ion].levels[level].spontrecombcoeff[TABLESIZE-1];

  return result;
}

// double interpolate_spontrecombcoeff_E(int element, int ion, int level, double T)
// {
//   /*int lowerindex = floor((T-MINTEMP)/T_step);
//   int upperindex = lowerindex + 1;
//   double T_upper =  MINTEMP + upperindex*T_step;
//   double T_lower =  MINTEMP + lowerindex*T_step;*/
//   int lowerindex = floor(log(T/MINTEMP)/T_step_log);
//   int upperindex = lowerindex + 1;
//   double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
//   double T_upper =  MINTEMP * exp(upperindex*T_step_log);
//   
//   double f_upper = elements[element].ions[ion].levels[level].spontrecombcoeff_E[upperindex];
//   double f_lower = elements[element].ions[ion].levels[level].spontrecombcoeff_E[lowerindex];
//   
//   return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
// }


// double interpolate_photoioncoeff_below(int element, int ion, int level, double T)
// {
//   int lowerindex = floor(log(T/MINTEMP)/T_step_log);
//   int upperindex = lowerindex + 1;
//   double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
//   double T_upper =  MINTEMP * exp(upperindex*T_step_log);
//   
//   double f_upper = elements[element].ions[ion].levels[level].photoioncoeff_below[upperindex];
//   double f_lower = elements[element].ions[ion].levels[level].photoioncoeff_below[lowerindex];
//   
//   return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
// }
// 
// double interpolate_photoioncoeff_above(int element, int ion, int level, double T)
// {
//   int lowerindex = floor(log(T/MINTEMP)/T_step_log);
//   int upperindex = lowerindex + 1;
//   double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
//   double T_upper =  MINTEMP * exp(upperindex*T_step_log);
//   
//   double f_upper = elements[element].ions[ion].levels[level].photoioncoeff_above[upperindex];
//   double f_lower = elements[element].ions[ion].levels[level].photoioncoeff_above[lowerindex];
//   
//   return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
// }


double interpolate_corrphotoioncoeff(int element, int ion, int level, double T)
{
  double result;
/*  int lowerindex = floor((T-MINTEMP)/T_step);
  int upperindex = lowerindex + 1;
  double T_upper =  MINTEMP + upperindex*T_step;
  double T_lower =  MINTEMP + lowerindex*T_step;*/
  int lowerindex = floor(log(T/MINTEMP)/T_step_log);
  if (lowerindex < TABLESIZE-1)
  {
    int upperindex = lowerindex + 1;
    double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
    double T_upper =  MINTEMP * exp(upperindex*T_step_log);
    
    double f_upper = elements[element].ions[ion].levels[level].corrphotoioncoeff[upperindex];
    double f_lower = elements[element].ions[ion].levels[level].corrphotoioncoeff[lowerindex];
    
    result=f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
  }
  else result=elements[element].ions[ion].levels[level].corrphotoioncoeff[TABLESIZE-1];
  
  return result;
}

// double interpolate_corrphotoioncoeff_above(int element, int ion, int level, double T)
// {
//   int lowerindex = floor(log(T/MINTEMP)/T_step_log);
//   int upperindex = lowerindex + 1;
//   double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
//   double T_upper =  MINTEMP * exp(upperindex*T_step_log);
//   
//   double f_upper = elements[element].ions[ion].levels[level].corrphotoioncoeff_above[upperindex];
//   double f_lower = elements[element].ions[ion].levels[level].corrphotoioncoeff_above[lowerindex];
//   
//   return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
// }


double interpolate_bfheatingcoeff(int element, int ion, int level, double T) // double T_e, double T_R)
{
  double result;
/*  int lowerindex = floor((T-MINTEMP)/T_step);
  int upperindex = lowerindex + 1;
  double T_upper =  MINTEMP + upperindex*T_step;
  double T_lower =  MINTEMP + lowerindex*T_step;*/
  int lowerindex = floor(log(T/MINTEMP)/T_step_log);
  if (lowerindex < TABLESIZE-1)
  {
    int upperindex = lowerindex + 1;
    double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
    double T_upper =  MINTEMP * exp(upperindex*T_step_log);
    
    double f_upper = elements[element].ions[ion].levels[level].bfheating_coeff[upperindex];
    double f_lower = elements[element].ions[ion].levels[level].bfheating_coeff[lowerindex];
    
    result=f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
  }
  else result=elements[element].ions[ion].levels[level].bfheating_coeff[TABLESIZE-1];
  
  return result;
}

// double interpolate_bfheatingcoeff_above(int element, int ion, int level, double T) // double T_e, double T_R)
// {
//   int lowerindex = floor(log(T/MINTEMP)/T_step_log);
//   int upperindex = lowerindex + 1;
//   double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
//   double T_upper =  MINTEMP * exp(upperindex*T_step_log);
//   
//   double f_upper = elements[element].ions[ion].levels[level].bfheating_coeff_above[upperindex];
//   double f_lower = elements[element].ions[ion].levels[level].bfheating_coeff_above[lowerindex];
//   
//   return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
// }


double interpolate_bfcoolingcoeff(int element, int ion, int level, double T)
{
  double result;
/*  int lowerindex = floor((T-MINTEMP)/T_step);
  int upperindex = lowerindex + 1;
  double T_upper =  MINTEMP + upperindex*T_step;
  double T_lower =  MINTEMP + lowerindex*T_step;*/
  int lowerindex = floor(log(T/MINTEMP)/T_step_log);
  if (lowerindex < TABLESIZE-1)
  {
    int upperindex = lowerindex + 1;
    double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
    double T_upper =  MINTEMP * exp(upperindex*T_step_log);
    
    double f_upper = elements[element].ions[ion].levels[level].bfcooling_coeff[upperindex];
    double f_lower = elements[element].ions[ion].levels[level].bfcooling_coeff[lowerindex];
    
    result=f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
  }
  else result=elements[element].ions[ion].levels[level].bfcooling_coeff[TABLESIZE-1];
  
  return result;
}

/*
double interpolate_stimulated_bfcoolingcoeff(int element, int ion, int level, double T)
{
  int lowerindex = floor(log(T/MINTEMP)/T_step_log);
  int upperindex = lowerindex + 1;
  double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
  double T_upper =  MINTEMP * exp(upperindex*T_step_log);
  
  double f_upper = elements[element].ions[ion].levels[level].stimulated_bfcooling_coeff[upperindex];
  double f_lower = elements[element].ions[ion].levels[level].stimulated_bfcooling_coeff[lowerindex];
  
  return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
}

double interpolate_stimulated_recomb(int element, int ion, int level, double T)
{
  int lowerindex = floor(log(T/MINTEMP)/T_step_log);
  int upperindex = lowerindex + 1;
  double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
  double T_upper =  MINTEMP * exp(upperindex*T_step_log);
  
  double f_upper = elements[element].ions[ion].levels[level].stimulated_recomb_coeff[upperindex];
  double f_lower = elements[element].ions[ion].levels[level].stimulated_recomb_coeff[lowerindex];
  
  return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
}
*/

/*double interpolate_zeta(int element, int ion, double T)
{
  int lowerindex = floor(log(T/MINTEMP)/T_step_log);
  int upperindex = lowerindex + 1;
  double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
  double T_upper =  MINTEMP * exp(upperindex*T_step_log);
  
  double f_upper = elements[element].ions[ion].zeta[upperindex];
  double f_lower = elements[element].ions[ion].zeta[lowerindex];
  
  return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
}*/

double interpolate_ions_spontrecombcoeff(int element, int ion, double T)
{
  double result;
  /*  int lowerindex = floor((T-MINTEMP)/T_step);
  int upperindex = lowerindex + 1;
  double T_upper =  MINTEMP + upperindex*T_step;
  double T_lower =  MINTEMP + lowerindex*T_step;*/
  int lowerindex = floor(log(T/MINTEMP)/T_step_log);
  if (lowerindex < TABLESIZE-1)
  {
    int upperindex = lowerindex + 1;
    double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
    double T_upper =  MINTEMP * exp(upperindex*T_step_log);
    
    double f_upper = elements[element].ions[ion].Alpha_sp[upperindex];
    double f_lower = elements[element].ions[ion].Alpha_sp[lowerindex];
    
    result=f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
  }
  else result=elements[element].ions[ion].Alpha_sp[TABLESIZE-1];
  
  return result;
}





///***************************************************************************/
double get_spontrecombcoeff(int element, int ion, int level, int modelgridindex)
/// Returns the rate coefficient for spontaneous recombination. Only needed during
/// packet propagation, therefore the value is taken from the 
/// cell history if known.
/// For ionisation to other levels than the ground level this must be adapted.
{
  double interpolate_spontrecombcoeff(int element, int ion, int level, double T);
  double alpha_sp;
  
  double T_e = get_Te(modelgridindex);
  
  if (use_cellhist >= 0)
  {
    alpha_sp = cellhistory[tid].chelements[element].chions[ion].chlevels[level].spontaneousrecombrate;
    /// Interpolate alpha_sp out of precalculated values
    if (alpha_sp < 0)
    {
      alpha_sp = interpolate_spontrecombcoeff(element,ion,level,T_e);
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].spontaneousrecombrate = alpha_sp;
    }
    /// Integrate alpha_sp directly. SLOW!!!
    /// ...
  }
  else
  {
    alpha_sp = interpolate_spontrecombcoeff(element,ion,level,T_e);
  }
  
  return alpha_sp;
}



///***************************************************************************/
double get_corrphotoioncoeff(int element, int ion, int level, int modelgridindex)
/// Returns the for stimulated emission corrected photoionisation rate coefficient.
/// Only needed during packet propagation, therefore the value is taken from the 
/// cell history if known.
#ifdef FORCE_LTE
{
  double interpolate_corrphotoioncoeff(int element, int ion, int level, double T);
  double gammacorr;
  
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  double T_R = get_TR(modelgridindex);
  
  if (use_cellhist >= 0)
  {
    gammacorr = cellhistory[tid].chelements[element].chions[ion].chlevels[level].corrphotoioncoeff;
    /// Interpolate gammacorr out of precalculated values
    if (gammacorr < 0)
    {
      gammacorr = interpolate_corrphotoioncoeff(element,ion,level,T_R);
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].corrphotoioncoeff = gammacorr;
    }
    /// Integrate gammacorr directly. SLOW!!!
    /// ...
  }
  else
  {
    gammacorr = interpolate_corrphotoioncoeff(element,ion,level,T_R);
  }
  
  return gammacorr;
}
#else
{
  double interpolate_corrphotoioncoeff(int element, int ion, int level, double T);
  double gammacorr;
  int index_in_groundlevelcontestimor;
  
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  double T_R = get_TR(modelgridindex);
  double W = get_W(modelgridindex);
  
  if (use_cellhist >= 0)
  {
    gammacorr = cellhistory[tid].chelements[element].chions[ion].chlevels[level].corrphotoioncoeff;
    /// Interpolate gammacorr out of precalculated values
    if (gammacorr < 0)
    {
      gammacorr = W*interpolate_corrphotoioncoeff(element,ion,level,T_R);
      index_in_groundlevelcontestimor = elements[element].ions[ion].levels[level].closestgroundlevelcont;
      if (index_in_groundlevelcontestimor >= 0) gammacorr *= corrphotoionrenorm[modelgridindex*nelements*maxion+index_in_groundlevelcontestimor];
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].corrphotoioncoeff = gammacorr;
    }
    /// Integrate gammacorr directly. SLOW!!!
    /// ...
  }
  else
  {
    gammacorr = W*interpolate_corrphotoioncoeff(element,ion,level,T_R);
    index_in_groundlevelcontestimor = elements[element].ions[ion].levels[level].closestgroundlevelcont;
    if (index_in_groundlevelcontestimor >= 0) gammacorr *= corrphotoionrenorm[modelgridindex*nelements*maxion+index_in_groundlevelcontestimor];
  }
  
  //if (gammacorr == 0) printout("histindex %d, element %d, ion %d, level %d, interpol %g, W %g, renorm %g, gammacorr %g\n",histindex,element,ion,level,interpolate_corrphotoioncoeff(element,ion,level,T_R),W,corrphotoionrenorm[modelgridindex*nelements*maxion+index_in_groundlevelcontestimor],gammacorr);
  return gammacorr;
  //double get_corrphotoioncoeff_ana(int element, int ion, int level, int cellnumber);
  //return get_corrphotoioncoeff_ana(element, ion, level, cellnumber);
}
#endif


#ifndef FORCE_LTE
///***************************************************************************/
double get_corrphotoioncoeff_ana(int element, int ion, int level, int modelgridindex)
/// Returns the for stimulated emission corrected photoionisation rate coefficient.
/// Only needed during packet propagation, therefore the value is taken from the 
/// cell history if known.
{
  double interpolate_corrphotoioncoeff(int element, int ion, int level, double T);
  double gammacorr;
  
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  double T_R = get_TR(modelgridindex);
  double W = get_W(modelgridindex);
  
  gammacorr = W*interpolate_corrphotoioncoeff(element,ion,level,T_R);
  
  return gammacorr;
}




///***************************************************************************/
double get_bfheatingcoeff(int element, int ion, int level, int modelgridindex)
{
  double interpolate_bfheatingcoeff(int element, int ion, int level, double T);
  double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
  double bfheating,nnlevel;
  
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  double T_R = get_TR(modelgridindex);
  double W = get_W(modelgridindex);
  int index_in_groundlevelcontestimor;
  
  /*nnlevel = calculate_exclevelpop(cellnumber,element,ion,level);
  bfheating = nnlevel * W * interpolate_bfheatingcoeff_below(element,ion,level,T_R);
  index_in_groundlevelcontestimor = elements[element].ions[ion].levels[level].closestgroundlevelcont;
  if (index_in_groundlevelcontestimor >= 0) bfheating *= bfheatingestimator[cellnumber*nelements*maxion+index_in_groundlevelcontestimor];*/
  bfheating = W * interpolate_bfheatingcoeff(element,ion,level,T_R);
  index_in_groundlevelcontestimor = elements[element].ions[ion].levels[level].closestgroundlevelcont;
  if (index_in_groundlevelcontestimor >= 0) bfheating *= bfheatingestimator[modelgridindex*nelements*maxion+index_in_groundlevelcontestimor];
  
  return bfheating;
  //double get_bfheatingcoeff_ana(int element, int ion, int level, int cellnumber);
  //return get_bfheatingcoeff_ana(element, ion, level, cellnumber);  
}


///***************************************************************************/
double get_bfheatingcoeff_ana(int element, int ion, int level, int modelgridindex)
{
  double interpolate_bfheatingcoeff(int element, int ion, int level, double T);
  double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
  double bfheating,nnlevel;
  
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  double T_R = get_TR(modelgridindex);
  double W = get_W(modelgridindex);
  int index_in_groundlevelcontestimor;
  
  /*nnlevel = calculate_exclevelpop(cellnumber,element,ion,level);
  bfheating = nnlevel * W * interpolate_bfheatingcoeff_below(element,ion,level,T_R);*/
  bfheating = W * interpolate_bfheatingcoeff(element,ion,level,T_R);
  
  return bfheating;
}


#endif



///***************************************************************************/
double get_bfcooling(int element, int ion, int level, int modelgridindex)
/// Returns the rate for bfheating. This can be called during packet propagation
/// or update_grid. Therefore we need to decide whether a cell history is
/// known or not.
{
  double interpolate_bfcoolingcoeff(int element, int ion, int level, double T);
  double get_groundlevelpop(int modelgridindex, int element, int ion);
  double ionstagepop(int modelgridindex, int element, int ion);
  double bfcooling;
  
  double T_e = get_Te(modelgridindex);
  double nne = get_nne(modelgridindex);
  //double nnionlevel = get_groundlevelpop(modelgridindex,element,ion+1);
  double nnion = ionstagepop(modelgridindex,element,ion+1);
  
  if (use_cellhist >= 0)
  {
    bfcooling = cellhistory[tid].chelements[element].chions[ion].chlevels[level].bfcooling;
    /// Interpolate bfcooling out of precalculated values
    if (bfcooling < 0)
    {
      //bfcooling = interpolate_bfcoolingcoeff(element,ion,level,T_e) * nnionlevel*nne;
      bfcooling = interpolate_bfcoolingcoeff(element,ion,level,T_e) * nnion*nne;
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].bfcooling = bfcooling;
    }
    /// Direct integration
    /// ...
  }
  else
  {
    /// Interpolate bfcooling out of precalculated values
    //bfcooling = interpolate_bfcoolingcoeff(element,ion,level,T_e) * nnionlevel*nne;
    bfcooling = interpolate_bfcoolingcoeff(element,ion,level,T_e) * nnion*nne;
    /// Direct integration
    /// ...
  }
  
  #ifdef DEBUG_ON
    if (!isfinite(bfcooling))
    {
      printout("[fatal] get_bfcooling: bfcooling infinite (%g) for element %d, ion %d, level %d in modelgridcell %d\n",bfcooling,element,ion,level,modelgridindex);
      printout("[fatal] get_bfcooling: bfcoolingcoeff %g, nnion %g, nne %g, T_e %g\n",interpolate_bfcoolingcoeff(element,ion,level,T_e),nnion,nne,T_e);
    }
  #endif
  
  return bfcooling;
}


///***************************************************************************/
/*double get_spontrecomb(int element, int ion, int level, int cellnumber)
/// Returns the rate coefficient for spontaneous recombination. Only needed during
/// packet propagation, therefore the value is taken from the 
/// cell history if known.
/// For ionisation to other levels than the ground level this must be adapted.
{
  double interpolate_spontrecombcoeff(int element, int ion, int level, double T);
  double get_groundlevelpop(int cellnumber, int element, int ion);
  double alpha_sp;
  
  double T_e = cell[cellnumber].T_e;
  double nne = cell[cellnumber].nne;
  double nnionlevel = get_groundlevelpop(cellnumber,element,ion+1);
  
  alpha_sp = interpolate_spontrecombcoeff(element,ion,level,T_e) * nnionlevel*nne;
  
  return alpha_sp;
}*/





///***************************************************************************/
/*
void check_interpolation(double T_min, double T_max)
/// Function writes "exact" and interpolated values of the rate coefficients
/// to alpha_sp_file, alpha_st_file and gamma_file to judge the quality
/// of the interpolation. Both integrator types are used.
{
  double calculate_sahafact(int element, int ion, int level, double T, double E_threshold);
  double epsilon(int element, int ion, int level);

  double interpolate_alpha_sp(int element, int ion, int level, double T_e);
  double interpolate_modified_alpha_sp(int element, int ion, int level, double T_e);
  double interpolate_alpha_st(int element, int ion, int level, double T_e);
  double interpolate_modified_alpha_st(int element, int ion, int level, double T_e);
  double interpolate_gamma(int element, int ion, int level, double T_e);
  double interpolate_modified_gamma(int element, int ion, int level, double T_e);
  
  double alpha_sp_integrand_gsl(double nu, void * paras);
  double modified_alpha_sp_integrand_gsl(double nu, void * paras);
  double alpha_st_integrand_gsl(double nu, void * paras);
  double modified_alpha_st_integrand_gsl(double nu, void * paras);
  double gamma_integrand_gsl(double nu, void * paras);
  double modified_gamma_integrand_gsl(double nu, void * paras);
  
  gslintegration_paras intparas; 
  gsl_function F1,F2,F3,F4,F5,F6;
  F1.function = &alpha_sp_integrand_gsl;
  F2.function = &modified_alpha_sp_integrand_gsl;
  F3.function = &alpha_st_integrand_gsl;
  F4.function = &modified_alpha_st_integrand_gsl;
  F5.function = &gamma_integrand_gsl;
  F6.function = &modified_gamma_integrand_gsl;
  
  float qromb(float (*func)(float), float a, float b);
  float alpha_sp_integrand(float nu);
  float modified_alpha_sp_integrand(float nu);
  float alpha_st_integrand(float nu);
  float modified_alpha_st_integrand(float nu);
  float gamma_integrand(float nu);
  float modified_gamma_integrand(float nu);
  
  
  double alpha_sp,modified_alpha_sp,alpha_st,modified_alpha_st,gamma,modified_gamma;
  double alpha_sp_gsl,modified_alpha_sp_gsl,alpha_st_gsl,modified_alpha_st_gsl,gamma_gsl,modified_gamma_gsl;
  double alpha_sp_nr,modified_alpha_sp_nr,alpha_st_nr,modified_alpha_st_nr,gamma_nr,modified_gamma_nr;
  double error;
  size_t neval;
  double T_e,sf;
  int level;
  int iter;
  
  
  FILE *alpha_sp_file;
  FILE *alpha_st_file;
  FILE *gamma_file;
  if ((alpha_sp_file = fopen("alpha_sp.out", "w")) == NULL)
  {
    printf("Cannot open alpha_sp.out.\n");
    exit(0);
  }
  setvbuf(alpha_sp_file, NULL, _IOLBF, 1);
  if ((alpha_st_file = fopen("alpha_st.out", "w")) == NULL)
  {
    printf("Cannot open alpha_st.out.\n");
    exit(0);
  }
  setvbuf(alpha_st_file, NULL, _IOLBF, 1);
  if ((gamma_file = fopen("gamma.out", "w")) == NULL)
  {
    printf("Cannot open gamma.out.\n");
    exit(0);
  }
  setvbuf(gamma_file, NULL, _IOLBF, 1);
  
  /// Works so far only for hydrogen or the first ionisation stage of any element!
  mastate[tid].element = 0;
  mastate[tid].ion = 0;
  for (level = 0; level < get_nlevels(0,0); level++)
  {
    mastate[tid].level = level;
    double E_threshold = epsilon(0,1,0) - epsilon(0,0,level);
    double nu_threshold = E_threshold/H;
    intparas.nu_edge = nu_threshold;
    alpha_sp_integrand_parameters.nu_edge = nu_threshold;
    double tstep = (T_max-T_min)/99;
    for (iter = 0; iter < 100; iter++)
    {
      T_e = T_min + iter*tstep;
      sf = calculate_sahafact(0,0,level,T_e,E_threshold);
      
      /// calculate from tabulated values
      alpha_sp = interpolate_alpha_sp(0,0,level,T_e);
      modified_alpha_sp = interpolate_modified_alpha_sp(0,0,level,T_e);
      alpha_st = interpolate_alpha_st(0,0,level,T_e)*sf;  /// The sahafactor was taken out of the precalculation
      modified_alpha_st = interpolate_modified_alpha_st(0,0,level,T_e)*sf; /// The sahafactor was taken out of the precalculation
      gamma = interpolate_gamma(0,0,level,T_e);
      modified_gamma = interpolate_modified_gamma(0,0,level,T_e);
      
      /// calculate with gsl integrators
      intparas.T = T_e;
      F1.params = &intparas;
      F2.params = &intparas;
      F3.params = &intparas;
      F4.params = &intparas;
      F5.params = &intparas;
      F6.params = &intparas;
      gsl_integration_qng(&F1, nu_threshold, 10*nu_threshold, 0, 1e-2, &alpha_sp_gsl, &error, &neval);
      gsl_integration_qng(&F2, nu_threshold, 10*nu_threshold, 0, 1e-2, &modified_alpha_sp_gsl, &error, &neval);
      gsl_integration_qng(&F3, nu_threshold, 10*nu_threshold, 0, 1e-2, &alpha_st_gsl, &error, &neval);
      gsl_integration_qng(&F4, nu_threshold, 10*nu_threshold, 0, 1e-2, &modified_alpha_st_gsl, &error, &neval);
      gsl_integration_qng(&F5, nu_threshold, 10*nu_threshold, 0, 1e-2, &gamma_gsl, &error, &neval);
      gsl_integration_qng(&F6, nu_threshold, 10*nu_threshold, 0, 1e-2, &modified_gamma_gsl, &error, &neval);
      alpha_sp_gsl *= FOURPI * sf;
      modified_alpha_sp_gsl *= FOURPI * sf;
      alpha_st_gsl *= FOURPI * sf;
      modified_alpha_st_gsl *= FOURPI * sf;
      gamma_gsl *= FOURPI;
      modified_gamma_gsl *= FOURPI;
      
      /// calculate with qromb integrator of NR
      alpha_sp_integrand_parameters.T = T_e;
      alpha_sp_nr = qromb(alpha_sp_integrand,nu_threshold,10*nu_threshold);
      modified_alpha_sp_nr = qromb(modified_alpha_sp_integrand,nu_threshold,10*nu_threshold);
      alpha_st_nr = qromb(alpha_st_integrand,nu_threshold,10*nu_threshold);
      modified_alpha_st_nr = qromb(modified_alpha_st_integrand,nu_threshold,10*nu_threshold);
      gamma_nr = qromb(gamma_integrand,nu_threshold,10*nu_threshold);
      modified_gamma_nr = qromb(modified_gamma_integrand,nu_threshold,10*nu_threshold);
      alpha_sp_nr *= FOURPI * sf;
      modified_alpha_sp_nr *= FOURPI * sf;
      alpha_st_nr *= FOURPI * sf;
      modified_alpha_st_nr *= FOURPI * sf;
      gamma_nr *= FOURPI;
      modified_gamma_nr *= FOURPI;
      
      fprintf(alpha_sp_file,"%g %g %g %g %g %g %g\n", T_e,alpha_sp,alpha_sp_gsl,alpha_sp_nr,modified_alpha_sp,modified_alpha_sp_gsl,modified_alpha_sp_nr);
      fprintf(alpha_st_file,"%g %g %g %g %g %g %g\n", T_e,alpha_st,alpha_st_gsl,alpha_st_nr,modified_alpha_st,modified_alpha_st_gsl,modified_alpha_st_nr);
      fprintf(gamma_file,"%g %g %g %g %g %g %g\n",T_e,gamma,gamma_gsl,gamma_nr,modified_gamma,modified_gamma_gsl,modified_gamma_nr);
    }
  }
  fclose(alpha_sp_file);
  fclose(alpha_st_file);
  fclose(gamma_file);
  exit(0);
}
*/


