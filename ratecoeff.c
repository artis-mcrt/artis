#include "sn3d.h"
#include "atomic.h"
#include "ltepop.h"
#include "macroatom.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "grid_init.h"
#include <gsl/gsl_integration.h>


typedef struct
{
  int modelgridindex;
  double nu_edge;
  float *photoion_xs;
} gsl_integral_paras_gammacorr;


// private functions
void calculate_rate_coefficients(void);
void write_ratecoeff_dat(void);
void calculate_ion_alpha_sp(void);
double gamma_integrand_gsl(double nu, void *paras);
double gammacorr_integrand_gsl(double nu, void *paras);
double approx_bfheating_integrand_gsl(double nu, void *paras);
double bfcooling_integrand_gsl(double nu, void *paras);
double bfcooling_integrand_gsl_2(double nu, void *paras);
double stimulated_bfcooling_integrand_gsl(double nu, void *paras);
double stimulated_recomb_integrand_gsl(double nu, void *paras);
double calculate_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double gammacorr_integrand_gsl_radfield(double nu, void *paras);
double calculate_bfheatingcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double bfheating_integrand_gsl_radfield(double nu, void *paras);


///****************************************************************************
void tabulate_ratecoefficients_gsl(void)
/// Precalculates the rate coefficients for stimulated and spontaneous
/// recombination and photoionisation on a given temperature grid using
/// libgsl integrators.
/// NB: with the nebular approximation they only depend on T_e, T_R and W.
/// W is easily factored out. For stimulated recombination we must assume
/// T_e = T_R for this precalculation.
{
  //size_t neval;   /// for qng integrator
  int calculate;

  /// Determine the temperture grids gridsize
  T_step = (1.*MAXTEMP-MINTEMP) / (TABLESIZE-1.);               /// global variables
  T_step_log = (log(MAXTEMP)-log(MINTEMP)) / (TABLESIZE-1.);

  /// Then readin the precalculated rate coefficients from file
  FILE *ratecoeff_file;
  if ((ratecoeff_file = fopen("ratecoeff.dat", "r")) != NULL)
  {
    /// Check whether current temperature range and atomic data match
    /// the precalculated rate coefficients
    int check;
    float T_min,T_max;
    int dum1;
    fscanf(ratecoeff_file,"%g %g %d",&T_min,&T_max,&dum1);
    printout("Tmin %g Tmax %g TABLESIZE %d \n",T_min,T_max,dum1);
    if (T_min == MINTEMP && T_max == MAXTEMP && dum1 == TABLESIZE)
    {
      check = 1;
      for (int element = 0; element < nelements; element++)
      {
        int nions = get_nions(element);
        for (int ion = 0; ion < nions; ion++)
        {
          int dum2,dum3,dum4;
          fscanf(ratecoeff_file,"%d %d %d %d",&dum1,&dum2,&dum3,&dum4);
          int nlevels = get_nlevels(element,ion);
          int ionisinglevels = get_ionisinglevels(element,ion);
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
      for (int element = 0; element < nelements; element++)
      {
        int nions = get_nions(element) - 1;
        for (int ion = 0; ion < nions; ion++)
        {
          //nlevels = get_nlevels(element,ion);
          int nlevels = get_ionisinglevels(element,ion); /// number of ionising levels associated with current ion
          int nbfcont = get_bfcontinua(element,ion);     /// number of ionising levels of the current ion which are used in the simulation
          for (int level = 0; level < nlevels; level++)
          {
            /// Loop over the phixs target states
            int nphixstargets = get_nphixstargets(element,ion,level);
            for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++)
            {
              /// Loop over the temperature grid
              for (int iter = 0; iter < TABLESIZE; iter++)
              {
                double alpha_sp,bfcooling_coeff,gammacorr,bfheating_coeff;
                fscanf(ratecoeff_file,"%lg %lg %lg %lg\n", &alpha_sp,&bfcooling_coeff,&gammacorr,&bfheating_coeff);

                if (level < nbfcont)
                {
                  elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].spontrecombcoeff[iter] = alpha_sp;
                  elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[iter] = bfcooling_coeff;
                  elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].corrphotoioncoeff[iter] = gammacorr;
                  elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfheating_coeff[iter] = bfheating_coeff;
                }
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
    calculate_rate_coefficients();
  }
  calculate_ion_alpha_sp();
}

void calculate_rate_coefficients(void)
{
  double intaccuracy = 1e-2;        /// Fractional accuracy of the integrator

  /// Calculate the rate coefficients for each level of each ion of each element
  for (int element = 0; element < nelements; element++)
  {
    int nions = get_nions(element) - 1;
    #ifdef _OPENMP
      #pragma omp parallel for
    #endif
    for (int ion = 0; ion < nions; ion++)
    {
      //nlevels = get_nlevels(element,ion);
      int nlevels = get_ionisinglevels(element,ion);
      /// That's only an option for pure LTE
      //if (TAKE_N_BFCONTINUA < nlevels) nlevels = TAKE_N_BFCONTINUA;
      printout("Performing rate integrals for Z = %d, ionstage %d...\n",elements[element].anumber,ion+1);

      gsl_integration_workspace *w;
      w = gsl_integration_workspace_alloc(1000);
      mastate[tid].element = element;   /// Global variable which passes the current element to all subfunctions of macroatom.c
      mastate[tid].ion = ion;   /// Global variable which passes the current ion to all subfunctions of macroatom.c
      for (int level = 0; level < nlevels; level++)
      {
        if ((level > 0) && (level % 10 == 0))
          printout("  completed up to level %d of %d\n",level,nlevels);

        int nphixstargets = get_nphixstargets(element,ion,level);
        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++)
        {
          int upperlevel = get_phixsupperlevel(element,ion,level,phixstargetindex);
          float phixstargetprobability = get_phixsprobability(element,ion,level,phixstargetindex);

          //printout("element %d, ion %d, level %d, upperlevel %d, epsilon %g, continuum %g, nlevels %d\n",element,ion,level,upperlevel,epsilon(element,ion,level),epsilon(element,ion+1,upperlevel),nlevels);

          mastate[tid].level = level;                   // Global variable which passes the current level to all subfunctions of macroatom.c
          double E_threshold = epsilon(element,ion+1,upperlevel) - epsilon(element,ion,level);
          double nu_threshold = E_threshold / H;
          double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table
          gslintegration_paras intparas;
          intparas.nu_edge = nu_threshold;              // Global variable which passes the threshold to the integrator
                                                        // the threshold of the first target gives nu of the first phixstable point
          // Loop over the temperature grid
          for (int iter = 0; iter < TABLESIZE; iter++)
          {
            double error;
            double T_e = MINTEMP * exp(iter * T_step_log);
            //T_e = MINTEMP + iter*T_step;
            double sfac = calculate_sahafact(element,ion,level,phixstargetindex,T_e,E_threshold);
            //printout("%d %g\n",iter,T_e);

            intparas.T = T_e;

            //gsl_function F_gamma;
            //F_gamma.function = &gamma_integrand_gsl;
            //F_gamma.params = &intparas;
            //gsl_function F_alpha_sp_E;
            //F_alpha_sp_E.function = &alpha_sp_E_integrand_gsl;
            //F_alpha_sp_E.params = &intparas;
            //F_stimulated_bfcooling.function = &stimulated_bfcooling_integrand_gsl;
            //F_stimulated_bfcooling.params = &intparas;
            //F_stimulated_recomb.function = &stimulated_recomb_integrand_gsl;
            //F_stimulated_recomb.params = &intparas;

            //TODO: These integrals will be zero for photoionization processes with all zero cross sections.
            //We could speed this up by detecting this case and skipping the integrals,
            //although the integrator is probably pretty fast in these cases anyway

            /// Spontaneous recombination and bf-cooling coefficient don't depend on the cutted radiation field
            double alpha_sp = 0.0;
            gsl_function F_alpha_sp;
            F_alpha_sp.function = &alpha_sp_integrand_gsl;
            F_alpha_sp.params = &intparas;
            gsl_integration_qag(&F_alpha_sp, nu_threshold, nu_max_phixs, 0, intaccuracy, 1000, 6, w, &alpha_sp, &error);
            alpha_sp *= FOURPI * sfac * phixstargetprobability;

            //if (iter == 0)
            //  printout("alpha_sp: element %d ion %d level %d upper level %d at temperature %g, alpha_sp is %g (integral %g, sahafac %g)\n", element, ion, level, upperlevel, T_e, alpha_sp, alpha_sp/(FOURPI * sfac * phixstargetprobability),sfac);

            double gammacorr = 0.0;
            gsl_function F_gammacorr;
            F_gammacorr.function = &gammacorr_integrand_gsl;
            F_gammacorr.params = &intparas;
            gsl_integration_qag(&F_gammacorr, nu_threshold, nu_max_phixs, 0, intaccuracy, 1000, 6, w, &gammacorr, &error);
            gammacorr *= FOURPI * phixstargetprobability;

            double bfheating_coeff = 0.0;
            gsl_function F_bfheating;
            F_bfheating.function = &approx_bfheating_integrand_gsl;
            F_bfheating.params = &intparas;
            gsl_integration_qag(&F_bfheating, nu_threshold, nu_max_phixs, 0, intaccuracy, 1000, 6, w, &bfheating_coeff, &error);
            bfheating_coeff *= FOURPI * phixstargetprobability;

            double bfcooling_coeff = 0.0;
            gsl_function F_bfcooling;
            F_bfcooling.function = &bfcooling_integrand_gsl;
            F_bfcooling.params = &intparas;
            gsl_integration_qag(&F_bfcooling, nu_threshold, nu_max_phixs, 0, intaccuracy, 1000, 6, w, &bfcooling_coeff, &error);
            bfcooling_coeff *= FOURPI * sfac * phixstargetprobability;

            /// Save the calculated coefficients to memory
            elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].spontrecombcoeff[iter] = alpha_sp;
            elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[iter] = bfcooling_coeff;
            elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].corrphotoioncoeff[iter] = gammacorr;
            elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfheating_coeff[iter] = bfheating_coeff;
          }
        }
      }
      gsl_integration_workspace_free(w);
    }
  }

  /// And the master process writes them to file in a serial operation
  if (rank_global == 0)
  {
    write_ratecoeff_dat();
  }
}

void write_ratecoeff_dat(void)
{
  FILE *ratecoeff_file;
  if ((ratecoeff_file = fopen("ratecoeff.dat", "w")) == NULL)
  {
    printout("Cannot open ratecoeff.dat\n");
    exit(0);
  }
  fprintf(ratecoeff_file,"%g %g %d\n",MINTEMP,MAXTEMP,TABLESIZE);
  for (int element = 0; element < nelements; element++)
  {
    int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      fprintf(ratecoeff_file,"%d %d %d %d\n",get_element(element), get_ionstage(element,ion), get_nlevels(element,ion),  get_ionisinglevels(element,ion));
    }
  }

  for (int element = 0; element < nelements; element++)
  {
    int nions = get_nions(element) - 1;
    for (int ion = 0; ion < nions; ion++)
    {
      //nlevels = get_nlevels(element,ion);
      int nlevels = get_ionisinglevels(element,ion);
      for (int level = 0; level < nlevels; level++)
      {
        /// Loop over the phixs targets
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
        {
          /// Loop over the temperature grid
          for (int iter = 0; iter < TABLESIZE; iter++)
          {
            double alpha_sp = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].spontrecombcoeff[iter];
            double bfcooling_coeff = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[iter];
            double gammacorr = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].corrphotoioncoeff[iter];
            double bfheating_coeff = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfheating_coeff[iter];
            //fprintf(ratecoeff_file,"%g %g %g %g %g %g %g %g %g\n", alpha_sp,alpha_sp_E,bfcooling_coeff,gamma_below,gamma_above,gammacorr_below,gammacorr_above,bfheating_coeff_below,bfheating_coeff_above);
            fprintf(ratecoeff_file,"%g %g %g %g\n", alpha_sp,bfcooling_coeff,gammacorr,bfheating_coeff);
          }
        }
      }
    }
  }
  fclose(ratecoeff_file);
}

//calculate the ion total recombination coefficients
void calculate_ion_alpha_sp(void)
{
  for (int iter = 0; iter < TABLESIZE; iter++)
  {
    double T_e = MINTEMP * exp(iter*T_step_log);
    //T_e = MINTEMP + iter*T_step;
    for (int element = 0; element < nelements; element++)
    {
      int nions = get_nions(element);
      for (int ion = 0; ion < nions-1; ion++)
      {
        //nlevels = get_nlevels(element,ion);
        //nlevels = get_ionisinglevels(element,ion); ///number of levels of the current ion which have an associated photoion cross section
        int nlevels = get_bfcontinua(element,ion); /// number of ionising levels used in the simulation
        int nlevelsupperion = get_nlevels(element,ion+1);
        double zeta = 0.;
        double pfunc = 0.;
        for (int upperlevel = 0; upperlevel < nlevelsupperion; upperlevel++)
        {
          double bfac = exp(-(epsilon(element,ion+1,upperlevel) - epsilon(element,ion+1,0))/KB/T_e); // Boltzmann factor
          pfunc += stat_weight(element,ion+1,upperlevel) * bfac; // partition function
        }
        for (int level = 0; level < nlevels; level++)
        {
          int nphixstargets = get_nphixstargets(element,ion,level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++)
          {
            int upperlevel = get_phixsupperlevel(element,ion,level,phixstargetindex);
            double bfac = exp(-(epsilon(element,ion+1,upperlevel) - epsilon(element,ion+1,0))/KB/T_e); // Boltzmann factor
            zeta += interpolate_spontrecombcoeff(element,ion,level,phixstargetindex,T_e) * bfac *
                    stat_weight(element,ion+1,upperlevel);
          }
        }
        zeta /= pfunc;
        elements[element].ions[ion].Alpha_sp[iter] = zeta;
//        zeta = interpolate_spontrecombcoeff(element,ion,0,T_e) / zeta;
//        elements[element].ions[ion].zeta[iter] = zeta;
        //printout("Alpha result: element %d ion %d temperature %g Alpha %g\n",element,ion,T_e,zeta);
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
  double sigma_bf = photoionization_crosssection(nu_edge,nu);
  double x = TWOOVERCLIGHTSQUARED * sigma_bf * pow(nu,2) * exp(-HOVERKB*nu/T);
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
  double sigma_bf = photoionization_crosssection(nu_edge,nu);
  double x = TWOOVERCLIGHTSQUARED * sigma_bf * pow(nu,3)/nu_edge * exp(-HOVERKB*nu/T);
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
  double T = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;

  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection(nu_edge,nu);

  /// Dependence on dilution factor W is linear. This allows to set it here to
  /// 1. and scale to its actual value later on.
  double x = sigma_bf / H / nu * radfield2(nu,T,1.);
  //x = sigma_bf/H/nu * radfield2(nu,T,1.);
  //if (HOVERKB*nu/T < 1e-2) x = sigma_bf * pow(nu,2)/(HOVERKB*nu/T);
  //else if (HOVERKB*nu/T >= 1e2) x = sigma_bf * pow(nu,2)*exp(-HOVERKB*nu/T);
  //else x = sigma_bf * pow(nu,2)/(exp(HOVERKB*nu/T)-1);

  return x;
}

double gammacorr_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the rate coefficient for photoionization
/// using gsl integrators. Corrected for stimulated recombination.
{
  double T = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;

  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection(nu_edge,nu);

  /// Dependence on dilution factor W is linear. This allows to set it here to
  /// 1. and scale to its actual value later on.
  /// Assumption T_e = T_R makes n_kappa/n_i * (n_i/n_kappa)* = 1
  return sigma_bf / H / nu * radfield2(nu,T,1.) * (1 - exp(-H*nu/KB/T));
}

double approx_bfheating_integrand_gsl(double nu, void *paras)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximative way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;
  double sigma_bf = photoionization_crosssection(nu_edge,nu);

  /// Precalculation for T_e=T_R and W=1
  float T  = ((gslintegration_paras *) paras)->T;
  double x = sigma_bf * (1-nu_edge/nu) * radfield2(nu,T,1) * (1-exp(-H*nu/KB/T));

  /// Precalculation for a (T_R,T_e)-grid, but still W is assumed to be 1.
  /// The radfield part can be corrected later because of its linear dependence.
  /// But not the W in the stimulated correction term!
  /*float T_e  = ((gslintegration_paras *) paras)->T;
  float T_R  = ((gslintegration_paras *) paras)->T2;
  int element = mastate[tid].element;
  int ion  = mastate[tid].ion;
  int level = mastate[tid].level;
  double E_threshold = nu_edge*H;
  double sf_Te = calculate_sahafact(element,ion,level,phixstargetindex,T_e,E_threshold);
  double sf_TR = calculate_sahafact(element,ion,level,phixstargetindex,T_R,E_threshold);
  x = sigma_bf*(1-nu_edge/nu)*radfield2(nu,T_R,1) * (1 - sqrt(T_e/T_R) * sf_Te/sf_TR * exp(-H*nu/KB/T_e));*/

  return x;
}

/*double bfheating_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the modified rate coefficient for photoionization
/// using gsl integrators.
{
  double get_groundlevelpop(int cellnumber, int element, int ion);
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
  double sfac = calculate_sahafact(element,ion,level,phixstargetindex,T_e,E_threshold);
  double nnionlevel = get_groundlevelpop(cellnumber,element,ion+1);

  x = sigma_bf*(1-nu_edge/nu)*radfield2(nu,T_R,W) * (1-nnionlevel*nne/nnlevel*sf*exp(-H*nu/KB/T_e));
  return x;
}*/

/*double ffheating_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the free-free heating rate using gsl integrators.
{
  double ionstagepop(int cellnumber, int element, int ion);

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
  float T = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;

  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection(nu_edge,nu);

  return sigma_bf * (1-nu_edge/nu) * TWOHOVERCLIGHTSQUARED * pow(nu,3) * exp(-HOVERKB*nu/T);
}

double bfcooling_integrand_gsl_2(double nu, void *paras)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximative way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  float T  = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;

  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection(nu_edge,nu);

  return sigma_bf*(1/nu_edge-1/nu) * TWOOVERCLIGHTSQUARED*pow(nu,3) * exp(-HOVERKB*nu/T);
}


double stimulated_bfcooling_integrand_gsl(double nu, void *paras)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximate way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  float T  = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;

  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection(nu_edge,nu);

  return sigma_bf * (1-nu_edge/nu) * radfield2(nu, T, 1) * exp(-HOVERKB*nu/T);
}


double stimulated_recomb_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the rate coefficient for spontaneous recombination
/// using gsl integrators.
{
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
  double sigma_bf = photoionization_crosssection(nu_edge,nu);
  double x = sigma_bf / H / nu * radfield2(nu,T,1) * exp(-HOVERKB*nu/T);
  ///in formula this looks like
  ///x = sigma_bf/H/nu * 2*H*pow(nu,3)/pow(CLIGHT,2) * exp(-H*nu/KB/T);

  ///set contributions from Lyman continuum artificially to zero to overcome it's large opacity
  //if (exchangepkt_ptr->MA_level == 0) x = 0;
  return x;
}




///***************************************************************************/
/// The following functions are used to interpolate the rate coefficients
/// for a given temperature out of the precalculated values.
double interpolate_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, double T)
{
  /*int lowerindex = floor((T-MINTEMP)/T_step);
  int upperindex = lowerindex + 1;
  double T_upper =  MINTEMP + upperindex*T_step;
  double T_lower =  MINTEMP + lowerindex*T_step;*/
  int lowerindex = floor(log(T / MINTEMP) / T_step_log);
  if (lowerindex < TABLESIZE - 1)
  {
    int upperindex = lowerindex + 1;
    double T_lower =  MINTEMP * exp(lowerindex*T_step_log);
    double T_upper =  MINTEMP * exp(upperindex*T_step_log);

    double f_upper = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].spontrecombcoeff[upperindex];
    double f_lower = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].spontrecombcoeff[lowerindex];
    //printout("interpolate_spontrecombcoeff element %d, ion %d, level %d, upper %g, lower %g\n",element,ion,level,f_upper,f_lower);

    return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
  }
  else
    return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].spontrecombcoeff[TABLESIZE-1];
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


double interpolate_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, double T)
{
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

    double f_upper = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].corrphotoioncoeff[upperindex];
    double f_lower = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].corrphotoioncoeff[lowerindex];

    return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
  }
  else
    return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].corrphotoioncoeff[TABLESIZE-1];
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


double interpolate_bfheatingcoeff(int element, int ion, int level, int phixstargetindex, double T) // double T_e, double T_R)
{
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

    double f_upper = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfheating_coeff[upperindex];
    double f_lower = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfheating_coeff[lowerindex];

    return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
  }
  else
    return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfheating_coeff[TABLESIZE-1];
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


double interpolate_bfcoolingcoeff(int element, int ion, int level, int phixstargetindex, double T)
{
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

    double f_upper = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[upperindex];
    double f_lower = elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[lowerindex];

    return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
  }
  else
    return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].bfcooling_coeff[TABLESIZE-1];
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

/* multiply by LTE ion population and nne electron density to get # recombs/sec */
double interpolate_ions_spontrecombcoeff(int element, int ion, double T)
{
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

    return f_lower + (f_upper-f_lower)/(T_upper-T_lower) * (T-T_lower);
  }
  else
    return elements[element].ions[ion].Alpha_sp[TABLESIZE-1];
}



///***************************************************************************/
double get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex)
/// Returns the rate coefficient for spontaneous recombination. Only needed during
/// packet propagation, therefore the value is taken from the
/// cell history if known.
/// For ionisation to other levels than the ground level this must be adapted.
{
  double alpha_sp = -1.0;

  if (use_cellhist)
  {
    alpha_sp = cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].spontaneousrecombrate;
    /// Interpolate alpha_sp out of precalculated values
  }

  if (alpha_sp < 0. || !use_cellhist)
  {
    double T_e = get_Te(modelgridindex);
    alpha_sp = interpolate_spontrecombcoeff(element,ion,level,phixstargetindex,T_e);

    if (use_cellhist)
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].spontaneousrecombrate = alpha_sp;
  }

  return alpha_sp;
}



///***************************************************************************/
double get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex)
/// Returns the for stimulated emission corrected photoionisation rate coefficient.
/// Only needed during packet propagation, therefore the value is taken from the
/// cell history if known.
{
  double gammacorr = -1.0;
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  if (use_cellhist)
  {
    gammacorr = cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].corrphotoioncoeff;
  }

  if (gammacorr < 0 || !use_cellhist)
  {
  #ifdef FORCE_LTE
    /// Interpolate gammacorr out of precalculated values
    double T_R = get_TR(modelgridindex);
    gammacorr = interpolate_corrphotoioncoeff(element,ion,level,phixstargetindex,T_R);
  #else
  # ifdef NO_LUT_PHOTOION
    //double W = get_W(modelgridindex);
    gammacorr = calculate_corrphotoioncoeff(element,ion,level,phixstargetindex,modelgridindex);
    //double oldgammacorr = W * interpolate_corrphotoioncoeff(element,ion,level,phixstargetindex,T_R);
    //if (fabs(gammacorr/oldgammacorr - 1.0) > 0.5)
    //  printout("LUKETEST: New corrphotoion coeff: %g, old value: %g\n",gammacorr,oldgammacorr);
  # else
    double W = get_W(modelgridindex);
    double T_R = get_TR(modelgridindex);

    gammacorr = W * interpolate_corrphotoioncoeff(element,ion,level,phixstargetindex,T_R);
    int index_in_groundlevelcontestimator = elements[element].ions[ion].levels[level].closestgroundlevelcont;
    if (index_in_groundlevelcontestimator >= 0)
      gammacorr *= corrphotoionrenorm[modelgridindex*nelements*maxion+index_in_groundlevelcontestimator];

  # endif
  #endif
    if (use_cellhist)
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].corrphotoioncoeff = gammacorr;
  }

  return gammacorr;
  //return get_corrphotoioncoeff_ana(element, ion, level, cellnumber);
}


#ifndef FORCE_LTE
///***************************************************************************/
double get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex)
/// Returns the for stimulated emission corrected photoionisation rate coefficient.
/// Only needed during packet propagation, therefore the value is taken from the
/// cell history if known.
{
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  double T_R = get_TR(modelgridindex);
  double W = get_W(modelgridindex);

  double gammacorr = W * interpolate_corrphotoioncoeff(element,ion,level,phixstargetindex,T_R);
  return gammacorr;
}


///***************************************************************************/
double get_bfheatingcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex)
{
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  double bfheating = -1.;

  if (use_cellhist)
  {
    bfheating = cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].bfheating;
  }

  if (bfheating < 0 || !use_cellhist)
  {

    /*double nnlevel = calculate_exclevelpop(cellnumber,element,ion,level);
    bfheating = nnlevel * W * interpolate_bfheatingcoeff_below(element,ion,level,T_R);
    index_in_groundlevelcontestimator = elements[element].ions[ion].levels[level].closestgroundlevelcont;
    if (index_in_groundlevelcontestimator >= 0) bfheating *= bfheatingestimator[cellnumber*nelements*maxion+index_in_groundlevelcontestimator];*/

  #ifdef NO_LUT_BFHEATING
    bfheating = calculate_bfheatingcoeff(element,ion,level,phixstargetindex,modelgridindex);
  #else
    double T_R = get_TR(modelgridindex);
    double W = get_W(modelgridindex);
    bfheating = W * interpolate_bfheatingcoeff(element,ion,level,phixstargetindex,T_R);
    int index_in_groundlevelcontestimator = elements[element].ions[ion].levels[level].closestgroundlevelcont;
    if (index_in_groundlevelcontestimator >= 0)
      bfheating *= bfheatingestimator[modelgridindex*nelements*maxion + index_in_groundlevelcontestimator];

    if (!isfinite(bfheating))
    {
      printout("[fatal] get_bfheatingcoeff returns a NaN! W %g interpolate_bfheatingcoeff(element,ion,level,phixstargetindex,T_R) %g index_in_groundlevelcontestimator %d bfheatingestimator[modelgridindex*nelements*maxion+index_in_groundlevelcontestimator] %g",
               W,interpolate_bfheatingcoeff(element,ion,level,phixstargetindex,T_R),index_in_groundlevelcontestimator,bfheatingestimator[modelgridindex*nelements*maxion+index_in_groundlevelcontestimator]);
      abort();
    }
  #endif

    if (use_cellhist)
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].bfheating = bfheating;
  }

  return bfheating;
  //return get_bfheatingcoeff_ana(element, ion, level, cellnumber);
}


///***************************************************************************/
double get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex)
{
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  double T_R = get_TR(modelgridindex);
  double W = get_W(modelgridindex);

  /*double nnlevel = calculate_exclevelpop(cellnumber,element,ion,level);
  bfheating = nnlevel * W * interpolate_bfheatingcoeff_below(element,ion,level,T_R);*/
  return W * interpolate_bfheatingcoeff(element,ion,level,phixstargetindex,T_R);
}

#endif /* IFNDEF FORCE_LTE */


///***************************************************************************/
double get_bfcooling(int element, int ion, int level, int phixstargetindex, int modelgridindex)
/// Returns the rate for bfheating. This can be called during packet propagation
/// or update_grid. Therefore we need to decide whether a cell history is
/// known or not.
{
  double bfcooling = -99.;

  if (use_cellhist)
    bfcooling = cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].bfcooling;

  if (bfcooling < 0)
  {
    /// Interpolate bfcooling out of precalculated values
    //int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
    double T_e = get_Te(modelgridindex);
    double nnion = ionstagepop(modelgridindex,element,ion+1);
    double nne = get_nne(modelgridindex);
    //double nnupperlevel = calculate_exclevelpop(modelgridindex,element,ion+1,upper);
    //bfcooling = interpolate_bfcoolingcoeff(element,ion,level,T_e) * nnionlevel * nne;
    bfcooling = interpolate_bfcoolingcoeff(element,ion,level,phixstargetindex,T_e) * nnion * nne;

    if (use_cellhist)
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].bfcooling = bfcooling;

    #ifdef DEBUG_ON
      if (!isfinite(bfcooling))
      {
        printout("[fatal] get_bfcooling: bfcooling infinite (%g) for element %d, ion %d, level %d in modelgridcell %d\n",bfcooling,element,ion,level,modelgridindex);
        //printout("[fatal] get_bfcooling: bfcoolingcoeff %g, nnion %g, nne %g, T_e %g\n",interpolate_bfcoolingcoeff(element,ion,level,phixstargetindex,T_e),nnion,nne,T_e);
      }
    #endif
  }

  return bfcooling;
}


///***************************************************************************/
/*double get_spontrecomb(int element, int ion, int level, int cellnumber)
/// Returns the rate coefficient for spontaneous recombination. Only needed during
/// packet propagation, therefore the value is taken from the
/// cell history if known.
/// For ionisation to other levels than the ground level this must be adapted.
{
  double alpha_sp;

  double T_e = cell[cellnumber].T_e;
  double nne = cell[cellnumber].nne;
  double nnionlevel = get_groundlevelpop(cellnumber,element,ion+1);

  alpha_sp = interpolate_spontrecombcoeff(element,ion,level,T_e) * nnionlevel*nne;

  return alpha_sp;
}*/


double calculate_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex)
{
  double integratorrelaccuracy = 1e-2;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  int upperlevel = get_phixsupperlevel(element,ion,level,phixstargetindex);
  float phixstargetprobability = get_phixsprobability(element,ion,level,phixstargetindex);

  double E_threshold = epsilon(element,ion+1,upperlevel) - epsilon(element,ion,level);
  double nu_threshold = E_threshold / H;
  double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table

  gsl_integral_paras_gammacorr intparas;
  intparas.nu_edge = nu_threshold;
  intparas.modelgridindex = modelgridindex;
  intparas.photoion_xs = elements[element].ions[ion].levels[level].photoion_xs;

  double gammacorr = 0.0;
  gsl_function F_gammacorr;
  F_gammacorr.function = &gammacorr_integrand_gsl_radfield;
  F_gammacorr.params = &intparas;
  double error = 0.0;
  gsl_integration_qag(&F_gammacorr, nu_threshold, nu_max_phixs, 0,
                      integratorrelaccuracy, 1000, 4, w, &gammacorr, &error);
  gammacorr *= FOURPI * phixstargetprobability;

  gsl_integration_workspace_free(w);

  return gammacorr;
}


double gammacorr_integrand_gsl_radfield(double nu, void *paras)
/// Integrand to calculate the rate coefficient for photoionization
/// using gsl integrators. Corrected for stimulated recombination.
{
  gsl_integral_paras_gammacorr localparas = *((gsl_integral_paras_gammacorr*) paras);
  int modelgridindex = localparas.modelgridindex;
  double nu_edge = localparas.nu_edge;
  float *photoion_xs = localparas.photoion_xs;

  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!

  double T_R = get_TR(modelgridindex);

  int i = (int) ((nu/nu_edge - 1.0) / NPHIXSNUINCREMENT);

  #ifdef DEBUG_ON
  /*if (i > NPHIXSPOINTS-1)
  {
    printout("gammacorr_integrand_gsl_radfield called with nu > nu_edge * %g",last_phixs_nuovernuedge);
    abort();
  }*/
  #endif

  //TODO: MK thesis page 41, use population ratios and Te
  return ONEOVERH * photoion_xs[i] / nu * radfield(nu,modelgridindex) *
         (1 - exp(-HOVERKB * nu / T_R));
}


double calculate_bfheatingcoeff(int element, int ion, int level,
                                int phixstargetindex, int modelgridindex)
{
  double error = 0.0;
  double integratoraccuracy = 1e-2;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  int upperlevel = get_phixsupperlevel(element,ion,level,phixstargetindex);
  float phixstargetprobability = get_phixsprobability(element,ion,level,phixstargetindex);
  mastate[tid].level = element;  //passed to photoionization_crosssection
  mastate[tid].level = ion;
  mastate[tid].level = level;

  double E_threshold = epsilon(element,ion+1,upperlevel) - epsilon(element,ion,level);
  double nu_threshold = E_threshold / H;
  double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table

  gsl_integral_paras_gammacorr intparas;
  intparas.nu_edge = nu_threshold;
  intparas.modelgridindex = modelgridindex;

  double bfheating = 0.0;
  gsl_function F_bfheating;
  F_bfheating.function = &bfheating_integrand_gsl_radfield;
  F_bfheating.params = &intparas;
  gsl_integration_qag(&F_bfheating, nu_threshold, nu_max_phixs, 0,
                      integratoraccuracy, 1000, 6, w, &bfheating, &error);
  bfheating *= FOURPI * phixstargetprobability;

  gsl_integration_workspace_free(w);

  return bfheating;
}


double bfheating_integrand_gsl_radfield(double nu, void *paras)
/// Integrand to calculate the rate coefficient for bfheating
/// using gsl integrators.
{
  int modelgridindex = ((gsl_integral_paras_gammacorr *) paras)->modelgridindex;
  double nu_edge = ((gsl_integral_paras_gammacorr *) paras)->nu_edge;

  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection(nu_edge, nu);

  double T_R = get_TR(modelgridindex);

  return sigma_bf * (1-nu_edge/nu) * radfield(nu,modelgridindex) *
         (1 - exp(-HOVERKB*nu/T_R));
}

///***************************************************************************/
/*
void check_interpolation(double T_min, double T_max)
/// Function writes "exact" and interpolated values of the rate coefficients
/// to alpha_sp_file, alpha_st_file and gamma_file to judge the quality
/// of the interpolation. Both integrator types are used.
{
  double interpolate_alpha_sp(int element, int ion, int level, double T_e);
  double interpolate_modified_alpha_sp(int element, int ion, int level, double T_e);
  double interpolate_alpha_st(int element, int ion, int level, double T_e);
  double interpolate_modified_alpha_st(int element, int ion, int level, double T_e);
  double interpolate_gamma(int element, int ion, int level, double T_e);
  double interpolate_modified_gamma(int element, int ion, int level, double T_e);

  double modified_alpha_sp_integrand_gsl(double nu, void * paras);
  double alpha_st_integrand_gsl(double nu, void * paras);
  double modified_alpha_st_integrand_gsl(double nu, void * paras);
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
      sfac = calculate_sahafact(0,0,level,phixstargetindex,T_e,E_threshold);

      /// calculate from tabulated values
      alpha_sp = interpolate_alpha_sp(0,0,level,T_e);
      modified_alpha_sp = interpolate_modified_alpha_sp(0,0,level,T_e);
      alpha_st = interpolate_alpha_st(0,0,level,T_e)*sf;  /// The sahafactor was taken out of the precalculation
      modified_alpha_st = interpolate_modified_alpha_st(0,0,level,T_e)*sf; /// The sahafactor was taken out of the precalculation
      gamma = interpolate_gamma(0,0,level,T_e);
      modified_gamma = interpolate_modified_gamma(0,0,level,T_e);

      /// calculate with gsl integrators
      nu_max_phixs = nu_threshold * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table
      intparas.T = T_e;
      F1.params = &intparas;
      F2.params = &intparas;
      F3.params = &intparas;
      F4.params = &intparas;
      F5.params = &intparas;
      F6.params = &intparas;
      gsl_integration_qng(&F1, nu_threshold, nu_max_phixs, 0, 1e-2, &alpha_sp_gsl, &error, &neval);
      gsl_integration_qng(&F2, nu_threshold, nu_max_phixs, 0, 1e-2, &modified_alpha_sp_gsl, &error, &neval);
      gsl_integration_qng(&F3, nu_threshold, nu_max_phixs, 0, 1e-2, &alpha_st_gsl, &error, &neval);
      gsl_integration_qng(&F4, nu_threshold, nu_max_phixs, 0, 1e-2, &modified_alpha_st_gsl, &error, &neval);
      gsl_integration_qng(&F5, nu_threshold, nu_max_phixs, 0, 1e-2, &gamma_gsl, &error, &neval);
      gsl_integration_qng(&F6, nu_threshold, nu_max_phixs, 0, 1e-2, &modified_gamma_gsl, &error, &neval);
      alpha_sp_gsl *= FOURPI * sf;
      modified_alpha_sp_gsl *= FOURPI * sf;
      alpha_st_gsl *= FOURPI * sf;
      modified_alpha_st_gsl *= FOURPI * sf;
      gamma_gsl *= FOURPI;
      modified_gamma_gsl *= FOURPI;

      /// calculate with qromb integrator of NR
      alpha_sp_integrand_parameters.T = T_e;
      alpha_sp_nr = qromb(alpha_sp_integrand,nu_threshold,nu_max_phixs);
      modified_alpha_sp_nr = qromb(modified_alpha_sp_integrand,nu_threshold,nu_max_phixs);
      alpha_st_nr = qromb(alpha_st_integrand,nu_threshold,nu_max_phixs);
      modified_alpha_st_nr = qromb(modified_alpha_st_integrand,nu_threshold,nu_max_phixs);
      gamma_nr = qromb(gamma_integrand,nu_threshold,nu_max_phixs);
      modified_gamma_nr = qromb(modified_gamma_integrand,nu_threshold,nu_max_phixs);
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
