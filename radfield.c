#include "grid_init.h"
#include "radfield.h"
#include "sn3d.h"
#include <math.h>
#include <stdbool.h>

//#const int RADFIELDBINCOUNT = 1000;
#define RADFIELDBINCOUNT 10

#define nu_lower_first_initial (CLIGHT / (10000e-8))
#define nu_upper_last_initial (CLIGHT / (50e-8))

double nu_lower_first = nu_lower_first_initial;

double J_normfactor[MMODELGRID + 1];

bool radfield_initialized = false;

typedef enum
{
  FIT_DILUTED_BLACKBODY = 0,
  FIT_CONSTANT = 1,
} enum_bin_fit_type;

typedef struct
{
  double nu_upper;   //lower wavelength boundary of this bin
  double J_raw;      //value needs to be multipled by J_normfactor
                     //to get the true value
  double nuJ_raw;
  double W;          // scaling factor
  double T_R;        // radiation temperature
  int contribcount;
  enum_bin_fit_type fit_type;
} radfieldbin;

radfieldbin *radfieldbins[MMODELGRID + 1];

typedef enum
{
  ONE = 0,
  TIMES_NU = 1,
  TIMES_E = 2
} enum_prefactor;

typedef struct
{
  double T_R;
  enum_prefactor prefactor;
} gsl_planck_integral_paras;

typedef struct
{
  int modelgridindex;
  int binindex;
} gsl_T_R_solver_paras;

FILE *radfieldfile = NULL;

// private functions
double find_T_R(int modelgridindex, int binindex);
double delta_nu_bar(double T_R, void *paras);
double integrate_planck(double T_R, double nu_lower, double nu_upper, enum_prefactor prefactor, double *error);
double gsl_integrand_planck(double nu, void *paras);
double radfield_get_bin_J(int modelgridindex, int binindex);
int radfield_get_bin_contribcount(int modelgridindex, int binindex);
double radfield_get_bin_W(int modelgridindex, int binindex);
double radfield_get_bin_T_R(int modelgridindex, int binindex);
int radfield_select_bin(int modelgridindex, double nu);
double radfield_get_bin_nu_lower(int modelgridindex, int binindex);
double radfield_get_bin_nu_upper(int modelgridindex, int binindex);

void radfield_init(void)
{
  if (radfield_initialized == false)
  {
    char filename[100] = "radfield.out";
    radfieldfile = fopen(filename, "w");
    if (radfieldfile == NULL)
    {
      printout("Cannot open %s.\n",filename);
      exit(0);
    }
    fprintf(radfieldfile,"%8s %15s %8s %11s %11s %9s %9s %9s %9s %9s %11s\n",
            "timestep","modelgridindex","bin_num","nu_lower","nu_upper",
            "nuJ","J","J_nu_avg","ncontrib","T_R","W");
    fflush(radfieldfile);

    for (int modelgridindex = 0; modelgridindex < MMODELGRID + 1; modelgridindex++)
    {
      radfieldbins[modelgridindex] = (radfieldbin *) calloc(RADFIELDBINCOUNT, sizeof(radfieldbin));

      radfield_set_J_normfactor(modelgridindex, -1.0);

      double prev_nu_upper = nu_lower_first_initial;
      //double delta_nu = (nu_upper_last_initial - nu_lower_first_initial) / RADFIELDBINCOUNT; // upper limit if no edges are crossed
      double delta_lambda = ((1 / nu_lower_first_initial) - (1 / nu_upper_last_initial)) / RADFIELDBINCOUNT; // upper limit if no edges are crossed

      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      {
        double delta_nu = pow(prev_nu_upper,2) * delta_lambda; // equally spaced in wavelength
        radfieldbins[modelgridindex][binindex].nu_upper = prev_nu_upper + delta_nu;
        prev_nu_upper = radfieldbins[modelgridindex][binindex].nu_upper; // importantly, the below part doesn't change this

        // Align the bin edges with bound-free edges
        if (binindex < RADFIELDBINCOUNT - 1)
        {
          /*for (int i = 0; i < nbfcontinua_ground; i++)
          {
            double nu_edge = phixslist[tid].groundcont[i].nu_edge;
            //printout("bf edge at %g, nu_lower_first %g, nu_upper_last %g\n",nu_edge,nu_lower_first,nu_upper_last);
            if ((nu_edge < nu_lower_first) || (nu_edge > nu_upper_last))
            {
              printout("MISSED bf edge at %g, nu_lower_first %g, nu_upper_last %g\n",nu_edge,nu_lower_first,nu_upper_last);
            }
            if ((nu_edge > prev_nu_upper) &&
                (nu_edge < radfieldbins[modelgridindex][binindex].nu_upper))
              radfieldbins[modelgridindex][binindex].nu_upper = nu_edge;
          }*/
        }
        else
        {
          radfieldbins[modelgridindex][binindex].nu_upper = nu_upper_last_initial;
        }

        radfieldbins[modelgridindex][binindex].W = -1.0;
        radfieldbins[modelgridindex][binindex].T_R = -1.0;
        //radfieldbins[modelgridindex][binindex].fit_type = FIT_CONSTANT;
        radfieldbins[modelgridindex][binindex].fit_type = FIT_DILUTED_BLACKBODY;
      }
    }
    radfield_initialized = true;
  }
}


void radfield_write_to_file(int modelgridindex, int timestep)
{
  if (!radfield_initialized)
  {
    printout("Call to write_to_radfield_file before init_radfield");
    abort();
  }

  int totalcontribs = 0;
  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    totalcontribs += radfield_get_bin_contribcount(modelgridindex, binindex);

  for (int binindex = -1; binindex < RADFIELDBINCOUNT; binindex++)
  {
    double nu_lower = 0.0;
    double nu_upper = 0.0;
    double nuJ_out = 0.0;
    double J_out = 0.0;
    double T_R = 0.0;
    double W = 0.0;
    double J_nu_bar = 0.0;
    int contribcount = 0;

    if (binindex >= 0)
    {
      nu_lower = radfield_get_bin_nu_lower(modelgridindex, binindex);
      nu_upper = radfieldbins[modelgridindex][binindex].nu_upper;
      nuJ_out = radfieldbins[modelgridindex][binindex].nuJ_raw *
                   J_normfactor[modelgridindex];
      J_out = radfield_get_bin_J(modelgridindex, binindex);
      T_R = radfield_get_bin_T_R(modelgridindex, binindex);
      W = radfield_get_bin_W(modelgridindex, binindex);
      J_nu_bar = J_out / (nu_upper - nu_lower);
      contribcount = radfield_get_bin_contribcount(modelgridindex, binindex);
    }
    else
    {
      nuJ_out = nuJ[modelgridindex];
      J_out = J[modelgridindex];
      T_R = get_TR(modelgridindex);
      W = get_W(modelgridindex);
      contribcount = totalcontribs;
    }

    fprintf(radfieldfile,"%8d %15d %8d %11.5e %11.5e %9.3e %9.3e %9.3e %9d %9.1f %11.5e\n",
            timestep,modelgridindex,binindex,nu_lower,nu_upper,nuJ_out,J_out,J_nu_bar,contribcount,T_R,W);
  }
  fflush(radfieldfile);
}


void radfield_close_file(void)
{
  fclose(radfieldfile);

  for (int dmgi = 0; dmgi < MMODELGRID + 1; dmgi++)
    free(radfieldbins[dmgi]);

  //free(radfieldbins);
}

// set up the new bins and clear the estimators in preparation
// for a timestep
void radfield_zero_estimators(int modelgridindex)
{
  nuJ[modelgridindex] = 0.;

  if (!radfield_initialized)
    radfield_init();

  printout("radfield: zeroing estimators in %d bins in cell %d...",RADFIELDBINCOUNT,modelgridindex);

  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
  {
    radfieldbins[modelgridindex][binindex].J_raw = 0.0;
    radfieldbins[modelgridindex][binindex].nuJ_raw = 0.0;
    radfieldbins[modelgridindex][binindex].contribcount = 0;
  }
  radfield_set_J_normfactor(modelgridindex, -1.0);
  printout("done.\n");
}

inline
void radfield_update_estimators(int modelgridindex, double distance,
                                double e_cmf, double nu_cmf)
{
  int binindex = 0;
  if (nu_cmf <= radfield_get_bin_nu_lower(modelgridindex,binindex))
  {
    #ifdef DEBUG_ON
    printout("radfield: Extending nu_lower_first from %g down to %g\n",nu_lower_first,nu_cmf);
    #endif
    nu_lower_first = nu_cmf;
  }
  else if (nu_cmf > radfieldbins[modelgridindex][RADFIELDBINCOUNT - 1].nu_upper)
  {
    binindex = RADFIELDBINCOUNT - 1;
    #ifdef DEBUG_ON
    printout("radfield: Extending nu_upper_last from %g up to %g\n",radfield_get_bin_nu_upper(modelgridindex,binindex),nu_cmf);
    #endif
    radfieldbins[modelgridindex][binindex].nu_upper = nu_cmf;
  }
  else
  {
    binindex = radfield_select_bin(modelgridindex,nu_cmf);
  }

  if (binindex >= 0)
  {
    #ifdef _OPENMP
      #pragma omp critical
    #endif
    {
      radfieldbins[modelgridindex][binindex].J_raw += distance * e_cmf;
      radfieldbins[modelgridindex][binindex].nuJ_raw += distance * e_cmf *
                                                        nu_cmf;
      radfieldbins[modelgridindex][binindex].contribcount += 1;
    }
  }
  else
  {
    // dropping the contribution of this packet
    printout("WARNING: radfield_update_estimators dropping packet contribution for nu_cmf %g\n",
             nu_cmf);
    printout("           modelgridindex %d binindex %d nu_lower_first %g nu_upper_last %g \n",
             modelgridindex, binindex, nu_lower_first, radfield_get_bin_nu_upper(modelgridindex,RADFIELDBINCOUNT - 1));
  }
}


double radfield(double nu, int modelgridindex)
{
  double T_R_fullspec = get_TR(modelgridindex);
  double W_fullspec   = get_W(modelgridindex);

  #ifdef USE_MULTIBIN_RADFIELD_MODEL
  if (radfield_initialized) // && radfieldbins[modelgridindex] != NULL
  {
    int binindex = radfield_select_bin(modelgridindex,nu);
    if (binindex >= 0)
    {
      double W_bin = radfieldbins[modelgridindex][binindex].W;
      if (W_bin >= 0.)
      {
        if (radfieldbins[modelgridindex][binindex].fit_type == FIT_CONSTANT)
        {
            return W_bin;
        }
        else
        {
          double T_R_bin = radfieldbins[modelgridindex][binindex].T_R;
          if (T_R_bin > 0.)
          {
            double J_nu = radfield2(nu, T_R_bin, W_bin);
            /*if (fabs(J_nu / J_nu_fullspec - 1.0) > 0.5)
            {
              printout("WARNING: radfield: significant discrepancy. J_nu_fullspec %g, J_nu %g, nu %g W_bin %g T_R_bin %g\n",
                       J_nu_fullspec, J_nu, nu, W_bin, T_R_bin);
            }*/
            return J_nu;
          }
          //else
          //{
          //  printout("WARNING: Radfield modelgridindex %d binindex %d has W %g T_R=%g<=0, using W %g T_R %g nu %g\n",
          //           modelgridindex, binindex, W_bin, T_R_bin, W_fullspec, T_R_fullspec, nu);
          //}
        }
      }
      else
      {
        //printout("WARNING: Radfield modelgridindex %d binindex %d has W_bin=%g<0, using W %g T_R %g nu %g\n",
        //         modelgridindex, binindex, W_bin, W_fullspec, T_R_fullspec, nu);
      }
    }
    else
    {
      //printout("WARNING: Radfield modelgridindex %d binindex %d nu %g nu_lower_first %g nu_upper_last %g \n",
      //         modelgridindex, binindex, nu, nu_lower_first, nu_upper_last);
    }
  }
  /*else
  {
    printout("radfield: WARNING: Radfield called before initialized. Using global T_R %g W %g nu %g modelgridindex %d\n",
             W_fullspec, T_R_fullspec, nu, modelgridindex);
  }*/
  #endif

  double J_nu_fullspec = radfield2(nu, T_R_fullspec, W_fullspec);
  return J_nu_fullspec;
}


void radfield_fit_parameters(int modelgridindex)
// finds the best fitting W and temperature parameters in each spectral bin
// using J and nuJ
{
  double T_J, T_R, W;
  get_radfield_params_fullspec(J[modelgridindex],nuJ[modelgridindex],modelgridindex,&T_J,&T_R,&W);
  modelgrid[modelgridindex].TJ = T_J;
  modelgrid[modelgridindex].TR = T_R;
  modelgrid[modelgridindex].W = W;

  if (J_normfactor[modelgridindex] <= 0)
  {
    printout("radfield: FATAL J_normfactor = %g in cell %d at call to radfield_fit_parameters",J_normfactor[modelgridindex],modelgridindex);
    abort();
  }

  float T_R_fullspec = get_TR(modelgridindex);
  double J_fullspec = J[modelgridindex];
  //double plank_integral_zero_inf = STEBO * pow(T_R_fullspec,4) / PI;

  printout("Full-spectrum fit radfield params for mgi %d: J %g, T_R %g, W %g\n",
           modelgridindex, J_fullspec, T_R_fullspec, get_W(modelgridindex));

  printout("radfield: Finding parameters for %d bins...\n",RADFIELDBINCOUNT);

  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
  {
    double J_bin = radfield_get_bin_J(modelgridindex,binindex);
    double nu_lower = radfield_get_bin_nu_lower(modelgridindex, binindex);
    double nu_upper = radfieldbins[modelgridindex][binindex].nu_upper;
    double T_R_bin = -1.0;
    double W_bin = -1.0;
    if (radfieldbins[modelgridindex][binindex].contribcount > 10)
    {
      if (radfieldbins[modelgridindex][binindex].fit_type == FIT_DILUTED_BLACKBODY)
      {
        T_R_bin = find_T_R(modelgridindex, binindex);

        double integerror = 0.0;
        double plank_integral = integrate_planck(T_R_bin, nu_lower,
                                                 nu_upper, ONE, &integerror);

        W_bin = J_bin / plank_integral;

        if (W_bin > 1e2)
        {
          printout("W %g too high, try setting T_R of bin %d to %g. J_bin %g planck_integral %g +/- %g\n",W_bin,binindex,MAXTEMP,plank_integral,integerror);
          double plank_integral = integrate_planck(MAXTEMP, nu_lower,
                                                   nu_upper, ONE, &integerror);
          W_bin = J_bin / plank_integral;
          if (W_bin > 1e2)
          {
            printout("W still too high %g\n",W_bin);
            T_R_bin = -99.0;
            W_bin = 0.;
            abort();
          }
        }
      }
      else //if (radfieldbins[modelgridindex][binindex].fit_type == FIT_CONSTANT)
      {
        double delta_nu = nu_upper - nu_lower;

        T_R_bin = -1;
        W_bin = J_bin / delta_nu;
      }
    }
    radfieldbins[modelgridindex][binindex].T_R = T_R_bin;
    radfieldbins[modelgridindex][binindex].W = W_bin;
  }

  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
  {
    double J_bin = radfield_get_bin_J(modelgridindex,binindex);
    double T_R_bin = radfield_get_bin_T_R(modelgridindex,binindex);
    double W_bin = radfield_get_bin_W(modelgridindex,binindex);
    printout("bin %d: J %g, T_R %g, W %g\n",
           binindex, J_bin, T_R_bin, W_bin);
  }
}


void get_radfield_params_fullspec(double J, double nuJ, int modelgridindex, double *T_J, double *T_R, double *W)
{
  double nubar = nuJ/J;
  if (!isfinite(nubar) || nubar == 0.)
  {
    /// Return old T_R
    printout("[warning] update_grid: T_R estimator infinite in cell %d, use value of last timestep\n",modelgridindex);
    *T_J = modelgrid[modelgridindex].TJ;
    *T_R = modelgrid[modelgridindex].TR;
    *W = modelgrid[modelgridindex].W;
  }
  else
  {
    *T_J = pow(PI/STEBO*J,1./4.);
    if (*T_J > MAXTEMP)
    {
      printout("[warning] update_grid: temperature estimator T_J=%g exceeds T_max=%g in cell %d. Set T_J = T_max!\n",*T_J,MAXTEMP,modelgridindex);
      *T_J = MAXTEMP;
    }
    if (*T_J < MINTEMP)
    {
      printout("[warning] update_grid: temperature estimator T_J=%g below T_min %g in cell %d. Set T_J = T_min!\n",*T_J,MINTEMP,modelgridindex);
      *T_J = MINTEMP;
    }

    *T_R = H*nubar/KB/3.832229494;
    if (*T_R > MAXTEMP)
    {
      printout("[warning] update_grid: temperature estimator T_R=%g exceeds T_max=%g in cell %d. Set T_R = T_max!\n",*T_R,MAXTEMP,modelgridindex);
      *T_R = MAXTEMP;
    }
    if (*T_R < MINTEMP)
    {
      printout("[warning] update_grid: temperature estimator T_R=%g below T_min %g in cell %d. Set T_R = T_min!\n",*T_R,MINTEMP,modelgridindex);
      *T_R = MINTEMP;
    }

    *W = PI*J/STEBO/pow(*T_R,4);
  }
}


double find_T_R(int modelgridindex, int binindex)
{
  double fractional_accuracy = 1e-4;
  int maxit = 100;
  double T_R = 0.0;

  double T_R_min = MINTEMP;
  double T_R_max = MAXTEMP * 4;

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
  gsl_root_fsolver *T_R_solver;
  solvertype = gsl_root_fsolver_brent;

  gsl_T_R_solver_paras paras;
  paras.modelgridindex = modelgridindex;
  paras.binindex = binindex;

  gsl_function find_T_R_f;
  find_T_R_f.function = &delta_nu_bar;
  find_T_R_f.params = &paras;

  /// Check whether the equation has a root in [T_min,T_max]
  double delta_nu_bar_min = delta_nu_bar(T_R_min,find_T_R_f.params);
  double delta_nu_bar_max = delta_nu_bar(T_R_max,find_T_R_f.params);

  //printout(
  //  "call_T_R_finder: bin %d delta_nu_bar(T_R_min) %g, delta_nu_bar(T_R_max) %g\n",
  //  binindex, delta_nu_bar_min,delta_nu_bar_max);

  if (!isfinite(delta_nu_bar_min) || !isfinite(delta_nu_bar_max))
    delta_nu_bar_max = delta_nu_bar_min = -1;

  if (delta_nu_bar_min * delta_nu_bar_max < 0)
  {
    /// If it has, then solve for the root T_R

    /// now start so iterate on solution
    T_R_solver = gsl_root_fsolver_alloc(solvertype);
    gsl_root_fsolver_set(T_R_solver, &find_T_R_f, T_R_min, T_R_max);
    int iteration_num = 0;
    int status = GSL_CONTINUE;
    do
    {
      iteration_num++;
      status = gsl_root_fsolver_iterate(T_R_solver);
      T_R = gsl_root_fsolver_root(T_R_solver);
      //cell[cellnumber].T_e = T_e;

      T_R_min = gsl_root_fsolver_x_lower(T_R_solver);
      T_R_max = gsl_root_fsolver_x_upper(T_R_solver);
      status = gsl_root_test_interval(T_R_min,T_R_max,0,fractional_accuracy);
      printout(
        "call_T_R_finder: bin %4d iter %d, T_R is between %7.1f and %7.1f, guess %7.1f, delta_nu_bar %g, status %d\n",
        binindex,iteration_num,T_R_min,T_R_max,T_R,delta_nu_bar(T_R,&paras),
        status);
    }
    while (status == GSL_CONTINUE && iteration_num < maxit);

    if (status == GSL_CONTINUE)
      printout(
        "[warning] call_T_R_finder: T_R did not converge within %d iterations\n",
        maxit);

    gsl_root_fsolver_free(T_R_solver);
  }
  /// Quick solver style: works if we can assume that there is either one or no
  /// solution on [MINTEM.MAXTEMP] (check that by doing a plot of heating-cooling
  /// vs. T_e using the tb_info switch)
  else if (delta_nu_bar_max < 0)
  {
    /// Thermal balance equation always negative ===> T_R = T_min
    /// Calculate the rates again at this T_e to print them to file
    T_R = MINTEMP;
    printout("call_T_R_finder: bin %4d no solution in interval, clamping to MINTEMP\n",
             binindex);
  }
  else
  {
    T_R = MAXTEMP;
    printout("call_T_R_finder: bin %4d no solution in interval, clamping to MAXTEMP\n",
             binindex);
  }

  return T_R;
}


double delta_nu_bar(double T_R, void *paras)
// difference between the average nu and the average nu of a planck function
// at temperature T_R, in the frequency range corresponding to a bin
{
  int modelgridindex = ((gsl_T_R_solver_paras *) paras)->modelgridindex;
  int binindex = ((gsl_T_R_solver_paras *) paras)->binindex;

  double nu_lower = radfield_get_bin_nu_lower(modelgridindex, binindex);
  double nu_upper = radfieldbins[modelgridindex][binindex].nu_upper;

  double nuJ_bin_raw = radfieldbins[modelgridindex][binindex].nuJ_raw;
  double J_bin_raw = radfieldbins[modelgridindex][binindex].J_raw;
  double nu_bar = nuJ_bin_raw / J_bin_raw;

  double integerror = 0.0;
  double nu_times_plank_integral = integrate_planck(T_R, nu_lower,
                                                    nu_upper, TIMES_NU, &integerror);

  double plank_integral = integrate_planck(T_R, nu_lower, nu_upper, ONE, &integerror);

  double nu_bar_plank = nu_times_plank_integral / plank_integral;

  //printout("nu_bar %g nu_bar_plank(T=%g) %g\n",nu_bar,T_R,nu_bar_plank);

  //double delta_nu_bar = nu_bar_plank - nu_bar;
  double delta_nu_bar = nu_bar_plank / nu_bar - 1.0;

  return delta_nu_bar;
}


double integrate_planck(double T_R, double nu_lower, double nu_upper, enum_prefactor prefactor, double *error)
{
  //double error = 0.0;
  double integratoraccuracy = 1e-10;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);

  gsl_planck_integral_paras intparas;
  intparas.T_R = T_R;
  intparas.prefactor = prefactor;

  double integral = 0.0;
  gsl_function F_plank;
  F_plank.function = &gsl_integrand_planck;
  F_plank.params = &intparas;
  gsl_integration_qag(&F_plank, nu_lower, nu_upper, 0.,
                      integratoraccuracy, 100000, 6, w, &integral,
                      error);

  gsl_integration_workspace_free(w);

  return integral;
}


double gsl_integrand_planck(double nu, void *paras)
{
  double T_R = ((gsl_planck_integral_paras *) paras)->T_R;
  enum_prefactor prefactor = ((gsl_planck_integral_paras *) paras)->prefactor;

  double integrand = TWOHOVERCLIGHTSQUARED * pow(nu,3) /
                     (expm1(HOVERKB * nu / T_R));

  if (prefactor == TIMES_NU)
    integrand *= nu;
  else if (prefactor == TIMES_E)
    integrand *= H * nu;

  return integrand;
}


void radfield_set_J_normfactor(int modelgridindex, double normfactor)
{
  J_normfactor[modelgridindex] = normfactor;
}


double radfield_get_bin_J(int modelgridindex, int binindex)
{
  if (J_normfactor[modelgridindex] <= 0.0)
  {
    printout(
      "radfield: Fatal error: radfield_get_bin_J called before J_normfactor set for modelgridindex %d",
      modelgridindex);
    abort();
  }
  return radfieldbins[modelgridindex][binindex].J_raw *
         J_normfactor[modelgridindex];
}


inline
double radfield_get_bin_W(int modelgridindex, int binindex)
{
  return radfieldbins[modelgridindex][binindex].W;
}


inline
double radfield_get_bin_T_R(int modelgridindex, int binindex)
{
  return radfieldbins[modelgridindex][binindex].T_R;
}


inline
int radfield_select_bin(int modelgridindex, double nu)
{
  // the slow way to search a list, one by one until found
  /*
  for (int dbinindex = 0; dbinindex < RADFIELDBINCOUNT; dbinindex++)
  {
    if (radfield_get_bin_nu_lower(modelgridindex,dbinindex) < nu &&
        radfieldbins[modelgridindex][dbinindex].nu_upper >= nu)
    {
      return dbinindex;
    }
  }
  return -1;
  */

  // binary search
  int low = 0;
  int high = RADFIELDBINCOUNT - 1;
  while (low <= high)
  {
    int mid = low + ((high - low) / 2);
    if (radfieldbins[modelgridindex][mid].nu_upper <= nu)
    {
      low = mid + 1;
    }
    else if (radfield_get_bin_nu_lower(modelgridindex, mid) > nu)
    {
      high = mid - 1;
    }
    else
    {
      return mid;
    }
   }
   return -1;
}


inline
double radfield_get_bin_nu_lower(int modelgridindex, int binindex)
{
  if (binindex > 0)
    return radfieldbins[modelgridindex][binindex - 1].nu_upper;
  else
    return nu_lower_first;
}


inline
double radfield_get_bin_nu_upper(int modelgridindex, int binindex)
{
  return radfieldbins[modelgridindex][binindex].nu_upper;
}


inline
int radfield_get_bin_contribcount(int modelgridindex, int binindex)
{
  return radfieldbins[modelgridindex][binindex].contribcount;
}
