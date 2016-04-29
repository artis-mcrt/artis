#include "grid_init.h"
#include "radfield.h"
#include "sn3d.h"
#include <math.h>
#include <stdbool.h>

#define RADFIELDBINCOUNT 10

double nu_lower_first = KB * MINTEMP / H / 2;
double nu_upper_last = KB * MAXTEMP / H;

double J_normfactor[MMODELGRID + 1];

typedef struct
{
  double nu_upper; //lower wavelength boundary of this bin
  double J_raw;        //value needs to be multipled by J_normfactor
                   //to get the true value
  double nuJ_raw;
  double W;
  double T_R;
} radfieldbin;

radfieldbin radfieldbins[MMODELGRID + 1][RADFIELDBINCOUNT];

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
double calculate_planck_integral(double T_R, double nu_lower, double nu_upper,
                                 enum_prefactor prefactor);
double gsl_integrand_planck(double nu, void *paras);


void init_radfield_file(void)
{
  char filename[100] = "radfield.out";
  radfieldfile = fopen(filename, "w");
  if (radfieldfile == NULL)
  {
    printout("Cannot open %s.\n",filename);
    exit(0);
  }
  fprintf(radfieldfile,"%8s %15s %8s %9s %9s %9s %9s %9s %9s \n",
          "timestep","modelgridindex","bin_num","nu_lower","nu_upper",
          "nuJ","J","T_R","W");
  fflush(radfieldfile);
}


void write_to_radfield_file(int modelgridindex, int timestep)
{
  if (radfieldfile == NULL)
  {
    printout("Call to write_to_radfield_file before init_radfield_file");
    abort();
  }
  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
  {
    double nu_lower = get_bin_nu_lower(modelgridindex, binindex);
    double nu_upper = radfieldbins[modelgridindex][binindex].nu_upper;
    double nuJ = radfieldbins[modelgridindex][binindex].nuJ_raw * J_normfactor[modelgridindex];
    double J = radfieldbins[modelgridindex][binindex].J_raw * J_normfactor[modelgridindex];
    double T_R = radfieldbins[modelgridindex][binindex].T_R;
    double W = radfieldbins[modelgridindex][binindex].W;

    fprintf(radfieldfile,"%8d %15d %8d %9.3e %9.3e %9.3e %9.3e %9.1f %9.3e \n",
            timestep,modelgridindex,binindex,nu_lower,nu_upper,nuJ,J,T_R,W);
  }
  fflush(radfieldfile);
}


void close_radfield_file(void)
{
  fclose(radfieldfile);
}


void zero_radfield_estimators(int modelgridindex)
{
  double delta_nu = (nu_upper_last - nu_lower_first) / RADFIELDBINCOUNT;

  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
  {
    radfieldbins[modelgridindex][binindex].nu_upper = nu_lower_first +
                                                      delta_nu +
                                                      (delta_nu * binindex);

    radfieldbins[modelgridindex][binindex].J_raw = 0.0;
    radfieldbins[modelgridindex][binindex].nuJ_raw = 0.0;
    radfieldbins[modelgridindex][binindex].W = 0.0;
    radfieldbins[modelgridindex][binindex].T_R = 0.0;
  }

  J_normfactor[modelgridindex] = -1.0;
}


void update_radfield(int modelgridindex, double distance, double e_cmf,
                     double nu_cmf)
{
  int binindex = get_frequency_bin(modelgridindex,nu_cmf);

  if (binindex >= 0)
  {
    radfieldbins[modelgridindex][binindex].J_raw += distance * e_cmf;
    radfieldbins[modelgridindex][binindex].nuJ_raw += distance * e_cmf * nu_cmf;
  }
  //else, dropping the contribution of this packet
}


double radfield(double nu, int modelgridindex)
{
  float T_R = get_TR(modelgridindex);
  float W   = get_W(modelgridindex);

  return W * TWOHOVERCLIGHTSQUARED *
         pow(nu,3) * 1.0 / (expm1(HOVERKB * nu / T_R));
}


void radfield_fit_parameters(int modelgridindex)
// finds the best fitting W and temperature parameters in each spectral bin
// using J and nuJ
{
  float T_R_old = get_TR(modelgridindex);
  double J_old = J[modelgridindex];
  //double plank_integral_zero_inf = STEBO * pow(T_R_old,4) / PI;

  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
  {
    double nu_lower = get_bin_nu_lower(modelgridindex, binindex);
    double nu_upper = radfieldbins[modelgridindex][binindex].nu_upper;

    double J_bin = radfieldbins[modelgridindex][binindex].J_raw * J_normfactor[modelgridindex];

    double T_R_bin = find_T_R(modelgridindex, binindex);

    radfieldbins[modelgridindex][binindex].T_R = T_R_bin;

    double plank_integral = calculate_planck_integral(T_R_bin, nu_lower, nu_upper,
                                                      ONE);

    double W_bin = J_bin / plank_integral;

    radfieldbins[modelgridindex][binindex].W = W_bin;
  }

  printout("T_R SOLVER RESULTS:\n");
  printout("old full-spectrum fit: J %g, T_R %g, W %g\n",
           J_old, T_R_old, get_W(modelgridindex));
  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
  {
    double J_bin = radfield_get_bin_J(modelgridindex,binindex);
    double T_R_bin = radfieldbins[modelgridindex][binindex].T_R;
    double W_bin = radfieldbins[modelgridindex][binindex].W;
    printout("bin %d: J %g, T_R %g, W %g\n",
             binindex, J_bin, T_R_bin, W_bin);
  }
}


double find_T_R(int modelgridindex, int binindex)
{
  double fractional_accuracy = 1e-4;
  int maxit = 100;
  int iteration_num,status;
  double T_R = 0.0;

  double T_R_min = MINTEMP;
  double T_R_max = MAXTEMP*4;

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

  printout(
    "call_T_R_finder: delta_nu_bar(T_R_min) %g, delta_nu_bar(T_R_max) %g\n",
    delta_nu_bar_min,delta_nu_bar_max);

  if (!isfinite(delta_nu_bar_min) || !isfinite(delta_nu_bar_max))
    delta_nu_bar_max = delta_nu_bar_min = -1;

  if (delta_nu_bar_min * delta_nu_bar_max < 0)
  {
    /// If it has, then solve for the root T_R

    /// now start so iterate on solution
    T_R_solver = gsl_root_fsolver_alloc(solvertype);
    gsl_root_fsolver_set(T_R_solver, &find_T_R_f, T_R_min, T_R_max);
    iteration_num = 0;
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
        "call_T_R_finder: bin %d iter %d, T_R is between %g and %g, guess %g, delta_nu_bar %g, status %d\n",
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
    /// Thermal balance equation always negative ===> T_R = T_min
    /// Calculate the rates again at this T_e to print them to file
    T_R = MINTEMP;
  else
    T_R = MAXTEMP;

  return T_R;
}


double delta_nu_bar(double T_R, void *paras)
// difference between the average nu and the average nu of a planck function
// at temperature T_R, in the frequency range corresponding to a bin
{
  int modelgridindex = ((gsl_T_R_solver_paras *) paras)->modelgridindex;
  int binindex = ((gsl_T_R_solver_paras *) paras)->binindex;

  double nu_lower = nu_lower_first;
  if (binindex > 0)
    nu_lower = radfieldbins[modelgridindex][binindex - 1].nu_upper;

  double nu_upper = radfieldbins[modelgridindex][binindex].nu_upper;

  double nuJ_bin_raw = radfieldbins[modelgridindex][binindex].nuJ_raw;
  double J_bin_raw = radfieldbins[modelgridindex][binindex].J_raw;
  double nu_bar = nuJ_bin_raw / J_bin_raw;

  double nu_times_plank_integral = calculate_planck_integral(T_R, nu_lower,
                                                             nu_upper,
                                                             TIMES_NU);

  double plank_integral = calculate_planck_integral(T_R, nu_lower, nu_upper,
                                                    ONE);

  double nu_bar_plank = nu_times_plank_integral / plank_integral;

  printout("nu_bar %g nu_bar_plank(T=%g) %g\n",nu_bar,T_R,nu_bar_plank);

  //double delta_nu_bar = nu_bar_plank - nu_bar;
  double delta_nu_bar = nu_bar_plank / nu_bar - 1.0;

  return delta_nu_bar;
}


double calculate_planck_integral(double T_R, double nu_lower, double nu_upper,
                                 enum_prefactor prefactor)
{
  double error = 0.0;
  double integratoraccuracy = 1e-10;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);

  gsl_planck_integral_paras intparas;
  intparas.T_R = T_R;
  intparas.prefactor = prefactor;

  double integral = 0.0;
  gsl_function F_plank;
  F_plank.function = &gsl_integrand_planck;
  F_plank.params = &intparas;
  gsl_integration_qag(&F_plank, nu_lower, nu_upper, 0,
                      integratoraccuracy, 100000, 6, w, &integral,
                      &error);

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
  if (J_normfactor[modelgridindex] < 0.0)
  {
    printout("Error: radfield_get_bin_J called before J_normfactor set for modelgridindex %d",
             modelgridindex);
    abort();
  }
  return radfieldbins[modelgridindex][binindex].J_raw * J_normfactor[modelgridindex];
}


int get_frequency_bin(int modelgridindex, double nu)
{
  int binindex = -1;

  for (int dbinindex = 0; dbinindex < RADFIELDBINCOUNT; dbinindex++)
  {
    if (get_bin_nu_lower(modelgridindex,dbinindex) <= nu &&
        radfieldbins[modelgridindex][dbinindex].nu_upper >= nu)
    {
      binindex = dbinindex;
      break;
    }
  }

  return binindex;
}


double get_bin_nu_lower(int modelgridindex, int binindex)
{
  if (binindex > 0)
    return radfieldbins[modelgridindex][binindex-1].nu_upper;
  else
    return nu_lower_first;
}
