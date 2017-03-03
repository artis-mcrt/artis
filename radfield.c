#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_debye.h>
#include <gsl/gsl_roots.h>
#include "atomic.h"
#include "grid_init.h"
#include "radfield.h"
#include "sn3d.h"

#define RADFIELDBINCOUNT 96
// static const int RADFIELDBINCOUNT = 96;

extern inline double radfield_dbb(double nu, float T, float W);

static const double nu_lower_first_initial = (CLIGHT / (20000e-8)); // in Angstroms
static const double nu_upper_last_initial = (CLIGHT /  (2000e-8));  // in Angstroms

static const double boost_region_nu_lower = (CLIGHT / (2500e-8)); // in Angstroms
static const double boost_region_nu_upper = (CLIGHT / (2100e-8));  // in Angstroms
static const double boost_region_factor = 1.0;
static const bool boost_region_on = false;

static double J_normfactor[MMODELGRID + 1];

static bool radfield_initialized = false;

static const double T_R_min = MINTEMP;
static const double T_R_max = MAXTEMP * 2;

// typedef enum
// {
//   FIT_DILUTED_BLACKBODY = 0,
//   FIT_CONSTANT = 1,
// } enum_bin_fit_type;

struct radfieldbin_previous
{
  double prev_J_normed;
  double prev_nuJ_normed;
  int prev_contribcount;
};

struct radfieldbin_current
{
  double J_raw;           // value needs to be multipled by J_normfactor to get the true value
  double nuJ_raw;
  int contribcount;
  float W;                // dilution (scaling) factor
  float T_R;              // radiation temperature
  // enum_bin_fit_type fit_type;
};

static double radfieldbin_nu_upper[RADFIELDBINCOUNT]; // array of upper frequency boundaries of bins
static struct radfieldbin_current *radfieldbin_current[MMODELGRID + 1];
static struct radfieldbin_previous *radfieldbin_previous[MMODELGRID + 1];

static double J[MMODELGRID + 1];
#ifdef DO_TITER
  static double J_reduced_save[MMODELGRID + 1];
#endif

#ifndef FORCE_LTE
  static double nuJ[MMODELGRID + 1];
  #ifdef DO_TITER
    static double nuJ_reduced_save[MMODELGRID + 1];
  #endif
#endif


typedef enum
{
  ONE = 0,
  TIMES_NU = 1,
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

static FILE *restrict radfieldfile = NULL;


static inline
double radfield_get_bin_nu_upper(int binindex)
{
  return radfieldbin_nu_upper[binindex];
}


void radfield_init(int my_rank)
{
  if (!MULTIBIN_RADFIELD_MODEL_ON)
    return;
  if (radfield_initialized)
    printout("WARNING: Tried to initialize radfield twice!\n");

  printout("Initialising radiation field with %d bins from (%6.2f eV, %6.1f A) to (%6.2f eV, %6.1f A)\n",
           RADFIELDBINCOUNT, H * nu_lower_first_initial / EV, 1e8 * CLIGHT / nu_lower_first_initial,
           H * nu_upper_last_initial / EV, 1e8 * CLIGHT / nu_upper_last_initial);
  char filename[100];
  sprintf(filename,"radfield_%.4d.out", my_rank);
  radfieldfile = fopen(filename, "w");
  if (radfieldfile == NULL)
  {
    printout("Cannot open %s.\n",filename);
    abort();
  }
  fprintf(radfieldfile,"%8s %15s %8s %11s %11s %9s %9s %9s %9s %9s %12s\n",
          "timestep","modelgridindex","bin_num","nu_lower","nu_upper",
          "nuJ","J","J_nu_avg","ncontrib","T_R","W");
  fflush(radfieldfile);

  double prev_nu_upper = nu_lower_first_initial;
  const double delta_nu = (nu_upper_last_initial - nu_lower_first_initial) / RADFIELDBINCOUNT; // upper limit if no edges are crossed
  //const double delta_lambda = ((1 / nu_lower_first_initial) - (1 / nu_upper_last_initial)) / RADFIELDBINCOUNT; // upper limit if no edges are crossed

  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
  {
    //const double delta_nu = pow(prev_nu_upper,2) * delta_lambda; // equally spaced in wavelength
    radfieldbin_nu_upper[binindex] = nu_lower_first_initial + (binindex + 1) * delta_nu;

    // Align the bin edges with bound-free edges, except for the last one
    if (binindex < RADFIELDBINCOUNT - 1)
    {
      for (int i = 0; i < nbfcontinua_ground; i++)
      {
        const double nu_edge = phixslist[tid].groundcont[i].nu_edge;
        const int Z = get_element(phixslist[tid].groundcont[i].element);
        const int ion_stage = get_ionstage(phixslist[tid].groundcont[i].element, phixslist[tid].groundcont[i].ion);
        const int level = phixslist[tid].groundcont[i].level;
        //printout("bf edge at %g, nu_lower_first %g, nu_upper_last %g\n",nu_edge,nu_lower_first,nu_upper_last);
        if (binindex == 0 && ((nu_edge < nu_lower_first_initial) || (nu_edge > nu_upper_last_initial)))
        {
          printout("Missed bf edge at %12.5e Hz (%6.2f eV, %6.1f A), nu_lower_first %11.5e Hz, nu_upper_last %11.5e Hz, Z=%d ion_stage %d level %d\n",
                   nu_edge, H * nu_edge / EV, 1e8 * CLIGHT / nu_edge, nu_lower_first_initial, nu_upper_last_initial, Z, ion_stage, level);
        }

        if ((nu_edge > prev_nu_upper) && (nu_edge < radfield_get_bin_nu_upper(binindex)))
        {
          radfieldbin_nu_upper[binindex] = nu_edge;
          printout("Shifting bin %d nu_upper to bf edge at %12.5e Hz (%6.2f eV, %6.1f A) for Z=%d ion_stage %d level %d\n",
                   binindex, nu_edge, H * nu_edge / EV, 1e8 * CLIGHT / nu_edge, Z, ion_stage, level);
        }
      }
    }
    prev_nu_upper = radfield_get_bin_nu_upper(binindex);
  }


  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    radfield_set_J_normfactor(modelgridindex, -1.0);
    // printout("DEBUGCELLS: cell %d associated_cells %d\n", modelgridindex, mg_associated_cells[modelgridindex]);
    if (mg_associated_cells[modelgridindex] > 0)
    {
      radfieldbin_previous[modelgridindex] = (struct radfieldbin_previous *) calloc(RADFIELDBINCOUNT, sizeof(struct radfieldbin_previous));
      radfieldbin_current[modelgridindex] = (struct radfieldbin_current *) calloc(RADFIELDBINCOUNT, sizeof(struct radfieldbin_current));

      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      {
        radfieldbin_previous[modelgridindex][binindex].prev_J_normed = -1.;
        radfieldbin_previous[modelgridindex][binindex].prev_nuJ_normed = -1.;
        radfieldbin_previous[modelgridindex][binindex].prev_contribcount = 0;
        radfieldbin_current[modelgridindex][binindex].J_raw = 0.;
        radfieldbin_current[modelgridindex][binindex].nuJ_raw = 0.;
        radfieldbin_current[modelgridindex][binindex].contribcount = 0;
        radfieldbin_current[modelgridindex][binindex].W = -1.;
        radfieldbin_current[modelgridindex][binindex].T_R = -1.;
        // radfieldbin_current[modelgridindex][binindex].fit_type = FIT_DILUTED_BLACKBODY;
      }

    }
  }
  radfield_initialized = true;
}


static double radfield_get_bin_J(int modelgridindex, int binindex, bool averaged)
// get the normalised J_nu, from the current or current and previous timestep (averaged)
{
  if (J_normfactor[modelgridindex] <= 0.0)
  {
    printout("radfield: Fatal error: radfield_get_bin_J called before J_normfactor set for modelgridindex %d, = %g",modelgridindex,J_normfactor[modelgridindex]);
    abort();
  }
  else if (modelgridindex >= MMODELGRID)
  {
    printout("radfield: Fatal error: radfield_get_bin_J called before on modelgridindex %d >= MMODELGRID",modelgridindex);
    abort();
  }
  const double J_current = radfieldbin_current[modelgridindex][binindex].J_raw * J_normfactor[modelgridindex];
  if (!averaged || radfieldbin_previous[modelgridindex][binindex].prev_J_normed < 0.)
    return J_current;
  else
    return (J_current + radfieldbin_previous[modelgridindex][binindex].prev_J_normed) / 2;
}


static void radfield_set_bin_J(int modelgridindex, int binindex, double value)
// get the normalised J_nu, from the current or current and previous timestep (averaged)
{
  if (J_normfactor[modelgridindex] <= 0.0)
  {
    printout("radfield: Fatal error: radfield_set_bin_J called before J_normfactor set for modelgridindex %d",modelgridindex);
    abort();
  }
  else if (modelgridindex >= MMODELGRID)
  {
    printout("radfield: Fatal error: radfield_set_bin_J called before on modelgridindex %d >= MMODELGRID",modelgridindex);
    abort();
  }
  radfieldbin_current[modelgridindex][binindex].J_raw = value / J_normfactor[modelgridindex];
}


static double radfield_get_bin_nuJ(int modelgridindex, int binindex, bool averaged)
{
  if (J_normfactor[modelgridindex] <= 0.0)
  {
    printout("radfield: Fatal error: radfield_get_bin_nuJ called before J_normfactor set for modelgridindex %d",modelgridindex);
    abort();
  }
  else if (modelgridindex >= MMODELGRID)
  {
    printout("radfield: Fatal error: radfield_get_bin_nuJ called before on modelgridindex %d >= MMODELGRID",modelgridindex);
    abort();
  }
  const double nuJ_current = radfieldbin_current[modelgridindex][binindex].nuJ_raw * J_normfactor[modelgridindex];
  if (!averaged || radfieldbin_previous[modelgridindex][binindex].prev_nuJ_normed < 0.)
    return nuJ_current;
  else
    return (nuJ_current + radfieldbin_previous[modelgridindex][binindex].prev_nuJ_normed) / 2;
}


static inline
double radfield_get_bin_nu_bar(int modelgridindex, int binindex)
// importantly, this is average beween the current and previous timestep
{
  const double nuJ_sum = radfield_get_bin_nuJ(modelgridindex, binindex, true);
  const double J_sum = radfield_get_bin_J(modelgridindex, binindex, true);
  return nuJ_sum / J_sum;
}


static inline
double radfield_get_bin_nu_lower(int binindex)
{
  if (binindex > 0)
    return radfieldbin_nu_upper[binindex - 1];
  else
    return nu_lower_first_initial;
}


static inline
int radfield_get_bin_contribcount(int modelgridindex, int binindex, bool averaged)
// averaged with the previous timestep
{
  const int contribcount = radfieldbin_current[modelgridindex][binindex].contribcount;
  if (!averaged)
    return contribcount;
  else
    return contribcount + radfieldbin_previous[modelgridindex][binindex].prev_contribcount;
}


static inline
float radfield_get_bin_W(int modelgridindex, int binindex)
{
  return radfieldbin_current[modelgridindex][binindex].W;
}


static inline
float radfield_get_bin_T_R(int modelgridindex, int binindex)
{
  return radfieldbin_current[modelgridindex][binindex].T_R;
}


static inline
int radfield_select_bin(double nu)
{
  // linear search one by one until found
  if (nu >= radfieldbin_nu_upper[RADFIELDBINCOUNT - 1])
    return -1; // out of range, nu higher than highest bin
  else if (nu < radfield_get_bin_nu_lower(0))
    return -2; // out of range, nu lower than lowest bin
  else
  {
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      if (radfieldbin_nu_upper[binindex] > nu)
      {
        return binindex;
      }
    }

    // binary search for bin with nu_lower <= nu > nu_upper
    // int low = 0;
    // int high = RADFIELDBINCOUNT - 1;
    // while (low <= high)
    // {
    //   int mid = low + ((high - low) / 2);
    //   if (radfieldbin_nu_upper[mid] <= nu)
    //   {
    //     low = mid + 1;
    //   }
    //   else if (radfield_get_bin_nu_lower(mid) > nu)
    //   {
    //     high = mid - 1;
    //   }
    //   else
    //   {
    //     return mid;
    //   }
    //  }
    assert(false);
  }
}


void radfield_write_to_file(int modelgridindex, int timestep)
{
  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
# ifdef _OPENMP
# pragma omp critical (radfield_out_file)
  {
# endif
    if (!radfield_initialized)
    {
      printout("Call to radfield_write_to_file before radfield_init\n");
      abort();
    }

    int totalcontribs = 0;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      totalcontribs += radfield_get_bin_contribcount(modelgridindex, binindex, false);

    for (int binindex = -1; binindex < RADFIELDBINCOUNT; binindex++)
    {
      double nu_lower = 0.0;
      double nu_upper = 0.0;
      double nuJ_out = 0.0;
      double J_out = 0.0;
      float T_R = 0.0;
      float W = 0.0;
      double J_nu_bar = 0.0;
      int contribcount = 0;

      if (binindex >= 0)
      {
        nu_lower = radfield_get_bin_nu_lower(binindex);
        nu_upper = radfield_get_bin_nu_upper(binindex);
        nuJ_out = radfield_get_bin_nuJ(modelgridindex, binindex, true);
        J_out = radfield_get_bin_J(modelgridindex, binindex, true);
        T_R = radfield_get_bin_T_R(modelgridindex, binindex);
        W = radfield_get_bin_W(modelgridindex, binindex);
        J_nu_bar = J_out / (nu_upper - nu_lower);
        contribcount = radfield_get_bin_contribcount(modelgridindex, binindex, true);
      }
      else
      {
        nuJ_out = nuJ[modelgridindex];
        J_out = J[modelgridindex];
        T_R = get_TR(modelgridindex);
        W = get_W(modelgridindex);
        contribcount = totalcontribs;
      }

      fprintf(radfieldfile,"%8d %15d %8d %11.5e %11.5e %9.3e %9.3e %9.3e %9d %9.1f %12.5e\n",
              timestep,modelgridindex,binindex,nu_lower,nu_upper,nuJ_out,J_out,J_nu_bar,contribcount,T_R,W);
    }
    fflush(radfieldfile);
# ifdef _OPENMP
  }
# endif
  }
}


void radfield_close_file(void)
{
  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    fclose(radfieldfile);

    for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
    {
      if (mg_associated_cells[modelgridindex] > 0)
      {
        free(radfieldbin_current[modelgridindex]);
        free(radfieldbin_previous[modelgridindex]);
      }
    }
  }
}


void radfield_zero_estimators(int modelgridindex)
// set up the new bins and clear the estimators in preparation
// for a timestep
{
  J[modelgridindex] = 0.; // this is required even if FORCE_LTE is on
#ifndef FORCE_LTE
  nuJ[modelgridindex] = 0.;

  if (MULTIBIN_RADFIELD_MODEL_ON && radfield_initialized && mg_associated_cells[modelgridindex] > 0)
  {
    // printout("radfield: zeroing estimators in %d bins in cell %d\n",RADFIELDBINCOUNT,modelgridindex);

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      if (J_normfactor[modelgridindex] >= 0.)
      {
        radfieldbin_previous[modelgridindex][binindex].prev_J_normed = radfield_get_bin_J(modelgridindex, binindex, false);
        radfieldbin_previous[modelgridindex][binindex].prev_nuJ_normed = radfield_get_bin_nuJ(modelgridindex, binindex, false);
        radfieldbin_previous[modelgridindex][binindex].prev_contribcount = radfield_get_bin_contribcount(modelgridindex, binindex, false);
      }
      radfieldbin_current[modelgridindex][binindex].J_raw = 0.0;
      radfieldbin_current[modelgridindex][binindex].nuJ_raw = 0.0;
      radfieldbin_current[modelgridindex][binindex].contribcount = 0;
    }
    radfield_set_J_normfactor(modelgridindex, -1.0);
  }
#endif
}


inline
void radfield_update_estimators(int modelgridindex, double distance_e_cmf, double nu_cmf)
{
  #ifdef _OPENMP
    #pragma omp atomic
  #endif
  J[modelgridindex] += distance_e_cmf;
  #ifdef DEBUG_ON
    if (!isfinite(J[modelgridindex]))
    {
      printout("[fatal] update_estimators: estimator becomes non finite: distance_e_cmf %g, nu_cmf %g ... abort\n",distance_e_cmf,nu_cmf);
      abort();
    }
  #endif

#ifndef FORCE_LTE
  #ifdef _OPENMP
    #pragma omp atomic
  #endif
  nuJ[modelgridindex] += distance_e_cmf * nu_cmf;
  #ifdef DEBUG_ON
    if (!isfinite(nuJ[modelgridindex]))
    {
      printout("[fatal] update_estimators: estimator becomes non finite: distance_e_cmf %g, nu_cmf %g ... abort\n",distance_e_cmf,nu_cmf);
      abort();
    }
  #endif

  // int binindex = 0;
  // if (nu_cmf <= radfield_get_bin_nu_lower(modelgridindex,binindex))
  // {
  //   #ifdef DEBUG_ON
  //   printout("radfield: Extending nu_lower_first from %g down to %g\n",nu_lower_first,nu_cmf);
  //   #endif
  //   nu_lower_first = nu_cmf;
  // }
  // else if (nu_cmf > radfieldbin_nu_upper[modelgridindex][RADFIELDBINCOUNT - 1])
  // {
  //   binindex = RADFIELDBINCOUNT - 1;
  //   #ifdef DEBUG_ON
  //   printout("radfield: Extending nu_upper_last from %g up to %g\n",radfield_get_bin_nu_upper(modelgridindex,binindex),nu_cmf);
  //   #endif
  //   radfield_get_bin_nu_upper(modelgridindex, binindex) = nu_cmf;
  // }
  // else
  // {
  //   binindex = radfield_select_bin(modelgridindex,nu_cmf);
  // }
  const int binindex = radfield_select_bin(nu_cmf);

  if (binindex >= 0)
  {
    {
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      radfieldbin_current[modelgridindex][binindex].J_raw += distance_e_cmf;
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      radfieldbin_current[modelgridindex][binindex].nuJ_raw += distance_e_cmf * nu_cmf;
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      radfieldbin_current[modelgridindex][binindex].contribcount += 1;
    }
  }
  // else
  // {
  //   printout("WARNING: radfield_update_estimators dropping packet contribution for nu_cmf %g\n",
  //            nu_cmf);
  //   printout("           modelgridindex %d binindex %d nu_lower_first %g nu_upper_last %g \n",
  //            modelgridindex, binindex, nu_lower_first, radfield_get_bin_nu_upper(modelgridindex,RADFIELDBINCOUNT - 1));
  // }
#endif
}


double radfield(double nu, int modelgridindex)
// mean intensity J_nu
{
  if (MULTIBIN_RADFIELD_MODEL_ON && (nts_global >= FIRST_NLTE_RADFIELD_TIMESTEP))
  {
    const int binindex = radfield_select_bin(nu);
    if (binindex >= 0)
    {
      const struct radfieldbin_current *restrict const bin = &radfieldbin_current[modelgridindex][binindex];
      if (bin->W >= 0.)
      {
        // if (bin->fit_type == FIT_DILUTED_BLACKBODY)
        {
          const double J_nu = radfield_dbb(nu, bin->T_R, bin->W);
          /*if (fabs(J_nu / J_nu_fullspec - 1.0) > 0.5)
          {
            printout("WARNING: radfield: significant discrepancy. J_nu_fullspec %g, J_nu %g, nu %g bin->W %g bin->T_R %g\n",
                     J_nu_fullspec, J_nu, nu, bin->W, bin->T_R);
          }*/
          return J_nu;
        }
        // else
        // {
        //   return bin->W;
        // }
      }
      else //W < 0
      {
        //printout("WARNING: Radfield modelgridindex %d binindex %d has W_bin=%g<0, using W %g T_R %g nu %g\n",
        //         modelgridindex, binindex, W_bin, W_fullspec, T_R_fullspec, nu);
      }
    }
    else //binindex < 0
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

  const float T_R_fullspec = get_TR(modelgridindex);
  const float W_fullspec   = get_W(modelgridindex);
  const double J_nu_fullspec = radfield_dbb(nu, T_R_fullspec, W_fullspec);
  return J_nu_fullspec;
}


static double gsl_integrand_planck(double nu, void *restrict paras)
{
  double T_R = ((gsl_planck_integral_paras *) paras)->T_R;
  enum_prefactor prefactor = ((gsl_planck_integral_paras *) paras)->prefactor;

  double integrand = TWOHOVERCLIGHTSQUARED * pow(nu,3) / (expm1(HOVERKB * nu / T_R));

  if (prefactor == TIMES_NU)
    integrand *= nu;

  return integrand;
}


static double planck_integral(double T_R, double nu_lower, double nu_upper, enum_prefactor prefactor)
{
  double integral = 0.;

  double error = 0.;
  const double epsrel = 1e-10;
  const double epsabs = 0.;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(65536);

  gsl_planck_integral_paras intparas;
  intparas.T_R = T_R;
  intparas.prefactor = prefactor;

  gsl_function F_planck;
  F_planck.function = &gsl_integrand_planck;
  F_planck.params = &intparas;

  gsl_error_handler_t *gsl_error_handler = gsl_set_error_handler_off();
  int status = gsl_integration_qag(&F_planck, nu_lower, nu_upper, epsabs, epsrel, 65536, GSL_INTEG_GAUSS61, w, &integral, &error);
  if (status != 0)
  {
    printout("planck_integral integrator status %d, GSL_FAILURE= %d. Integral value %g, setting to zero.\n", status,GSL_FAILURE,integral);
    integral = 0.;
  }
  gsl_set_error_handler(gsl_error_handler);

  gsl_integration_workspace_free(w);

  return integral;
}


static double planck_integral_analytic(double T_R, double nu_lower, double nu_upper, enum_prefactor prefactor)
{
  double integral = 0.;

  if (prefactor == TIMES_NU)
  {
    double debye_upper = gsl_sf_debye_4(HOVERKB * nu_upper / T_R) * pow(nu_upper,4);
    double debye_lower = gsl_sf_debye_4(HOVERKB * nu_lower / T_R) * pow(nu_lower,4);
    integral = TWOHOVERCLIGHTSQUARED * (debye_upper - debye_lower) * T_R / HOVERKB / 4;
  }
  else
  {
    double debye_upper = gsl_sf_debye_3(HOVERKB * nu_upper / T_R) * pow(nu_upper,3);
    double debye_lower = gsl_sf_debye_3(HOVERKB * nu_lower / T_R) * pow(nu_lower,3);
    integral = TWOHOVERCLIGHTSQUARED * (debye_upper - debye_lower) * T_R / HOVERKB / 3;

    if (integral == 0.)
    {
      /*double upperexp = exp(HOVERKB * nu_upper / T_R);
      double upperint = - pow(nu_upper,4) / 4
                        + pow(nu_upper,3) * log(1 - upperexp) / HOVERKB
                        + 3 * pow(nu_upper,2) * polylog(2,upperexp) / pow(HOVERKB,2)
                        - 6 * nu_upper * polylog(3,upperexp) / pow(HOVERKB,3)
                        + 6 * polylog(4,upperexp) / pow(HOVERKB,4);
      double lowerexp = exp(HOVERKB * nu_lower / T_R);
      double lowerint = - pow(nu_lower,4) / 4
                        + pow(nu_lower,3) * log(1 - lowerexp) / HOVERKB
                        + 3 * pow(nu_lower,2) * polylog(2,lowerexp) / pow(HOVERKB,2)
                        - 6 * nu_lower * polylog(3,lowerexp) / pow(HOVERKB,3)
                        + 6 * polylog(4,lowerexp) / pow(HOVERKB,4);
      double integral2 = TWOHOVERCLIGHTSQUARED * (upperint - lowerint);

      printout("planck_integral_analytic is zero. debye_upper %g debye_lower %g. Test alternative %g\n",
               debye_upper,debye_lower,integral2);*/
    }
  }

  return integral;
}


static double delta_nu_bar(double T_R, void *restrict paras)
// difference between the average nu and the average nu of a planck function
// at temperature T_R, in the frequency range corresponding to a bin
{
  int modelgridindex = ((gsl_T_R_solver_paras *) paras)->modelgridindex;
  int binindex = ((gsl_T_R_solver_paras *) paras)->binindex;

  double nu_lower = radfield_get_bin_nu_lower(binindex);
  double nu_upper = radfield_get_bin_nu_upper(binindex);

  double nu_bar_estimator = radfield_get_bin_nu_bar(modelgridindex, binindex);

  double nu_times_planck_numerical = planck_integral(T_R, nu_lower, nu_upper, TIMES_NU);
  double planck_integral_numerical = planck_integral(T_R, nu_lower, nu_upper, ONE);
  double nu_bar_planck_T_R = nu_times_planck_numerical / planck_integral_numerical;

  /*double nu_times_planck_integral = planck_integral_analytic(T_R, nu_lower, nu_upper, TIMES_NU);
  double planck_integral_result = planck_integral_analytic(T_R, nu_lower, nu_upper, ONE);
  double nu_bar_planck = nu_times_planck_integral / planck_integral_result;

  //printout("nu_bar %g nu_bar_planck(T=%g) %g\n",nu_bar,T_R,nu_bar_planck);

  if (!isfinite(nu_bar_planck))
  {
    double nu_times_planck_numerical = planck_integral(T_R, nu_lower, nu_upper, TIMES_NU);
    double planck_integral_numerical = planck_integral(T_R, nu_lower, nu_upper, ONE);
    double nu_bar_planck_numerical = nu_times_planck_numerical / planck_integral_numerical;

    printout("planck_integral_analytic is %g. Replacing with numerical result of %g.\n",nu_bar_planck,nu_bar_planck_numerical);
    nu_bar_planck = nu_bar_planck_numerical;
  }*/

  double delta_nu_bar = nu_bar_planck_T_R - nu_bar_estimator;

  if (!isfinite(delta_nu_bar))
  {
    printout("delta_nu_bar is %g. nu_bar_planck_T_R %g nu_times_planck_numerical %g planck_integral_numerical %g nu_bar_estimator %g\n",
             nu_bar_planck_T_R, nu_times_planck_numerical, planck_integral_numerical, nu_bar_estimator);
  }

  //double delta_nu_bar = nu_bar_planck_T_R / nu_bar_estimator - 1.0;

  //printout("delta_nu_bar %g nu_bar_planck %g\n",delta_nu_bar,nu_bar_planck);

  return delta_nu_bar;
}


static float find_T_R(int modelgridindex, int binindex)
{
  double T_R = 0.0;

  gsl_T_R_solver_paras paras;
  paras.modelgridindex = modelgridindex;
  paras.binindex = binindex;

  /// Check whether the equation has a root in [T_min,T_max]
  double delta_nu_bar_min = delta_nu_bar(T_R_min,&paras);
  double delta_nu_bar_max = delta_nu_bar(T_R_max,&paras);

  // printout("find_T_R: bin %4d delta_nu_bar(T_R_min) %g, delta_nu_bar(T_R_max) %g\n",
  //          binindex, delta_nu_bar_min,delta_nu_bar_max);

  if (!isfinite(delta_nu_bar_min) || !isfinite(delta_nu_bar_max))
    delta_nu_bar_max = delta_nu_bar_min = -1;

  if (delta_nu_bar_min * delta_nu_bar_max < 0)
  {
    /// If there is a root in the interval, solve for T_R

    const double epsrel = 1e-4;
    const double epsabs = 0.;
    const int maxit = 100;

    gsl_function find_T_R_f;
    find_T_R_f.function = &delta_nu_bar;
    find_T_R_f.params = &paras;

    ///one dimensional gsl root solver, bracketing type
    gsl_root_fsolver *T_R_solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(T_R_solver, &find_T_R_f, T_R_min, T_R_max);
    int iteration_num = 0;
    int status;
    do
    {
      iteration_num++;
      gsl_root_fsolver_iterate(T_R_solver);
      T_R = gsl_root_fsolver_root(T_R_solver);

      const double T_R_lower = gsl_root_fsolver_x_lower(T_R_solver);
      const double T_R_upper = gsl_root_fsolver_x_upper(T_R_solver);
      status = gsl_root_test_interval(T_R_lower,T_R_upper,epsabs,epsrel);

      //printout("find_T_R: bin %4d iter %d, T_R is between %7.1f and %7.1f, guess %7.1f, delta_nu_bar %g, status %d\n",
      //         binindex,iteration_num,T_R_lower,T_R_upper,T_R,delta_nu_bar(T_R,&paras),status);
    }
    while (status == GSL_CONTINUE && iteration_num < maxit);

    if (status == GSL_CONTINUE)
      printout("[warning] find_T_R: T_R did not converge within %d iterations\n", maxit);

    gsl_root_fsolver_free(T_R_solver);
  }
  else if (delta_nu_bar_max < 0)
  {
    /// Thermal balance equation always negative ===> T_R = T_min
    /// Calculate the rates again at this T_e to print them to file
    T_R = T_R_max;
    printout("find_T_R: cell %d bin %4d no solution in interval, clamping to T_R_max=%g\n",
             modelgridindex, binindex, T_R_max);
  }
  else
  {
    T_R = T_R_min;
    printout("find_T_R: cell %d bin %4d no solution in interval, clamping to T_R_min=%g\n",
             modelgridindex, binindex, T_R_min);
  }

  return T_R;
}


static void get_radfield_params_fullspec(double J, double nuJ, int modelgridindex, float *T_J, float *T_R, float *W)
{
  const double nubar = nuJ / J;
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
    *T_J = pow(PI / STEBO * J, 1 / 4.);
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

    *T_R = H * nubar / KB / 3.832229494;
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

    *W = PI * J / STEBO / pow(*T_R, 4);
  }
}


void radfield_fit_parameters(int modelgridindex, int timestep)
// finds the best fitting W and temperature parameters in each spectral bin
// using J and nuJ
{
  float T_J;
  float T_R;
  float W;
  get_radfield_params_fullspec(J[modelgridindex], nuJ[modelgridindex], modelgridindex, &T_J, &T_R, &W);
  set_TJ(modelgridindex, T_J);
  set_TR(modelgridindex, T_R);
  set_W(modelgridindex, W);

  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    if (J_normfactor[modelgridindex] <= 0)
    {
      printout("radfield: FATAL J_normfactor = %g in cell %d at call to radfield_fit_parameters", J_normfactor[modelgridindex], modelgridindex);
      abort();
    }

    const float T_R_fullspec = get_TR(modelgridindex);
    const double J_fullspec = J[modelgridindex];
    //double planck_integral_zero_inf = STEBO * pow(T_R_fullspec,4) / PI;

    printout("Full-spectrum fit radfield params for cell %d at timestep %d: J %g, T_R %g, W %g\n",
             modelgridindex, timestep, J_fullspec, T_R_fullspec, get_W(modelgridindex));

    double J_bin_sum = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      J_bin_sum += radfield_get_bin_J(modelgridindex, binindex, true);

    printout("radfield bins sum to J of %g (%.1f%% of total J).\n",
             J_bin_sum, 100. * J_bin_sum / J[modelgridindex]);
    printout("radfield: Finding parameters for %d bins...\n", RADFIELDBINCOUNT);

    double J_bin_max = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      const double J_bin = radfield_get_bin_J(modelgridindex, binindex, true);
      if (J_bin > J_bin_max)
        J_bin_max = J_bin;
    }

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      const double nu_lower = radfield_get_bin_nu_lower(binindex);
      const double nu_upper = radfield_get_bin_nu_upper(binindex);
      if (boost_region_on && boost_region_nu_lower < nu_upper && boost_region_nu_upper > nu_lower)
      {
        printout("Artificially boosting bin %d\n",binindex);
        radfield_set_bin_J(modelgridindex, binindex, J_bin_max * boost_region_factor);
      }
      const double J_bin = radfield_get_bin_J(modelgridindex, binindex, true);
      float T_R_bin = -1.0;
      float W_bin = -1.0;
      int contribcount = radfield_get_bin_contribcount(modelgridindex, binindex, true);

      if (contribcount > 10)
      {
        // // enum_bin_fit_type bin_fit_type = radfieldbin_current[modelgridindex][binindex].fit_type;
        // if (bin_fit_type == FIT_DILUTED_BLACKBODY)
        {
          T_R_bin = find_T_R(modelgridindex, binindex);

          double planck_integral_result = planck_integral(T_R_bin, nu_lower, nu_upper, ONE);

          W_bin = J_bin / planck_integral_result;

          if (W_bin > 1e4)
          {
            printout("W %g too high, trying setting T_R of bin %d to %g. J_bin %g planck_integral %g\n",
                     W_bin, binindex, T_R_max, planck_integral_result);
            planck_integral_result = planck_integral(T_R_max, nu_lower, nu_upper, ONE);
            W_bin = J_bin / planck_integral_result;
            if (W_bin > 1e4)
            {
              printout("W still very high, W=%g. Zeroing bin...\n", W_bin);
              T_R_bin = -99.0;
              W_bin = 0.;
            }
            else
            {
              printout("new W is %g. Continuing with this value\n", W_bin);
            }
          }

        }
        // else if (bin_fit_type == FIT_CONSTANT)
        // {
        //   T_R_bin = -1.;
        //   W_bin = J_bin / (nu_upper - nu_lower);
        // }
        // else
        // {
        //   printout("radfield_fit_parameters: unknown fit type %d for bin %d\n", bin_fit_type, binindex);
        //   T_R_bin = -1.;
        //   W_bin = -1.;
        // }
      }
      else if (contribcount == 0)
      {
        T_R_bin = 0.;
        W_bin = 0.;
      }
      else
      {
        T_R_bin = -1;
        W_bin = -1;
      }
      radfieldbin_current[modelgridindex][binindex].T_R = T_R_bin;
      radfieldbin_current[modelgridindex][binindex].W = W_bin;
    }

    /*for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      double J_bin = radfield_get_bin_J(modelgridindex,binindex);
      double T_R_bin = radfield_get_bin_T_R(modelgridindex,binindex);
      double W_bin = radfield_get_bin_W(modelgridindex,binindex);
      printout("bin %4d: J %g, T_R %7.1f, W %12.5e\n",
             binindex, J_bin, T_R_bin, W_bin);
    }*/
    if (timestep % 2 == 0)
      radfield_write_to_file(modelgridindex, timestep);
  }
}


void radfield_set_J_normfactor(int modelgridindex, double normfactor)
{
  J_normfactor[modelgridindex] = normfactor;
}


void radfield_normalise_J(const int modelgridindex, const double estimator_normfactor_over4pi)
{
  assert(isfinite(J[modelgridindex]));
  J[modelgridindex] *= estimator_normfactor_over4pi;
}


void radfield_normalise_nuJ(const int modelgridindex, const double estimator_normfactor_over4pi)
{
  assert(isfinite(nuJ[modelgridindex]));
  nuJ[modelgridindex] *= estimator_normfactor_over4pi;
}


double get_T_R_from_J(const int modelgridindex)
{
  double T_R = pow(PI / STEBO * J[modelgridindex], 1/4);
  if (isfinite(T_R))
  {
    /// Make sure that T is in the allowed temperature range.
    if (T_R > MAXTEMP)
      T_R = MAXTEMP;
    if (T_R < MINTEMP)
      T_R = MINTEMP;
  }
  else
  {
    /// keep old value of T_R
    printout("[warning] update_grid: T_R estimator infinite in cell %d, use value of last timestep\n", modelgridindex);
    T_R = modelgrid[modelgridindex].TR;
  }
  return T_R;
}


#ifdef DO_TITER
void radfield_titer_J(const int modelgridindex)
{
  if (J_reduced_save[modelgridindex] >= 0)
  {
    J[modelgridindex] = (J[modelgridindex] + J_reduced_save[modelgridindex]) / 2;
  }
  J_reduced_save[modelgridindex] = J[modelgridindex];
}


#ifndef FORCE_LTE
void radfield_titer_nuJ(const int modelgridindex)
{
  if (nuJ_reduced_save[modelgridindex] >= 0)
  {
    nuJ[modelgridindex] = (nuJ[modelgridindex] + nuJ_reduced_save[modelgridindex]) / 2;
  }
  nuJ_reduced_save[modelgridindex] = nuJ[modelgridindex];
}
#endif
#endif


#ifdef MPI_ON
void radfield_reduce_estimators(void)
// reduce and broadcast (allreduce) the estimators for J and nuJ in all bins
{
  MPI_Allreduce(MPI_IN_PLACE, &J, MMODELGRID, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #ifndef FORCE_LTE
    MPI_Allreduce(MPI_IN_PLACE, &nuJ, MMODELGRID, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif

  if (!MULTIBIN_RADFIELD_MODEL_ON)
    return;
  printout("Reducing binned radiation field estimators\n");

  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    // printout("DEBUGCELLS: cell %d associated_cells %d\n", modelgridindex, mg_associated_cells[modelgridindex]);
    if (mg_associated_cells[modelgridindex] > 0)
    {
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      {
        // printout("MPI: pre-MPI_Allreduce, process %d modelgrid %d binindex %d has a individual contribcount of %d\n",my_rank,modelgridindex,binindex,radfieldbins[modelgridindex][binindex].contribcount);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &radfieldbin_current[modelgridindex][binindex].J_raw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &radfieldbin_current[modelgridindex][binindex].nuJ_raw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &radfieldbin_current[modelgridindex][binindex].contribcount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // printout("MPI: After MPI_Allreduce: Process %d modelgrid %d binindex %d has a contribcount of %d\n",my_rank,modelgridindex,binindex,radfieldbins[modelgridindex][binindex].contribcount);
      }
    }
  }
}


void radfield_MPI_Bcast(int root, int my_rank, int nstart, int ndo)
// broadcast computed radfield results including parameters
// from the cells belonging to root process to all processes
{
  int sender_nstart;
  int sender_ndo;
  // double nu_lower_first;
  if (root == my_rank)
  {
    sender_nstart = nstart;
    sender_ndo = ndo;
    // if (ndo > 0)
      // printout("radfield_MPI_Bcast root process %d will send data for cells %d to %d\n", my_rank, sender_nstart, sender_nstart + sender_ndo - 1);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&sender_nstart, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&sender_ndo, 1, MPI_INT, root, MPI_COMM_WORLD);
  // if (my_rank != root && sender_ndo > 0)
    // printout("radfield_MPI_Bcast process %d will recieve data for cells %d to %d\n", my_rank, sender_nstart, sender_nstart + sender_ndo - 1);

  for (int modelgridindex = sender_nstart; modelgridindex < sender_nstart + sender_ndo; modelgridindex++)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&J_normfactor[modelgridindex], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    if (mg_associated_cells[modelgridindex] > 0)
    {
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      {
        MPI_Barrier(MPI_COMM_WORLD);
        // printout("radfield_MPI_Bcast bin %d T_R before: %g\n", binindex, radfieldbins[modelgridindex][binindex].T_R);
        MPI_Bcast(&radfieldbin_current[modelgridindex][binindex].W, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
        MPI_Bcast(&radfieldbin_current[modelgridindex][binindex].T_R, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
        MPI_Bcast(&radfieldbin_current[modelgridindex][binindex].J_raw, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Bcast(&radfieldbin_current[modelgridindex][binindex].nuJ_raw, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Bcast(&radfieldbin_current[modelgridindex][binindex].contribcount, 1, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(&radfieldbin_previous[modelgridindex][binindex].prev_J_normed, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Bcast(&radfieldbin_previous[modelgridindex][binindex].prev_nuJ_normed, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Bcast(&radfieldbin_previous[modelgridindex][binindex].prev_contribcount, 1, MPI_INT, root, MPI_COMM_WORLD);
        // printout("radfield_MPI_Bcast MPI_Bcast radfield bin %d for cell %d from process %d to %d\n", binindex, modelgridindex, root, my_rank);
        // printout("radfield_MPI_Bcast bin %d T_R after: %g\n", binindex, radfieldbins[modelgridindex][binindex].T_R);
      }
    }
  }
}
#endif


void radfield_write_restart_data(FILE *gridsave_file)
{
  if (!MULTIBIN_RADFIELD_MODEL_ON)
    return;

  printout("Writing restart data for binned radiation field\n");

  fprintf(gridsave_file, "%d\n", 30490824); // special number marking the beginning of radfield data
  fprintf(gridsave_file,"%d %lg %lg %lg %lg\n",
          RADFIELDBINCOUNT, nu_lower_first_initial, nu_upper_last_initial,
          T_R_min, T_R_max);

  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    fprintf(gridsave_file,"%d %lg\n", binindex, radfieldbin_nu_upper[binindex]);

  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    if (mg_associated_cells[modelgridindex] > 0)
    {
      fprintf(gridsave_file,"%d %lg\n", modelgridindex, J_normfactor[modelgridindex]);
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      {
        fprintf(gridsave_file,"%lg %lg %lg %lg %d %g %g %d\n",
                radfieldbin_current[modelgridindex][binindex].J_raw,
                radfieldbin_current[modelgridindex][binindex].nuJ_raw,
                radfieldbin_previous[modelgridindex][binindex].prev_J_normed,
                radfieldbin_previous[modelgridindex][binindex].prev_nuJ_normed,
                radfieldbin_previous[modelgridindex][binindex].prev_contribcount,
                radfieldbin_current[modelgridindex][binindex].W,
                radfieldbin_current[modelgridindex][binindex].T_R,
                radfieldbin_current[modelgridindex][binindex].contribcount);
                //radfieldbin_current[modelgridindex][binindex].fit_type
      }
    }
  }
}


void radfield_read_restart_data(FILE *gridsave_file)
{
  if (!MULTIBIN_RADFIELD_MODEL_ON)
    return;

  printout("Reading restart data for radiation field\n");

  if (!radfield_initialized)
  {
    printout("ERROR: Radiation field has not been initialised yet. Can't read saved state.\n");
    abort();
  }

  int code_check;
  fscanf(gridsave_file, "%d\n", &code_check);
  if (code_check != 30490824)
  {
    printout("ERROR: Beginning of radfield restart data not found!");
    abort();
  }

  int bincount_in;
  double T_R_min_in, T_R_max_in, nu_lower_first_initial_in, nu_upper_last_initial_in;
  fscanf(gridsave_file,"%d %lg %lg %lg %lg\n",
         &bincount_in, &nu_lower_first_initial_in, &nu_upper_last_initial_in,
         &T_R_min_in, &T_R_max_in);

  double nu_lower_first_ratio = nu_lower_first_initial_in / nu_lower_first_initial;
  if (nu_lower_first_ratio > 1.0) nu_lower_first_ratio = 1 / nu_lower_first_ratio;
  double nu_upper_last_ratio = nu_upper_last_initial_in / nu_upper_last_initial;
  if (nu_upper_last_ratio > 1.0) nu_upper_last_ratio = 1 / nu_upper_last_ratio;

  if (bincount_in != RADFIELDBINCOUNT || T_R_min_in != T_R_min || T_R_max_in != T_R_max ||
      nu_lower_first_ratio < 0.999 || nu_upper_last_ratio < 0.999)
  {
    printout("ERROR: gridsave file specifies %d bins, nu_lower_first_initial %lg nu_upper_last_initial %lg T_R_min %lg T_R_max %lg\n",
             bincount_in, nu_lower_first_initial_in, nu_upper_last_initial_in, T_R_min_in, T_R_max_in);
    printout("require %d bins, nu_lower_first_initial %lg nu_upper_last_initial %lg T_R_min %lg T_R_max %lg\n",
            RADFIELDBINCOUNT, nu_lower_first_initial, nu_upper_last_initial, T_R_min, T_R_max);
    abort();
  }

  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
  {
    int binindex_in;
    fscanf(gridsave_file,"%d %lg\n", &binindex_in, &radfieldbin_nu_upper[binindex]);
    assert(binindex_in == binindex);
  }

  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    if (mg_associated_cells[modelgridindex] > 0)
    {
      int mgi_in;
      fscanf(gridsave_file,"%d %lg\n", &mgi_in, &J_normfactor[modelgridindex]);
      if (mgi_in != modelgridindex)
      {
        printout("ERROR: expected data for cell %d but found cell %d\n", modelgridindex, mgi_in);
        abort();
      }
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      {
        fscanf(gridsave_file,"%lg %lg %lg %lg %d %g %g %d\n",
               &radfieldbin_current[modelgridindex][binindex].J_raw,
               &radfieldbin_current[modelgridindex][binindex].nuJ_raw,
               &radfieldbin_previous[modelgridindex][binindex].prev_J_normed,
               &radfieldbin_previous[modelgridindex][binindex].prev_nuJ_normed,
               &radfieldbin_previous[modelgridindex][binindex].prev_contribcount,
               &radfieldbin_current[modelgridindex][binindex].W,
               &radfieldbin_current[modelgridindex][binindex].T_R,
               &radfieldbin_current[modelgridindex][binindex].contribcount);
      }
    }
  }
}

inline
int radfield_integrate(
  const gsl_function *f, double nu_a, double nu_b, double epsabs,
  double epsrel, size_t limit, int key, gsl_integration_workspace *workspace,
  double *result, double *abserr)
{
  if (MULTIBIN_RADFIELD_MODEL_ON && (nts_global >= FIRST_NLTE_RADFIELD_TIMESTEP))
  {
    double *pts = malloc((RADFIELDBINCOUNT + 3) * sizeof(double));
    int binindex_a = radfield_select_bin(nu_a);
    int binindex_b = radfield_select_bin(nu_b);
    int npts = 0;
    pts[npts++] = nu_a;
    if (binindex_a == binindex_b) // both higher, both lower, or match the same bin
    {
      // region doesn't contain any bins
      pts[npts++] = nu_b;
    }
    else
    {
      if (binindex_a < 0) // a is below the first bin
      {
        binindex_a = 0;
        pts[npts++] = radfield_get_bin_nu_lower(0);
      }

      const int maxbinplusone = (binindex_b < 0) ? RADFIELDBINCOUNT : binindex_b;

      for (int binindex = binindex_a; binindex < maxbinplusone; binindex++)
        pts[npts++] = radfield_get_bin_nu_upper(binindex);

      pts[npts++] = nu_b;
    }
    // for (int e = 0; e < npts; e++)
    // {
    //   printout("radfield_integrate singular point number %d at nu %g, (nu_a %g, nu_b %g), radfield_low %g radfield_high %g\n",
    //            e, pts[e], nu_a, nu_b, radfield(pts[e] * 0.9999, 0), radfield(pts[e] * 1.0001, 0));
    // }
    const int status = gsl_integration_qagp(f, pts, npts, epsabs, epsrel, limit, workspace, result, abserr);
    free(pts);
    return status;
  }
  else
    return gsl_integration_qag(f, nu_a, nu_b, epsabs, epsrel, limit, key, workspace, result, abserr);
}