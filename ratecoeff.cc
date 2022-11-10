#include "ratecoeff.h"

#include <gsl/gsl_integration.h>

#include <algorithm>
#include <cmath>
#include <cstring>
// #define  _XOPEN_SOURCE
#define D_POSIX_SOURCE
#include <cstdio>

#include "artisoptions.h"
#include "atomic.h"
#include "grid.h"
#include "input.h"
#include "ltepop.h"
#include "macroatom.h"
#include "md5.h"
#include "radfield.h"
#include "sn3d.h"

// typedef struct gslintegration_ffheatingparas
// {
//   double T_e;
//   int cellnumber;
// } gslintegration_ffheatingparas;

// typedef struct gslintegration_bfheatingparas
// {
//   double nu_edge;
//   int cellnumber;
// } gslintegration_bfheatingparas;

static __managed__ double T_step;
__managed__ double T_step_log;

typedef struct {
  const double nu_edge;
  const double departure_ratio;
  const float *const photoion_xs;
  const float T_e;
  const int modelgridindex;
} gsl_integral_paras_gammacorr;

static char adatafile_hash[33];
static char compositionfile_hash[33];
static char phixsfile_hash[33];

static bool read_ratecoeff_dat(void)
/// Try to read in the precalculated rate coefficients from file
/// return true if successful or false otherwise
{
  FILE *ratecoeff_file = fopen("ratecoeff.dat", "r");
  if (ratecoeff_file != NULL) {
    /// Check whether current atomic data and temperature range match
    /// the precalculated rate coefficients
    bool fileisamatch = true;  // assume true until a mismatch is detected

    char adatafile_hash_in[33];
    assert_always(fscanf(ratecoeff_file, "%32s\n", adatafile_hash_in) == 1);
    printout("ratecoeff.dat: MD5 adata.txt = %s ", adatafile_hash_in);
    if (strcmp(adatafile_hash, adatafile_hash_in) == 0)
      printout("(pass)\n");
    else {
      printout("MISMATCH: MD5 adata.txt = %s\n", adatafile_hash);
      fileisamatch = false;
    }

    char compositionfile_hash_in[33];
    assert_always(fscanf(ratecoeff_file, "%32s\n", compositionfile_hash_in) == 1);
    printout("ratecoeff.dat: MD5 compositiondata.txt %s ", compositionfile_hash_in);
    if (strcmp(compositionfile_hash, compositionfile_hash_in) == 0)
      printout("(pass)\n");
    else {
      printout("\nMISMATCH: MD5 compositiondata.txt = %s\n", compositionfile_hash);
      fileisamatch = false;
    }

    char phixsfile_hash_in[33];
    assert_always(fscanf(ratecoeff_file, "%32s\n", phixsfile_hash_in) == 1);
    printout("ratecoeff.dat: MD5 %s = %s ", phixsdata_filenames[phixs_file_version], phixsfile_hash_in);
    if (strcmp(phixsfile_hash, phixsfile_hash_in) == 0)
      printout("(pass)\n");
    else {
      printout("\nMISMATCH: MD5 %s = %s\n", phixsdata_filenames[phixs_file_version], phixsfile_hash);
      fileisamatch = false;
    }

    if (fileisamatch) {
      float T_min, T_max;
      int in_tablesize;
      int in_nlines;
      assert_always(fscanf(ratecoeff_file, "%g %g %d %d\n", &T_min, &T_max, &in_tablesize, &in_nlines) == 4);
      printout("ratecoeff.dat: Tmin %g Tmax %g TABLESIZE %d nlines %d ", T_min, T_max, in_tablesize, in_nlines);

      if (T_min == MINTEMP && T_max == MAXTEMP && in_tablesize == TABLESIZE && in_nlines == globals::nlines) {
        printout("(pass)\n");
        // this is redundant if the adata and composition data matches, but have
        // to read through to maintain consistency with older files
        for (int element = 0; element < get_nelements(); element++) {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions; ion++) {
            int in_element, in_ionstage, in_levels, in_ionisinglevels;
            assert_always(fscanf(ratecoeff_file, "%d %d %d %d\n", &in_element, &in_ionstage, &in_levels,
                                 &in_ionisinglevels) == 4);
            const int nlevels = get_nlevels(element, ion);
            int ionisinglevels = get_ionisinglevels(element, ion);
            if (get_element(element) != in_element || get_ionstage(element, ion) != in_ionstage ||
                nlevels != in_levels || ionisinglevels != in_ionisinglevels) {
              printout(
                  "Levels or ionising levels count mismatch! element %d %d ionstage %d %d nlevels %d %d ionisinglevels "
                  "%d %d\n",
                  get_element(element), in_element, get_ionstage(element, ion), in_ionstage, nlevels, in_levels,
                  ionisinglevels, in_ionisinglevels);
              fileisamatch = false;
              break;
            }
          }
          if (!fileisamatch) {
            break;
          }
        }
      } else {
        printout("\nMISMATCH: this simulation has MINTEMP %g MAXTEMP %g TABLESIZE %d nlines %d\n", MINTEMP, MAXTEMP,
                 TABLESIZE, globals::nlines);
        fileisamatch = false;
      }
    }

    if (fileisamatch) {
      printout("Matching ratecoeff.dat file found. Readin this file ...\n");
      for (int element = 0; element < get_nelements(); element++) {
        int nions = get_nions(element) - 1;
        for (int ion = 0; ion < nions; ion++) {
          // nlevels = get_nlevels(element,ion);
          const int nlevels =
              get_ionisinglevels(element, ion);  /// number of ionising levels associated with current ion
          // int nbfcont = get_ionisinglevels(element,ion);     /// number of ionising levels of the current ion which
          // are used in the simulation
          for (int level = 0; level < nlevels; level++) {
            /// Loop over the phixs target states
            const int nphixstargets = get_nphixstargets(element, ion, level);
            for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
              /// Loop over the temperature grid
              for (int iter = 0; iter < TABLESIZE; iter++) {
                double alpha_sp, bfcooling_coeff, corrphotoioncoeff, bfheating_coeff;
                assert_always(fscanf(ratecoeff_file, "%lg %lg %lg %lg\n", &alpha_sp, &bfcooling_coeff,
                                     &corrphotoioncoeff, &bfheating_coeff) == 4);

                // assert_always(std::isfinite(alpha_sp) && alpha_sp >= 0);
                globals::spontrecombcoeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] = alpha_sp;

                // assert_always(std::isfinite(bfcooling_coeff) && bfcooling_coeff >= 0);
                globals::bfcooling_coeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] = bfcooling_coeff;

#if (!NO_LUT_PHOTOION)
                if (corrphotoioncoeff >= 0) {
                  globals::corrphotoioncoeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] =
                      corrphotoioncoeff;
                } else {
                  printout(
                      "ERROR: NO_LUT_PHOTOION is off, but there are no corrphotoioncoeff values in ratecoeff file\n");
                  abort();
                }
#endif
#if (!NO_LUT_BFHEATING)
                if (bfheating_coeff >= 0) {
                  globals::bfheating_coeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] =
                      bfheating_coeff;
                } else {
                  printout(
                      "ERROR: NO_LUT_BFHEATING is off, but there are no bfheating_coeff values in the ratecoeff "
                      "file\n");
                  abort();
                }
#endif
              }
            }
          }
        }
      }
      fclose(ratecoeff_file);
      return true;
    } else {
      printout("[info] ratecoefficients_init: ratecoeff.dat does not match current simulation. Recalculating...\n");
      fclose(ratecoeff_file);
      return false;
    }
  } else {
    printout("[info] ratecoefficients_init:  No ratecoeff.dat file available. Creating a new one...\n");
    return false;
  }
}

static void write_ratecoeff_dat(void) {
  FILE *ratecoeff_file = fopen_required("ratecoeff.dat", "w");
  fprintf(ratecoeff_file, "%32s\n", adatafile_hash);
  fprintf(ratecoeff_file, "%32s\n", compositionfile_hash);
  fprintf(ratecoeff_file, "%32s\n", phixsfile_hash);
  fprintf(ratecoeff_file, "%g %g %d %d\n", MINTEMP, MAXTEMP, TABLESIZE, globals::nlines);
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      fprintf(ratecoeff_file, "%d %d %d %d\n", get_element(element), get_ionstage(element, ion),
              get_nlevels(element, ion), get_ionisinglevels(element, ion));
    }
  }

  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element) - 1;
    for (int ion = 0; ion < nions; ion++) {
      // nlevels = get_nlevels(element,ion);
      const int nlevels = get_ionisinglevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        /// Loop over the phixs targets
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++) {
          /// Loop over the temperature grid
          for (int iter = 0; iter < TABLESIZE; iter++) {
            const double alpha_sp =
                globals::spontrecombcoeff[get_bflutindex(iter, element, ion, level, phixstargetindex)];
            const double bfcooling_coeff =
                globals::bfcooling_coeff[get_bflutindex(iter, element, ion, level, phixstargetindex)];
            fprintf(ratecoeff_file, "%g %g", alpha_sp, bfcooling_coeff);

#if (!NO_LUT_PHOTOION)
            const double corrphotoioncoeff =
                globals::corrphotoioncoeff[get_bflutindex(iter, element, ion, level, phixstargetindex)];
#else
            const double corrphotoioncoeff = -1.0;
#endif
#if (!NO_LUT_BFHEATING)
            const double bfheating_coeff =
                globals::bfheating_coeff[get_bflutindex(iter, element, ion, level, phixstargetindex)];
#else
            const double bfheating_coeff = -1.0;
#endif
            fprintf(ratecoeff_file, " %g %g\n", corrphotoioncoeff, bfheating_coeff);
          }
        }
      }
    }
  }
  fclose(ratecoeff_file);
}

///****************************************************************************
/// The following functions define the integrands for these rate coefficients
/// for use with libgsl integrators.
static double alpha_sp_integrand_gsl(const double nu, void *const voidparas)
/// Integrand to calculate the rate coefficient for spontaneous recombination
/// using gsl integrators.
{
  const gslintegration_paras *const params = static_cast<gslintegration_paras *>(voidparas);

  const float sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, params->nu_edge, nu);
  const double x = TWOOVERCLIGHTSQUARED * sigma_bf * pow(nu, 2) * exp(-HOVERKB * nu / params->T);
  /// in formula this looks like
  /// x = sigma_bf/H/nu * 2*H*pow(nu,3)/pow(CLIGHT,2) * exp(-H*nu/KB/T);

  /// set contributions from Lyman continuum artificially to zero to overcome it's large opacity
  return x;
}

static double alpha_sp_E_integrand_gsl(const double nu, void *const voidparas)
/// Integrand to calculate the rate coefficient for spontaneous recombination
/// using gsl integrators.
{
  const gslintegration_paras *const params = (gslintegration_paras *)voidparas;

  const float T = params->T;
  const double nu_edge = params->nu_edge;

  const float sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge, nu);
  const double x = TWOOVERCLIGHTSQUARED * sigma_bf * pow(nu, 3) / nu_edge * exp(-HOVERKB * nu / T);
  /// in formula this looks like
  /// x = sigma_bf/H/nu * 2*H*pow(nu,3)/pow(CLIGHT,2) * exp(-H*nu/KB/T);

  /// set contributions from Lyman continuum artificially to zero to overcome it's large opacity
  return x;
}

/*static double gamma_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the rate coefficient for photoionization
/// using gsl integrators.
{
  double T = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;

  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge,nu);

  /// Dependence on dilution factor W is linear. This allows to set it here to
  /// 1. and scale to its actual value later on.
  double x = sigma_bf / H / nu * radfield::dbb(nu,T,1.);
  //x = sigma_bf/H/nu * radfield::dbb(nu,T,1.);
  //if (HOVERKB*nu/T < 1e-2) x = sigma_bf * pow(nu,2)/(HOVERKB*nu/T);
  //else if (HOVERKB*nu/T >= 1e2) x = sigma_bf * pow(nu,2)*exp(-HOVERKB*nu/T);
  //else x = sigma_bf * pow(nu,2)/(exp(HOVERKB*nu/T)-1);

  return x;
}*/

#if (!NO_LUT_PHOTOION)
static double gammacorr_integrand_gsl(const double nu, void *const voidparas)
/// Integrand to calculate the rate coefficient for photoionization
/// using gsl integrators. Corrected for stimulated recombination.
{
  const gslintegration_paras *const params = (gslintegration_paras *)voidparas;

  const float T = params->T;
  const double nu_edge = params->nu_edge;

  const float sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge, nu);

  /// Dependence on dilution factor W is linear. This allows to set it here to
  /// 1. and scale to its actual value later on.
  /// Assumption T_e = T_R makes n_kappa/n_i * (n_i/n_kappa)* = 1
  return sigma_bf * ONEOVERH / nu * radfield::dbb(nu, T, 1) * (1 - exp(-HOVERKB * nu / T));
}
#endif

#if (!NO_LUT_BFHEATING)
static double approx_bfheating_integrand_gsl(const double nu, void *const voidparas)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximative way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  const gslintegration_paras *const params = (gslintegration_paras *)voidparas;

  const float T = params->T;
  const double nu_edge = params->nu_edge;

  const float sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge, nu);

  /// Precalculation for T_e=T_R and W=1
  const double x = sigma_bf * (1 - nu_edge / nu) * radfield::dbb(nu, T, 1) * (1 - exp(-HOVERKB * nu / T));

  /// Precalculation for a (T_R,T_e)-grid, but still W is assumed to be 1.
  /// The radfield part can be corrected later because of its linear dependence.
  /// But not the W in the stimulated correction term!
  /*double T_e  = ((gslintegration_paras *) paras)->T;
  double T_R  = ((gslintegration_paras *) paras)->T2;
  double E_threshold = nu_edge*H;
  double sf_Te = calculate_sahafact(element,ion,level,upperionlevel,T_e,E_threshold);
  double sf_TR = calculate_sahafact(element,ion,level,upperionlevel,T_R,E_threshold);
  x = sigma_bf*(1-nu_edge/nu)*radfield::dbb(nu,T_R,1) * (1 - sqrt(T_e/T_R) * sf_Te/sf_TR * exp(-H*nu/KB/T_e));*/

  return x;
}
#endif

/*double bfheating_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the modified rate coefficient for photoionization
/// using gsl integrators.
{
  double get_groundlevelpop(int cellnumber, int element, int ion);
  double x;

  int cellnumber = ((gslintegration_bfheatingparas *) paras)->cellnumber;
  double nu_edge = ((gslintegration_bfheatingparas *) paras)->nu_edge;

  float T_e = globals::cell[cellnumber].T_e;
  double T_R = globals::cell[cellnumber].T_R;
  double W = globals::cell[cellnumber].W;
  float nne = globals::cell[cellnumber].nne;

  double sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge,nu);
  double E_threshold = nu_edge*H;
  double sfac = calculate_sahafact(element,ion,level,upperionlevel,T_e,E_threshold);
  double nnionlevel = get_groundlevelpop(cellnumber,element,ion+1);

  x = sigma_bf*(1-nu_edge/nu)*radfield::dbb(nu,T_R,W) * (1-nnionlevel*nne/nnlevel*sf*exp(-H*nu/KB/T_e));
  return x;
}*/

static double bfcooling_integrand_gsl(const double nu, void *const voidparas)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximative way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  const gslintegration_paras *const params = (gslintegration_paras *)voidparas;

  const float T = params->T;
  const double nu_edge = params->nu_edge;

  const float sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge, nu);

  // return sigma_bf * (1-nu_edge/nu) * TWOHOVERCLIGHTSQUARED * pow(nu,3) * exp(-HOVERKB*nu/T);
  return sigma_bf * (nu - nu_edge) * TWOHOVERCLIGHTSQUARED * nu * nu * exp(-HOVERKB * nu / T);
}

/*static double bfcooling_integrand_gsl_2(double nu, void *paras)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximative way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  double T  = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;

  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge,nu);

  return sigma_bf*(1/nu_edge-1/nu) * TWOOVERCLIGHTSQUARED*pow(nu,3) * exp(-HOVERKB*nu/T);
}*/

/*static double stimulated_bfcooling_integrand_gsl(double nu, void *paras)
/// Integrand to precalculate the bound-free heating ratecoefficient in an approximate way
/// on a temperature grid using the assumption that T_e=T_R and W=1 in the ionisation
/// formula. The radiation fields dependence on W is taken into account by multiplying
/// the resulting expression with the correct W later on.
{
  double T  = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;

  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge,nu);

  return sigma_bf * (1-nu_edge/nu) * radfield::dbb(nu, T, 1) * exp(-HOVERKB*nu/T);
}*/

/*static double stimulated_recomb_integrand_gsl(double nu, void *paras)
/// Integrand to calculate the rate coefficient for spontaneous recombination
/// using gsl integrators.
{
  double T = ((gslintegration_paras *) paras)->T;
  double nu_edge = ((gslintegration_paras *) paras)->nu_edge;

  //double nu_edge = (epsilon(element,ion+1,0)-epsilon(element,ion,level))/H;
  /// Information about the current level is passed via the global variable
  /// mastate[tid] and its child values element, ion, level
  /// MAKE SURE THAT THESE ARE SET IN THE CALLING FUNCTION!!!!!!!!!!!!!!!!!
  double sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, nu_edge,nu);
  double x = sigma_bf / H / nu * radfield::dbb(nu,T,1) * exp(-HOVERKB*nu/T);
  ///in formula this looks like
  ///x = sigma_bf/H/nu * 2*H*pow(nu,3)/pow(CLIGHT,2) * exp(-H*nu/KB/T);

  ///set contributions from Lyman continuum artificially to zero to overcome it's large opacity
  return x;
}*/

static void precalculate_rate_coefficient_integrals(void) {
  // target fractional accuracy of the integrator //=1e-5 took 8 hours with Fe I to V!
  const double intaccuracy = RATECOEFF_INTEGRAL_ACCURACY;
  const double epsrelwarning = 1e-2;  // fractional error to emit a warning

  /// Calculate the rate coefficients for each level of each ion of each element
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element) - 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ion = 0; ion < nions; ion++) {
      // nlevels = get_nlevels(element,ion);
      const int atomic_number = get_element(element);
      const int ionstage = get_ionstage(element, ion);
      const int nlevels = get_ionisinglevels(element, ion);
      /// That's only an option for pure LTE
      // if (TAKE_N_BFCONTINUA < nlevels) nlevels = TAKE_N_BFCONTINUA;
      printout("Performing rate integrals for Z = %d, ion_stage %d...\n", atomic_number, ionstage);

      gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);

      for (int level = 0; level < nlevels; level++) {
        if ((level > 0) && (level % 50 == 0)) printout("  completed up to level %d of %d\n", level, nlevels);

        const int nphixstargets = get_nphixstargets(element, ion, level);
        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
          const int upperlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
          const double phixstargetprobability = get_phixsprobability(element, ion, level, phixstargetindex);

          // printout("element %d, ion %d, level %d, upperlevel %d, epsilon %g, continuum %g, nlevels
          // %d\n",element,ion,level,upperlevel,epsilon(element,ion,level),epsilon(element,ion+1,upperlevel),nlevels);

          // const double E_threshold = epsilon(element,ion+1,upperlevel) - epsilon(element,ion,level);
          const double E_threshold = get_phixs_threshold(element, ion, level, phixstargetindex);
          const double nu_threshold = E_threshold / H;
          const double nu_max_phixs =
              nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table
          // Loop over the temperature grid
          for (int iter = 0; iter < TABLESIZE; iter++) {
            double error;
            int status = 0;
            const float T_e = MINTEMP * exp(iter * T_step_log);
            // T_e = MINTEMP + iter*T_step;
            const double sfac = calculate_sahafact(element, ion, level, upperlevel, T_e, E_threshold);
            // printout("%d %g\n",iter,T_e);

            assert_always(globals::elements[element].ions[ion].levels[level].photoion_xs != NULL);
            // the threshold of the first target gives nu of the first phixstable point
            gslintegration_paras intparas = {
                .nu_edge = nu_threshold,
                .T = T_e,
                .photoion_xs = globals::elements[element].ions[ion].levels[level].photoion_xs};

            // gsl_function F_gamma;
            // F_gamma.function = &gamma_integrand_gsl;
            // F_gamma.params = &intparas;
            // gsl_function F_alpha_sp_E;
            // F_alpha_sp_E.function = &alpha_sp_E_integrand_gsl;
            // F_alpha_sp_E.params = &intparas;
            // F_stimulated_bfcooling.function = &stimulated_bfcooling_integrand_gsl;
            // F_stimulated_bfcooling.params = &intparas;
            // F_stimulated_recomb.function = &stimulated_recomb_integrand_gsl;
            // F_stimulated_recomb.params = &intparas;

            /// Spontaneous recombination and bf-cooling coefficient don't depend on the cutted radiation field
            double alpha_sp = 0.0;
            const gsl_function F_alpha_sp = {.function = &alpha_sp_integrand_gsl, .params = &intparas};

            status = gsl_integration_qag(&F_alpha_sp, nu_threshold, nu_max_phixs, 0, intaccuracy, GSLWSIZE,
                                         GSL_INTEG_GAUSS61, gslworkspace, &alpha_sp, &error);
            if (status != 0 && (status != 18 || (error / alpha_sp) > epsrelwarning)) {
              printout("alpha_sp integrator status %d. Integral value %9.3e +/- %9.3e\n", status, alpha_sp, error);
            }
            alpha_sp *= FOURPI * sfac * phixstargetprobability;

            if (!std::isfinite(alpha_sp) || alpha_sp < 0) {
              printout(
                  "WARNING: alpha_sp was negative or non-finite for level %d Te %g. alpha_sp %g sfac %g "
                  "phixstargetindex %d "
                  "phixstargetprobability %g\n",
                  level, T_e, alpha_sp, sfac, phixstargetindex, phixstargetprobability);
              alpha_sp = 0;
            }
            // assert_always(alpha_sp >= 0);
            globals::spontrecombcoeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] = alpha_sp;

            // if (atomic_number == 26 && ionstage == 3 && level < 5)
            // {
            //   const double E_threshold_b = get_phixs_threshold(element, ion, level, phixstargetindex);
            //   const double sfac_b = calculate_sahafact(element,ion,level,upperlevel,T_e,E_threshold_b);
            //   const double nu_threshold_b = E_threshold_b / H;
            //   const double nu_max_phixs_b = nu_threshold * last_phixs_nuovernuedge; //nu of the uppermost point in
            //   the phixs table intparas.nu_edge = nu_threshold_b;              // Global variable which passes the
            //   threshold to the integrator
            //                                                 // the threshold of the first target gives nu of the
            //                                                 first phixstable point
            //   double alpha_sp_new;
            //   status = gsl_integration_qag(&F_alpha_sp, nu_threshold_b, nu_max_phixs_b, 0, intaccuracy, GSLWSIZE,
            //   GSL_INTEG_GAUSS61, w, &alpha_sp_new, &error); alpha_sp_new *= FOURPI * sfac_b * phixstargetprobability;
            //   printout("recomb: T_e %6.1f Z=%d ionstage %d->%d upper+1 %5d lower+1 %5d sfac %7.2e sfac_new %7.2e
            //   alpha %7.2e alpha_new %7.2e threshold_ev %7.2e threshold_new_ev %7.2e\n",
            //           T_e, get_element(element), get_ionstage(element, ion + 1),
            //           get_ionstage(element, ion), upperlevel + 1, level + 1,
            //           sfac, sfac_b,
            //           alpha_sp, alpha_sp_new,
            //           E_threshold / EV,
            //           E_threshold_b / EV);
            // }

            // if (iter == 0)
            //   printout("alpha_sp: element %d ion %d level %d upper level %d at temperature %g, alpha_sp is %g
            //   (integral %g, sahafac %g)\n", element, ion, level, upperlevel, T_e, alpha_sp, alpha_sp/(FOURPI * sfac *
            //   phixstargetprobability),sfac);

#if (!NO_LUT_PHOTOION)
            double gammacorr = 0.0;
            const gsl_function F_gammacorr = {.function = &gammacorr_integrand_gsl, .params = &intparas};

            status = gsl_integration_qag(&F_gammacorr, nu_threshold, nu_max_phixs, 0, intaccuracy, GSLWSIZE,
                                         GSL_INTEG_GAUSS61, gslworkspace, &gammacorr, &error);
            if (status != 0 && (status != 18 || (error / gammacorr) > epsrelwarning)) {
              printout("gammacorr integrator status %d. Integral value %9.3e +/- %9.3e\n", status, gammacorr, error);
            }
            gammacorr *= FOURPI * phixstargetprobability;
            assert_always(gammacorr >= 0);
            if (gammacorr < 0) {
              printout("WARNING: gammacorr was negative for level %d\n", level);
              gammacorr = 0;
            }
            globals::corrphotoioncoeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] = gammacorr;
#endif

#if (!NO_LUT_BFHEATING)
            double bfheating_coeff = 0.0;
            const gsl_function F_bfheating = {.function = &approx_bfheating_integrand_gsl, .params = &intparas};

            status = gsl_integration_qag(&F_bfheating, nu_threshold, nu_max_phixs, 0, intaccuracy, GSLWSIZE,
                                         GSL_INTEG_GAUSS61, gslworkspace, &bfheating_coeff, &error);

            if (status != 0 && (status != 18 || (error / bfheating_coeff) > epsrelwarning)) {
              printout("bfheating_coeff integrator status %d. Integral value %9.3e +/- %9.3e\n", status,
                       bfheating_coeff, error);
            }
            bfheating_coeff *= FOURPI * phixstargetprobability;
            if (bfheating_coeff < 0) {
              printout("WARNING: bfheating_coeff was negative for level %d\n", level);
              bfheating_coeff = 0;
            }
            globals::bfheating_coeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] = bfheating_coeff;
#endif

            double bfcooling_coeff = 0.0;
            const gsl_function F_bfcooling = {.function = &bfcooling_integrand_gsl, .params = &intparas};

            status = gsl_integration_qag(&F_bfcooling, nu_threshold, nu_max_phixs, 0, intaccuracy, GSLWSIZE,
                                         GSL_INTEG_GAUSS61, gslworkspace, &bfcooling_coeff, &error);
            if (status != 0 && (status != 18 || (error / bfcooling_coeff) > epsrelwarning)) {
              printout("bfcooling_coeff integrator status %d. Integral value %9.3e +/- %9.3e\n", status,
                       bfcooling_coeff, error);
            }
            bfcooling_coeff *= FOURPI * sfac * phixstargetprobability;
            if (!std::isfinite(bfcooling_coeff) || bfcooling_coeff < 0) {
              printout(
                  "WARNING: bfcooling_coeff was negative or non-finite for level %d Te %g. bfcooling_coeff %g sfac %g "
                  "phixstargetindex %d phixstargetprobability %g\n",
                  level, T_e, bfcooling_coeff, sfac, phixstargetindex, phixstargetprobability);
              bfcooling_coeff = 0;
            }
            globals::bfcooling_coeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] = bfcooling_coeff;
          }
        }
      }
      gsl_set_error_handler(previous_handler);
    }
  }
}

double select_continuum_nu(int element, int lowerion, int lower, int upperionlevel, float T_e) {
  const int phixstargetindex = get_phixtargetindex(element, lowerion, lower, upperionlevel);
  const double E_threshold = get_phixs_threshold(element, lowerion, lower, phixstargetindex);
  const double nu_threshold = ONEOVERH * E_threshold;

  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

  const int npieces = globals::NPHIXSPOINTS;

  gslintegration_paras intparas = {.nu_edge = nu_threshold,
                                   .T = T_e,
                                   .photoion_xs = globals::elements[element].ions[lowerion].levels[lower].photoion_xs};

  const gsl_function F_alpha_sp = {.function = &alpha_sp_E_integrand_gsl, .params = &intparas};

  const double zrand = 1. - gsl_rng_uniform(rng);  // Make sure that 0 < zrand <= 1

  // printout("emitted bf photon Z=%2d ionstage %d->%d upper %4d lower %4d lambda %7.1f lambda_edge %7.1f ratio %g zrand
  // %g\n",
  //    get_element(element), get_ionstage(element, lowerion + 1), get_ionstage(element, lowerion), upperionlevel,
  //    lower, 1e8 * CLIGHT / nu_selected, 1e8 * CLIGHT / nu_threshold, nu_selected / nu_threshold, zrand);

  constexpr double intaccuracy = CONTINUUM_NU_INTEGRAL_ACCURACY;  /// Fractional accuracy of the integrator

  const double deltanu = (nu_max_phixs - nu_threshold) / npieces;
  double error;

  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);

  double total_alpha_sp = 0.;
  gsl_integration_qag(&F_alpha_sp, nu_threshold, nu_max_phixs, 0, intaccuracy, GSLWSIZE, GSL_INTEG_GAUSS31,
                      gslworkspace, &total_alpha_sp, &error);

  double alpha_sp_old = total_alpha_sp;
  double alpha_sp = total_alpha_sp;

  int i;
  for (i = 1; i < npieces; i++) {
    alpha_sp_old = alpha_sp;
    const double xlow = nu_threshold + i * deltanu;

    // Spontaneous recombination and bf-cooling coefficient don't depend on the cutted radiation field
    gsl_integration_qag(&F_alpha_sp, xlow, nu_max_phixs, 0, intaccuracy, GSLWSIZE, GSL_INTEG_GAUSS31, gslworkspace,
                        &alpha_sp, &error);

    if (zrand >= alpha_sp / total_alpha_sp) break;
  }

  gsl_set_error_handler(previous_handler);

  const double nuoffset = (total_alpha_sp * zrand - alpha_sp_old) / (alpha_sp - alpha_sp_old) * deltanu;
  const double nu_lower = nu_threshold + (i - 1) * deltanu + nuoffset;

  assert_testmodeonly(std::isfinite(nu_lower));

  return nu_lower;
}

__host__ __device__ double get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, float T_e)
/// Returns the rate coefficient for spontaneous recombination.
{
  /*int lowerindex = floor((T-MINTEMP)/T_step);
  int upperindex = lowerindex + 1;
  double T_upper =  MINTEMP + upperindex*T_step;
  double T_lower =  MINTEMP + lowerindex*T_step;*/
  double Alpha_sp;
  const int lowerindex = floor(log(T_e / MINTEMP) / T_step_log);
  assert_always(lowerindex >= 0);
  if (lowerindex < TABLESIZE - 1) {
    const int upperindex = lowerindex + 1;
    const double T_lower = MINTEMP * exp(lowerindex * T_step_log);
    const double T_upper = MINTEMP * exp(upperindex * T_step_log);

    const double f_upper = globals::spontrecombcoeff[get_bflutindex(upperindex, element, ion, level, phixstargetindex)];
    const double f_lower = globals::spontrecombcoeff[get_bflutindex(lowerindex, element, ion, level, phixstargetindex)];
    // printout("interpolate_spontrecombcoeff element %d, ion %d, level %d, upper %g, lower %g\n",
    //          element,ion,level,f_upper,f_lower);
    Alpha_sp = (f_lower + (f_upper - f_lower) / (T_upper - T_lower) * (T_e - T_lower));
  } else {
    Alpha_sp = globals::spontrecombcoeff[get_bflutindex(TABLESIZE - 1, element, ion, level, phixstargetindex)];
  }
  return Alpha_sp;
}

double calculate_ionrecombcoeff(const int modelgridindex, const float T_e, const int element, const int upperion,
                                const bool assume_lte, const bool collisional_not_radiative, const bool printdebug,
                                const bool lower_superlevel_only, const bool per_groundmultipletpop,
                                const bool stimonly)
// multiply by upper ion population (or ground population if per_groundmultipletpop is true) and nne to get a rate
{
  const int lowerion = upperion - 1;
  if (lowerion < 0) return 0.;

  double alpha = 0.;
  if (lowerion < get_nions(element) - 1) {
    // this gets divided and cancelled out in the radiative case anyway
    const double nne = (modelgridindex >= 0) ? grid::get_nne(modelgridindex) : 1.0;

    double nnupperion = 0;
    // nnupperion = get_groundmultiplet_pop(modelgridindex, T_e, element, upperion, assume_lte);
    int upper_nlevels;
    if (per_groundmultipletpop) {
      // assume that photoionisation of the ion below is only to the ground multiplet levels of the current ion
      // const int nphixstargets = get_nphixstargets(element, lowerion, 0);
      // upper_nlevels = get_phixsupperlevel(element, lowerion, 0, nphixstargets - 1) + 1;

      upper_nlevels = get_nlevels_groundterm(element, lowerion + 1);
    } else {
      upper_nlevels = get_nlevels(element, lowerion + 1);
    }

    for (int upper = 0; upper < upper_nlevels; upper++) {
      double nnupperlevel;
      if (assume_lte) {
        const double T_exc = T_e;
        const double E_level = epsilon(element, lowerion + 1, upper);
        const double E_ground = epsilon(element, lowerion + 1, 0);
        const double nnground = (modelgridindex >= 0) ? get_groundlevelpop(modelgridindex, element, lowerion + 1) : 1.0;

        nnupperlevel = (nnground * stat_weight(element, lowerion + 1, upper) / stat_weight(element, lowerion + 1, 0) *
                        exp(-(E_level - E_ground) / KB / T_exc));
      } else {
        nnupperlevel = get_levelpop(modelgridindex, element, lowerion + 1, upper);
      }
      nnupperion += nnupperlevel;
    }

    if (nnupperion <= 0.) {
      return 0.;
    }

    double nnupperlevel_so_far = 0.;
    const int maxrecombininglevel = get_maxrecombininglevel(element, lowerion + 1);
    for (int upper = 0; upper <= maxrecombininglevel; upper++) {
      double nnupperlevel;
      if (assume_lte) {
        const double T_exc = T_e;
        const double E_level = epsilon(element, lowerion + 1, upper);
        const double E_ground = epsilon(element, lowerion + 1, 0);
        const double nnground = (modelgridindex >= 0) ? get_groundlevelpop(modelgridindex, element, lowerion + 1) : 1.0;

        nnupperlevel = (nnground * stat_weight(element, lowerion + 1, upper) / stat_weight(element, lowerion + 1, 0) *
                        exp(-(E_level - E_ground) / KB / T_exc));
      } else {
        nnupperlevel = get_levelpop(modelgridindex, element, lowerion + 1, upper);
      }
      nnupperlevel_so_far += nnupperlevel;
      for (int lower = 0; lower < get_nlevels(element, lowerion); lower++) {
        if (lower_superlevel_only && (!level_isinsuperlevel(element, lowerion, lower))) continue;

        double recomb_coeff = 0.;
        if (collisional_not_radiative) {
          const double epsilon_trans = epsilon(element, lowerion + 1, upper) - epsilon(element, lowerion, lower);
          recomb_coeff += col_recombination_ratecoeff(modelgridindex, element, upperion, upper, lower, epsilon_trans);
        } else if (!stimonly) {
          recomb_coeff += rad_recombination_ratecoeff(T_e, nne, element, lowerion + 1, upper, lower, modelgridindex);
        } else {
          recomb_coeff += stim_recombination_ratecoeff(nne, element, upperion, upper, lower, modelgridindex);
        }

        const double alpha_level = recomb_coeff / nne;
        const double alpha_ion_contrib = alpha_level * nnupperlevel / nnupperion;
        alpha += alpha_ion_contrib;
        if (printdebug && alpha_ion_contrib > 0. && lower < 50) {
          printout(
              "recomb: Z=%d ionstage %d->%d upper+1 %5d lower+1 %5d alpha_level %7.2e alpha_ion_contrib %7.2e sum "
              "%7.2e nnlevel %7.2e nnionfrac %7.2e\n",
              get_element(element), get_ionstage(element, lowerion + 1), get_ionstage(element, lowerion), upper + 1,
              lower + 1, alpha_level, alpha_ion_contrib, alpha, nnupperlevel, nnupperlevel_so_far / nnupperion);
        }
      }
    }
  }
  if (printdebug) {
    printout("recomb: Z=%2d ionstage %d->%d upper+1 [all] lower+1 [all] Alpha %g\n\n", get_element(element),
             get_ionstage(element, lowerion + 1), get_ionstage(element, lowerion), alpha);
  }
  return alpha;
}

static void scale_level_phixs(const int element, const int ion, const int level, const double factor)
// multiply the cross sections associated with a level by some factor and
// also update the quantities integrated from (and proportional to) the cross sections
{
  // if we store the cross sections in node shared memory, then only one rank should update it
  if (globals::rank_in_node == 0) {
    for (int n = 0; n < globals::NPHIXSPOINTS; n++) {
      globals::elements[element].ions[ion].levels[level].photoion_xs[n] *= factor;
    }
  }

  const int nphixstargets = get_nphixstargets(element, ion, level);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
    for (int iter = 0; iter < TABLESIZE; iter++) {
      globals::spontrecombcoeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] *= factor;

#if (!NO_LUT_PHOTOION)
      globals::corrphotoioncoeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] *= factor;
#endif

#if (!NO_LUT_BFHEATING)
      globals::bfheating_coeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] *= factor;
#endif

      globals::bfcooling_coeff[get_bflutindex(iter, element, ion, level, phixstargetindex)] *= factor;
    }
  }
}

static void read_recombrate_file(void)
// calibrate the recombination rates to tabulated values by scaling the photoionisation cross sections
{
  use_cellhist = false;
  FILE *recombrate_file = fopen("recombrates.txt", "r");
  if (recombrate_file == NULL) {
    printout("No recombrates.txt file found. Skipping recombination rate scaling...\n");
    return;
  }

  printout("Reading recombination rate file (recombrates.txt)...\n");

  const double Te_estimate = RECOMBCALIBRATION_T_ELEC;
  const double log_Te_estimate = log10(Te_estimate);

  printout("Calibrating recombination rates for a temperature of %.1f K\n", Te_estimate);

  struct rrc_row {
    double log_Te;
    double rrc_low_n;
    double rrc_total;
  };

  int atomicnumber;
  int upperionstage;
  int tablerows;

  while (fscanf(recombrate_file, "%d %d %d\n", &atomicnumber, &upperionstage, &tablerows) > 0) {
    // printout("%d %d %d\n", atomicnumber, upperionstage, tablerows);

    struct rrc_row T_highestbelow = {0, 0, 0};
    struct rrc_row T_lowestabove = {0, 0, 0};
    T_highestbelow.log_Te = -1;
    T_lowestabove.log_Te = -1;
    for (int i = 0; i < tablerows; i++) {
      struct rrc_row row;
      assert_always(fscanf(recombrate_file, "%lg %lg %lg\n", &row.log_Te, &row.rrc_low_n, &row.rrc_total) == 3);
      if (row.log_Te < log_Te_estimate && row.log_Te > T_highestbelow.log_Te) {
        T_highestbelow.log_Te = row.log_Te;
        T_highestbelow.rrc_low_n = row.rrc_low_n;
        T_highestbelow.rrc_total = row.rrc_total;
      }

      if (row.log_Te > log_Te_estimate && (row.log_Te < T_lowestabove.log_Te || T_lowestabove.log_Te < 0)) {
        T_lowestabove.log_Te = row.log_Te;
        T_lowestabove.rrc_low_n = row.rrc_low_n;
        T_lowestabove.rrc_total = row.rrc_total;
      }
    }
    const int element = get_elementindex(atomicnumber);
    if (element >= 0) {
      const int ion = upperionstage - get_ionstage(element, 0);  // the index of the upper ion
      if (ion > 0 && ion < get_nions(element)) {
        printout("Z=%d ionstage %d->%d\n", atomicnumber, upperionstage, upperionstage - 1);
        assert_always(T_highestbelow.log_Te > 0);
        assert_always(T_lowestabove.log_Te > 0);

        const int nlevels = get_ionisinglevels(element, ion - 1);

        const double x = (log_Te_estimate - T_highestbelow.log_Te) / (T_lowestabove.log_Te - T_highestbelow.log_Te);
        const double input_rrc_low_n = x * T_highestbelow.rrc_low_n + (1 - x) * T_lowestabove.rrc_low_n;
        const double input_rrc_total = x * T_highestbelow.rrc_total + (1 - x) * T_lowestabove.rrc_total;

        const bool assume_lte = true;
        const bool printdebug = false;
        const bool per_groundmultipletpop = true;

        double rrc = calculate_ionrecombcoeff(-1, Te_estimate, element, ion, assume_lte, false, printdebug, false,
                                              per_groundmultipletpop, false);
        printout("              rrc: %10.3e\n", rrc);

        if (input_rrc_low_n >= 0)  // if it's < 0, ignore it
        {
          printout("  input_rrc_low_n: %10.3e\n", input_rrc_low_n);

          const double phixs_multiplier = input_rrc_low_n / rrc;
          if (phixs_multiplier < 0.05 || phixs_multiplier >= 2.0) {
            printout("    Not scaling phixs of all levels by %.3f (because < 0.05 or >= 2.0)\n", phixs_multiplier);
          } else {
            printout("    scaling phixs of all levels by %.3f\n", phixs_multiplier);

            for (int level = 0; level < nlevels; level++) scale_level_phixs(element, ion - 1, level, phixs_multiplier);

            rrc = calculate_ionrecombcoeff(-1, Te_estimate, element, ion, assume_lte, false, printdebug, false,
                                           per_groundmultipletpop, false);
            printout("              rrc: %10.3e\n", rrc);
          }
        }

        // hopefully the RRC now matches the low_n value well, if it was defined
        // Next, use the superlevel recombination rates to make up the excess needed to reach the total RRC

        printout("  input_rrc_total: %10.3e\n", input_rrc_total);

        if (rrc < input_rrc_total) {
          const double rrc_superlevel = calculate_ionrecombcoeff(-1, Te_estimate, element, ion, assume_lte, false,
                                                                 printdebug, true, per_groundmultipletpop, false);
          printout("  rrc(superlevel): %10.3e\n", rrc_superlevel);

          if (rrc_superlevel > 0) {
            const double phixs_multiplier_superlevel = 1.0 + (input_rrc_total - rrc) / rrc_superlevel;
            printout("    scaling phixs of levels in the superlevel by %.3f\n", phixs_multiplier_superlevel);
            assert_always(phixs_multiplier_superlevel >= 0);

            const int first_superlevel_level = get_nlevels_nlte(element, ion - 1) + 1;
            for (int level = first_superlevel_level; level < nlevels; level++)
              scale_level_phixs(element, ion - 1, level, phixs_multiplier_superlevel);
          } else {
            printout("There is no superlevel recombination, so multiplying all levels instead\n");
            const double phixs_multiplier = input_rrc_total / rrc;
            printout("    scaling phixs of all levels by %.3f\n", phixs_multiplier);
            assert_always(phixs_multiplier >= 0);

            for (int level = 0; level < nlevels; level++) scale_level_phixs(element, ion - 1, level, phixs_multiplier);
          }
        } else {
          printout("rrc >= input_rrc_total!\n");
          const double phixs_multiplier = input_rrc_total / rrc;
          printout("    scaling phixs of all levels by %.3f\n", phixs_multiplier);
          assert_always(phixs_multiplier >= 0);

          for (int level = 0; level < nlevels; level++) scale_level_phixs(element, ion - 1, level, phixs_multiplier);
        }

        rrc = calculate_ionrecombcoeff(-1, Te_estimate, element, ion, assume_lte, false, printdebug, false,
                                       per_groundmultipletpop, false);
        printout("              rrc: %10.3e\n", rrc);
      }
    }
  }
  fclose(recombrate_file);
}

static void precalculate_ion_alpha_sp(void) {
  for (int iter = 0; iter < TABLESIZE; iter++) {
    const float T_e = MINTEMP * exp(iter * T_step_log);
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element) - 1;
      for (int ion = 0; ion < nions; ion++) {
        // nlevels = get_nlevels(element, ion);
        // nlevels = get_ionisinglevels(element, ion); ///number of levels of the current ion which have an associated
        // photoion cross section
        const int nlevels = get_ionisinglevels(element, ion);  /// number of ionising levels used in the simulation
        double zeta = 0.;
        for (int level = 0; level < nlevels; level++) {
          for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level);
               phixstargetindex++) {
            const double zeta_level = get_spontrecombcoeff(element, ion, level, phixstargetindex, T_e);
            zeta += zeta_level;
          }
        }
        globals::elements[element].ions[ion].Alpha_sp[iter] = zeta;
      }
    }
  }
}

void ratecoefficients_init(void)
/// Precalculates the rate coefficients for stimulated and spontaneous
/// recombination and photoionisation on a given temperature grid using
/// libgsl integrators.
/// NB: with the nebular approximation they only depend on T_e, T_R and W.
/// W is easily factored out. For stimulated recombination we must assume
/// T_e = T_R for this precalculation.
{
  /// Determine the temperture grids gridsize
  T_step = (1. * MAXTEMP - MINTEMP) / (TABLESIZE - 1.);  /// global variables
  T_step_log = (log(MAXTEMP) - log(MINTEMP)) / (TABLESIZE - 1.);

  md5_file("adata.txt", adatafile_hash);
  md5_file("compositiondata.txt", compositionfile_hash);
  assert_always(phixs_file_version >= 0);  // check that it has been changed from the default value of -1
  md5_file(phixsdata_filenames[phixs_file_version], phixsfile_hash);

  /// Check if we need to calculate the ratecoefficients or if we were able to read them from file
  if (!read_ratecoeff_dat()) {
    precalculate_rate_coefficient_integrals();
    /// And the master process writes them to file in a serial operation
    if (globals::rank_global == 0) {
      write_ratecoeff_dat();
    }
  }

  read_recombrate_file();

  precalculate_ion_alpha_sp();
}

#if (!NO_LUT_PHOTOION)
double interpolate_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, double T) {
  const int lowerindex = floor(log(T / MINTEMP) / T_step_log);
  if (lowerindex < TABLESIZE - 1) {
    const int upperindex = lowerindex + 1;
    const double T_lower = MINTEMP * exp(lowerindex * T_step_log);
    const double T_upper = MINTEMP * exp(upperindex * T_step_log);

    const double f_upper =
        globals::corrphotoioncoeff[get_bflutindex(upperindex, element, ion, level, phixstargetindex)];
    const double f_lower =
        globals::corrphotoioncoeff[get_bflutindex(lowerindex, element, ion, level, phixstargetindex)];

    return (f_lower + (f_upper - f_lower) / (T_upper - T_lower) * (T - T_lower));
  } else
    return globals::corrphotoioncoeff[get_bflutindex(TABLESIZE - 1, element, ion, level, phixstargetindex)];
}

double get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex)
/// Returns the for stimulated emission corrected photoionisation rate coefficient.
{
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  const double W = grid::get_W(modelgridindex);
  const double T_R = grid::get_TR(modelgridindex);

  return W * interpolate_corrphotoioncoeff(element, ion, level, phixstargetindex, T_R);
}
#endif

static double integrand_stimrecombination_custom_radfield(const double nu, void *voidparas) {
  {
    const gsl_integral_paras_gammacorr *const params = (gsl_integral_paras_gammacorr *)voidparas;
    const int modelgridindex = params->modelgridindex;
    const float T_e = params->T_e;

    const float sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, params->nu_edge, nu);

    const double Jnu = radfield::radfield(nu, modelgridindex);

    // TODO: MK thesis page 41, use population ratios and Te?
    return ONEOVERH * sigma_bf / nu * Jnu * exp(-HOVERKB * nu / T_e);
  }
}

static double calculate_stimrecombcoeff_integral(int element, int lowerion, int level, int phixstargetindex,
                                                 int modelgridindex) {
  // if (nnlevel <= 1.1 * MINPOP)
  // {
  //   return 0.;
  // }

  const double epsrel = 1e-3;
  const double epsabs = 0.;

  const double E_threshold = get_phixs_threshold(element, lowerion, level, phixstargetindex);
  const double nu_threshold = ONEOVERH * E_threshold;
  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

  const float T_e = grid::get_Te(modelgridindex);
  gsl_integral_paras_gammacorr intparas = {
      .nu_edge = nu_threshold,
      .photoion_xs = globals::elements[element].ions[lowerion].levels[level].photoion_xs,
      .T_e = T_e,
      .modelgridindex = modelgridindex,
  };

  const int upperionlevel = get_phixsupperlevel(element, lowerion, level, phixstargetindex);
  const double sf = calculate_sahafact(element, lowerion, level, upperionlevel, T_e, H * nu_threshold);

  const gsl_function F_stimrecomb = {.function = &integrand_stimrecombination_custom_radfield, .params = &intparas};
  double error = 0.0;

  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);
  double stimrecombcoeff = 0.0;

  // const int status =
  gsl_integration_qag(&F_stimrecomb, nu_threshold, nu_max_phixs, epsabs, epsrel, GSLWSIZE, GSL_INTEG_GAUSS61,
                      gslworkspace, &stimrecombcoeff, &error);

  gsl_set_error_handler(previous_handler);

  stimrecombcoeff *= FOURPI * sf * get_phixsprobability(element, lowerion, level, phixstargetindex);

  // if (status != 0)
  // {
  //   error *= FOURPI * get_phixsprobability(element, ion, level, phixstargetindex);
  //   printout("stimrecombcoeff gsl integrator warning %d. modelgridindex %d Z=%d ionstage %d lower %d phixstargetindex
  //   %d gamma %g error %g\n",
  //            status, modelgridindex, get_element(element), get_ionstage(element, ion), level, phixstargetindex,
  //            gammacorr, error);
  // }

  return stimrecombcoeff;
}

double get_stimrecombcoeff(int element, int lowerion, int level, int phixstargetindex, int modelgridindex)
/// Returns the stimulated recombination rate coefficient
// multiple by upper level population and nne to get rate
{
  double stimrecombcoeff = -1.;
#if (SEPARATE_STIMRECOMB)
  if (use_cellhist) {
    stimrecombcoeff = globals::cellhistory[tid]
                          .chelements[element]
                          .chions[lowerion]
                          .chlevels[level]
                          .chphixstargets[phixstargetindex]
                          .stimrecombcoeff;
  }
#endif

  if (!use_cellhist || stimrecombcoeff < 0) {
#ifdef FORCE_LTE
    stimrecombcoeff = 0.;
#else
    stimrecombcoeff = calculate_stimrecombcoeff_integral(element, lowerion, level, phixstargetindex, modelgridindex);
#endif

#if (SEPARATE_STIMRECOMB)
    if (use_cellhist) {
      globals::cellhistory[tid]
          .chelements[element]
          .chions[lowerion]
          .chlevels[level]
          .chphixstargets[phixstargetindex]
          .stimrecombcoeff = stimrecombcoeff;
    }
#endif
  }

  return stimrecombcoeff;
}

static double integrand_corrphotoioncoeff_custom_radfield(const double nu, void *const voidparas)
/// Integrand to calculate the rate coefficient for photoionization
/// using gsl integrators. Corrected for stimulated recombination.
{
  const gsl_integral_paras_gammacorr *const params = (gsl_integral_paras_gammacorr *)voidparas;
  const int modelgridindex = params->modelgridindex;

#if (SEPARATE_STIMRECOMB)
  const double corrfactor = 1.0;
#else
  const float T_e = params->T_e;
  double corrfactor = 1. - params->departure_ratio * exp(-HOVERKB * nu / T_e);
  if (corrfactor < 0) corrfactor = 0.;
#endif

  const float sigma_bf = photoionization_crosssection_fromtable(params->photoion_xs, params->nu_edge, nu);

  // const double Jnu = use_cellhist ? radfield::radfield(nu, modelgridindex) : radfield::dbb_mgi(nu, modelgridindex);
  const double Jnu = radfield::radfield(nu, modelgridindex);

  // TODO: MK thesis page 41, use population ratios and Te?
  return ONEOVERH * sigma_bf / nu * Jnu * corrfactor;
}

static double calculate_corrphotoioncoeff_integral(int element, int ion, int level, int phixstargetindex,
                                                   int modelgridindex) {
  constexpr double epsrel = 1e-3;
  constexpr double epsrelwarning = 1e-1;
  constexpr double epsabs = 0.;

  const double E_threshold = get_phixs_threshold(element, ion, level, phixstargetindex);
  const double nu_threshold = ONEOVERH * E_threshold;
  const double nu_max_phixs = nu_threshold * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

  const float T_e = grid::get_Te(modelgridindex);

#if SEPARATE_STIMRECOMB
  const double departure_ratio = 0.;  // zero the stimulated recomb contribution
#else
  // stimulated recombination is negative photoionisation
  const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
  // if (nnlevel <= 1.1 * MINPOP)
  // {
  //   return 0.;
  // }
  const double nne = grid::get_nne(modelgridindex);
  const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
  const double sf = calculate_sahafact(element, ion, level, upperionlevel, T_e, H * nu_threshold);
  const double nnupperionlevel = get_levelpop(modelgridindex, element, ion + 1, upperionlevel);
  double departure_ratio = nnlevel > 0. ? nnupperionlevel / nnlevel * nne * sf : 1.0;  // put that to phixslist
  if (!std::isfinite(departure_ratio)) {
    departure_ratio = 0.;
  }
#endif
  gsl_integral_paras_gammacorr intparas = {
      .nu_edge = nu_threshold,
      .departure_ratio = departure_ratio,
      .photoion_xs = globals::elements[element].ions[ion].levels[level].photoion_xs,
      .T_e = T_e,
      .modelgridindex = modelgridindex,
  };

  const gsl_function F_gammacorr = {.function = &integrand_corrphotoioncoeff_custom_radfield, .params = &intparas};
  double error = 0.0;

  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);

  double gammacorr = 0.0;
  const int status = gsl_integration_qag(&F_gammacorr, nu_threshold, nu_max_phixs, epsabs, epsrel, GSLWSIZE,
                                         GSL_INTEG_GAUSS61, gslworkspace, &gammacorr, &error);

  gsl_set_error_handler(previous_handler);

  if (status != 0 && (status != 18 || (error / gammacorr) > epsrelwarning)) {
    printout(
        "corrphotoioncoeff gsl integrator warning %d. modelgridindex %d Z=%d ionstage %d lower %d phixstargetindex %d "
        "integral %g error %g\n",
        status, modelgridindex, get_element(element), get_ionstage(element, ion), level, phixstargetindex, gammacorr,
        error);
    if (!std::isfinite(gammacorr)) gammacorr = 0.;
  }

  gammacorr *= FOURPI * get_phixsprobability(element, ion, level, phixstargetindex);

  return gammacorr;
}

double get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex)
/// Returns the photoionisation rate coefficient (corrected for stimulated emission)
{
  /// The correction factor for stimulated emission in gammacorr is set to its
  /// LTE value. Because the T_e dependence of gammacorr is weak, this correction
  /// correction may be evaluated at T_R!
  double gammacorr = -1;

  if (DETAILED_BF_ESTIMATORS_ON && globals::nts_global >= DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP) {
    gammacorr = radfield::get_bfrate_estimator(element, ion, level, phixstargetindex, modelgridindex);
    // gammacorr will be -1 if no estimators available
    if (gammacorr > 0) {
      return gammacorr;
    }
  }

  if (use_cellhist) {
    gammacorr = globals::cellhistory[tid]
                    .chelements[element]
                    .chions[ion]
                    .chlevels[level]
                    .chphixstargets[phixstargetindex]
                    .corrphotoioncoeff;
  }

  if (!use_cellhist || gammacorr < 0) {
#ifdef FORCE_LTE
    {
      /// Interpolate gammacorr out of precalculated values
      const double T_R = grid::get_TR(modelgridindex);
      gammacorr = interpolate_corrphotoioncoeff(element, ion, level, phixstargetindex, T_R);
    }
#else
    {
#if (NO_LUT_PHOTOION)
      gammacorr = calculate_corrphotoioncoeff_integral(element, ion, level, phixstargetindex, modelgridindex);
#else
      const double W = grid::get_W(modelgridindex);
      const double T_R = grid::get_TR(modelgridindex);

      gammacorr = W * interpolate_corrphotoioncoeff(element, ion, level, phixstargetindex, T_R);
      const int index_in_groundlevelcontestimator =
          globals::elements[element].ions[ion].levels[level].closestgroundlevelcont;
      if (index_in_groundlevelcontestimator >= 0) {
        gammacorr *= globals::corrphotoionrenorm[modelgridindex * get_nelements() * get_max_nions() +
                                                 index_in_groundlevelcontestimator];
      }
#endif
    }
#endif
    if (use_cellhist) {
      globals::cellhistory[tid]
          .chelements[element]
          .chions[ion]
          .chlevels[level]
          .chphixstargets[phixstargetindex]
          .corrphotoioncoeff = gammacorr;
    }
  }

  return gammacorr;
}

static int get_nlevels_important(int modelgridindex, int element, int ion, bool assume_lte, float T_e,
                                 double *nnlevelsum_out)
// get the number of levels that make up a fraction of the ion population
// of at least IONGAMMA_POPFRAC_LEVELS_INCLUDED
{
  if (IONGAMMA_POPFRAC_LEVELS_INCLUDED >= 1.) {
    return get_nlevels(element, ion);
  }
  // get the stored ion population for comparison with the cumulative sum of level pops
  const double nnion_real = ionstagepop(modelgridindex, element, ion);

  double nnlevelsum = 0.;
  int nlevels_important = get_ionisinglevels(element, ion);  // levels needed to get majority of ion pop

  // debug: treat all ionising levels as important
  // *nnlevelsum_out = nnion_real;
  // return nlevels_important;

  for (int lower = 0;
       (nnlevelsum / nnion_real < IONGAMMA_POPFRAC_LEVELS_INCLUDED) && (lower < get_ionisinglevels(element, ion));
       lower++) {
    double nnlowerlevel;
    if (assume_lte) {
      const double T_exc = T_e;  // remember, other parts of the code in LTE mode use TJ, not T_e
      const double E_level = epsilon(element, ion, lower);
      const double E_ground = epsilon(element, ion, 0);
      const double nnground = (modelgridindex >= 0) ? get_groundlevelpop(modelgridindex, element, ion) : 1.0;

      nnlowerlevel = (nnground * stat_weight(element, ion, lower) / stat_weight(element, ion, 0) *
                      exp(-(E_level - E_ground) / KB / T_exc));
    } else {
      nnlowerlevel = get_levelpop(modelgridindex, element, ion, lower);
    }
    nnlevelsum += nnlowerlevel;
    nlevels_important = lower + 1;
  }
  *nnlevelsum_out = nnlevelsum;
  assert_always(nlevels_important <= get_nlevels(element, ion));
  // printout("mgi %d element %d ion %d nlevels_important %d popfrac %g\n", modelgridindex, element, ion,
  // nlevels_important, nnlevelsum / nnion_real);
  return nlevels_important;
}

double calculate_iongamma_per_gspop(const int modelgridindex, const int element, const int ion)
// ionisation rate coefficient. multiply by get_groundlevelpop to get a rate [s^-1]
{
  const int nions = get_nions(element);
  double Gamma = 0.;
  if (ion >= nions - 1) {
    return 0.;
  }

  const float T_e = grid::get_Te(modelgridindex);
  const float nne = grid::get_nne(modelgridindex);

  // double nnlowerion = 0.;
  // const int nlevels_important = get_nlevels_important(modelgridindex, element, ion, false, T_e, &nnlowerion);
  const int nlevels_important = get_nlevels(element, ion);

  double Col_ion = 0.;
  for (int level = 0; level < nlevels_important; level++) {
    const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
    const int nphixstargets = get_nphixstargets(element, ion, level);
    for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
      const int upperlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);

      Gamma += nnlevel * get_corrphotoioncoeff(element, ion, level, phixstargetindex, modelgridindex);

      const double epsilon_trans = epsilon(element, ion + 1, upperlevel) - epsilon(element, ion, level);
      // printout("%g %g %g\n", get_levelpop(n,element,ion,level),col_ionization(n,0,epsilon_trans),epsilon_trans);

      Col_ion += nnlevel * col_ionization_ratecoeff(T_e, nne, element, ion, level, phixstargetindex, epsilon_trans);
    }
  }
  // printout("element %d ion %d: col/gamma %g Te %g ne %g\n", element, ion, Col_ion/Gamma, grid::get_Te(n),
  // grid::get_nne(n));
  Gamma += Col_ion;
  Gamma /= get_groundlevelpop(modelgridindex, element, ion);
  return Gamma;
}

double calculate_iongamma_per_ionpop(const int modelgridindex, const float T_e, const int element, const int lowerion,
                                     const bool assume_lte, const bool collisional_not_radiative, const bool printdebug,
                                     const bool force_bfest, const bool force_bfintegral)
// ionisation rate coefficient. multiply by the lower ion pop to get a rate
// currently only used for the estimator output file, not the simulation
{
  assert_always(lowerion < get_nions(element) - 1);
  assert_always(!force_bfest || !force_bfintegral);

  const float nne = (modelgridindex >= 0) ? grid::get_nne(modelgridindex) : 1.0;

  double nnlowerion = 0.;
  const int nlevels_important = get_nlevels_important(modelgridindex, element, lowerion, assume_lte, T_e, &nnlowerion);

  if (nnlowerion <= 0.) {
    return 0.;
  }

  double gamma_ion = 0.;
  double gamma_ion_used = 0.;
  for (int lower = 0; lower < nlevels_important; lower++) {
    double nnlowerlevel;
    if (assume_lte) {
      const double T_exc = T_e;
      const double E_level = epsilon(element, lowerion, lower);
      const double E_ground = epsilon(element, lowerion, 0);
      const double nnground = get_groundlevelpop(modelgridindex, element, lowerion);

      nnlowerlevel = (nnground * stat_weight(element, lowerion, lower) / stat_weight(element, lowerion, 0) *
                      exp(-(E_level - E_ground) / KB / T_exc));
    } else {
      nnlowerlevel = get_levelpop(modelgridindex, element, lowerion, lower);
    }

    for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, lowerion, lower); phixstargetindex++) {
      const int upper = get_phixsupperlevel(element, lowerion, lower, phixstargetindex);

      double gamma_coeff_integral = 0.;
      double gamma_coeff_bfest = 0.;
      double gamma_coeff_used = 0.;
      if (collisional_not_radiative) {
        const double epsilon_trans = epsilon(element, lowerion + 1, upper) - epsilon(element, lowerion, lower);
        gamma_coeff_used +=
            col_ionization_ratecoeff(T_e, nne, element, lowerion, lower, phixstargetindex, epsilon_trans);
      } else {
        gamma_coeff_used = get_corrphotoioncoeff(element, lowerion, lower, phixstargetindex,
                                                 modelgridindex);  // whatever ARTIS uses internally

        if (force_bfest || printdebug) {
          gamma_coeff_bfest =
              radfield::get_bfrate_estimator(element, lowerion, lower, phixstargetindex, modelgridindex);
        }

        if (force_bfintegral || printdebug) {
          // use the cellhistory but not the detailed bf estimators
          // TODO: restore cell history part
          gamma_coeff_integral +=
              calculate_corrphotoioncoeff_integral(element, lowerion, lower, phixstargetindex, modelgridindex);
          // double gamma_coeff_integral_level_ch = globals::cellhistory[tid]
          //                                            .chelements[element]
          //                                            .chions[lowerion]
          //                                            .chlevels[lower]
          //                                            .chphixstargets[phixstargetindex]
          //                                            .corrphotoioncoeff;
          // if (gamma_coeff_integral_level_ch >= 0) {
          //   gamma_coeff_integral += gamma_coeff_integral_level_ch;
          // } else {
          //   gamma_coeff_integral +=
          //       calculate_corrphotoioncoeff_integral(element, lowerion, lower, phixstargetindex, modelgridindex);
          // }
        }
      }

      const double gamma_ion_contribution_used = gamma_coeff_used * nnlowerlevel / nnlowerion;
      const double gamma_ion_contribution_bfest = gamma_coeff_bfest * nnlowerlevel / nnlowerion;
      const double gamma_ion_contribution_integral = gamma_coeff_integral * nnlowerlevel / nnlowerion;
      gamma_ion_used += gamma_ion_contribution_used;
      if (force_bfest) {
        gamma_ion += gamma_ion_contribution_bfest;
      } else if (force_bfintegral) {
        gamma_ion += gamma_ion_contribution_integral;
      } else {
        gamma_ion += gamma_ion_contribution_used;
      }

      if (printdebug && (gamma_ion_contribution_integral > 0. || gamma_ion_contribution_used > 0.) && lower < 20) {
        const double threshold_angstroms =
            1e8 * CLIGHT / (get_phixs_threshold(element, lowerion, lower, phixstargetindex) / H);
        printout(
            "Gamma_R: Z=%d ionstage %d->%d lower+1 %5d upper+1 %5d lambda_threshold %7.1f Gamma_integral %7.2e "
            "Gamma_bfest %7.2e Gamma_used %7.2e Gamma_used_sum %7.2e\n",
            get_element(element), get_ionstage(element, lowerion), get_ionstage(element, lowerion + 1), lower + 1,
            upper + 1, threshold_angstroms, gamma_ion_contribution_integral, gamma_ion_contribution_bfest,
            gamma_ion_contribution_used, gamma_ion_used);

#if (DETAILED_BF_ESTIMATORS_BYTYPE)
        radfield::print_bfrate_contributions(element, lowerion, lower, phixstargetindex, modelgridindex, nnlowerlevel,
                                             nnlowerion);
#endif
      }
    }
  }
  if (printdebug) {
    printout("Gamma_R: Z=%d ionstage %d->%d lower+1 [all] upper+1 [all] Gamma_used_ion %7.2e\n", get_element(element),
             get_ionstage(element, lowerion), get_ionstage(element, lowerion + 1), gamma_ion_used);
  }

  return gamma_ion;
}
