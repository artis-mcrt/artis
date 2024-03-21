#pragma once
#ifndef RATECOEFF_H
#define RATECOEFF_H

#include "sn3d.h"

void ratecoefficients_init();

void setup_photoion_luts();

[[nodiscard]] auto select_continuum_nu(int element, int lowerion, int lower, int upperionlevel, float T_e) -> double;

[[nodiscard]] auto interpolate_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex,
                                                 double T) -> double;

[[nodiscard]] auto get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, float T_e) -> double;
[[nodiscard]] auto get_stimrecombcoeff(int element, int lowerion, int level, int phixstargetindex,
                                       int modelgridindex) -> double;

[[nodiscard]] auto get_bfcoolingcoeff(int element, int ion, int level, int phixstargetindex, float T_e) -> double;

[[nodiscard]] auto get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex,
                                         int modelgridindex) -> double;
[[nodiscard]] auto get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex,
                                             int modelgridindex) -> double;

[[nodiscard]] auto iongamma_is_zero(int nonemptymgi, int element, int ion) -> bool;

[[nodiscard]] auto calculate_iongamma_per_gspop(int modelgridindex, int element, int ion) -> double;
[[nodiscard]] auto calculate_iongamma_per_ionpop(int modelgridindex, float T_e, int element, int lowerion,
                                                 bool assume_lte, bool collisional_not_radiative, bool printdebug,
                                                 bool force_bfest, bool force_bfintegral) -> double;

[[nodiscard]] auto calculate_ionrecombcoeff(int modelgridindex, float T_e, int element, int upperion, bool assume_lte,
                                            bool collisional_not_radiative, bool printdebug, bool lower_superlevel_only,
                                            bool per_groundmultipletpop, bool stimonly) -> double;

extern double T_step_log;

#if defined(STDPAR_ON) && defined(GPU_ON)
template <double func_integrand(double, void *)>
constexpr auto simpson_integrator(auto &params, double a, double b, int samplecount) -> double {
  assert_always(samplecount % 2 == 1);

  const double deltax = (b - a) / samplecount;

  double integral = 0.;
  for (int i = 0; i < samplecount; i++) {
    // Simpson's rule integral (will later be divided by 3)
    // n must be odd
    // integral = (xn - x0) / 3 * {f(x_0) + 4 * f(x_1) + 2 * f(x_2) + ... + 4 * f(x_1) + f(x_n-1)}
    // weights e.g., 1,4,2,4,2,4,1
    double weight{1.};
    if (i == 0 || i == (samplecount - 1)) {
      weight = 1.;
    } else if (i % 2 == 0) {
      weight = 2.;
    } else {
      weight = 4.;
    }

    const double x = a + deltax * i;

    integral += weight * func_integrand(x, &params) * deltax;
  }
  integral /= 3.;

  return integral;
}
#else
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include "constants.h"
#endif

template <double func_integrand(double, void *)>
auto integrator(auto params, double a, double b, double epsabs, double epsrel, int key, double *result,
                double *abserr) {
#if defined(STDPAR_ON) && defined(GPU_ON)
  const int samplecount = globals::NPHIXSPOINTS * 16 + 1;  // need an odd number for Simpson rule
  *result = simpson_integrator<func_integrand>(params, a, b, samplecount);
  *abserr = 0.;
  return 0;
#else
  const gsl_function F = {.function = (func_integrand), .params = &(params)};
  return gsl_integration_qag(&F, a, b, epsabs, epsrel, GSLWSIZE, key, gslworkspace.get(), result, abserr);
#endif
}

#endif  // RATECOEFF_H