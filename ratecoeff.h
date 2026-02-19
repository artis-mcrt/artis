#ifndef RATECOEFF_H
#define RATECOEFF_H

#include "constants.h"
#include "sn3d.h"

#ifdef USE_SIMPSON_INTEGRATOR
#include <algorithm>

#include "globals.h"

#elifdef BOOST_OFF

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include <cstddef>
#include <memory>

constexpr size_t GSLWSIZE = 16384;  // GSL integration workspace size
inline thread_local auto gslworkspace =
    std::unique_ptr<gsl_integration_workspace, void (*)(gsl_integration_workspace*)>{
        gsl_integration_workspace_alloc(GSLWSIZE),
        [](gsl_integration_workspace* const w) -> void { gsl_integration_workspace_free(w); }};

inline void gsl_error_handler_printout(const char* reason, const char* file, int line, int gsl_errno) {
  if (gsl_errno != GSL_EROUND) {
    printlnlog("WARNING: gsl ({}:{}): {} (Error code {})", file, line, reason, gsl_errno);
  }
}

#else
#include <cstdlib>

#pragma clang unsafe_buffer_usage begin
#include <boost/math/quadrature/gauss_kronrod.hpp>
#pragma clang unsafe_buffer_usage end

#endif

void ratecoefficients_init();

void setup_photoion_luts();

[[nodiscard]] __host__ __device__ auto select_continuum_nu(int element, int lowerion, int lower, int upperionlevel,
                                                           float T_e) -> double;
[[nodiscard]] __host__ __device__ auto get_spontrecombcoeff(int uniquelevelindex, int phixstargetindex, float T_e)
    -> double;
[[nodiscard]] auto get_stimrecombcoeff(int element, int lowerion, int level, int phixstargetindex, int nonemptymgi)
    -> double;

[[nodiscard]] __host__ __device__ auto get_bfcoolingcoeff(int element, int lowerion, int lowerionlevel,
                                                          int phixstargetindex, float T_e) -> double;
[[nodiscard]] __host__ __device__ auto get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex,
                                                              double T_R, double W) -> double;

[[nodiscard]] __host__ __device__ auto get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex,
                                                             int nonemptymgi) -> double;
[[nodiscard]] auto get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, int nonemptymgi)
    -> double;

[[nodiscard]] auto iongamma_is_zero(int nonemptymgi, int element, int ion) -> bool;

[[nodiscard]] auto calculate_iongamma_per_gspop(int nonemptymgi, int element, int ion) -> double;
[[nodiscard]] auto calculate_iongamma_per_ionpop(int nonemptymgi, float T_e, int element, int lowerion, bool assume_lte,
                                                 bool collisional_not_radiative, bool printdebug, bool force_bfest,
                                                 bool force_bfintegral) -> double;

[[nodiscard]] auto calculate_ionrecombcoeff(int nonemptymgi, float T_e, int element, int upperion, bool assume_lte,
                                            bool collisional_not_radiative, bool lower_superlevel_only,
                                            bool per_groundmultipletpop) -> double;

[[gnu::pure]] [[nodiscard]] __host__ __device__ auto get_ion_spontrecombcoeff(int uniqueionindex, float T_e) -> double;

template <double func_integrand(double, void* const)>
constexpr auto simpson_integrator(auto& params, const double a, const double b, const int samplecount) -> double {
  assert_testmodeonly(samplecount % 2 == 1);

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

    const double x = a + (deltax * i);

    integral += weight * func_integrand(x, &params) * deltax;
  }
  integral /= 3.;

  return integral;
}

template <double func_integrand(double, void* const), int GKNPOINTS = 61>
auto integrator(auto params, const double a, const double b, const double epsrel, double* abserr) -> double {
  static_assert(GKNPOINTS == 15 || GKNPOINTS == 31 || GKNPOINTS == 41 || GKNPOINTS == 51 || GKNPOINTS == 61,
                "Unsupported GKNPOINTS value");
  double result{0.};
#ifdef USE_SIMPSON_INTEGRATOR
  // need an odd number for Simpson rule
  const int samplecount = (std::max(1, static_cast<int>((b / a) / globals::NPHIXSNUINCREMENT)) * 4) + 1;

  result = simpson_integrator<func_integrand>(params, a, b, samplecount);
  *abserr = 0.;

#else
#ifdef BOOST_OFF

  // GSL's QAG adaptive integrator
  constexpr auto key = []() {
    switch (GKNPOINTS) {
      case 15:
        return GSL_INTEG_GAUSS15;
      case 21:
        return GSL_INTEG_GAUSS21;
      case 31:
        return GSL_INTEG_GAUSS31;
      case 41:
        return GSL_INTEG_GAUSS41;
      case 51:
        return GSL_INTEG_GAUSS51;
      default:
        return GSL_INTEG_GAUSS61;
    }
  }();

  gsl_error_handler_t* previous_handler = gsl_set_error_handler(gsl_error_handler_printout);
  const auto gslfunc = gsl_function{.function = func_integrand, .params = &params};
  const auto status =
      gsl_integration_qag(&gslfunc, a, b, 0., epsrel, GSLWSIZE, key, gslworkspace.get(), &result, abserr);
  gsl_set_error_handler(previous_handler);
  if (status != 0 && (status != 18 || (*abserr / std::abs(result)) > 0.1)) {
    printlnlog("[warning] integrator status {}. Integral value {:9.3e} +/- {:9.3e}", status, result, *abserr);
  }

#else

  // Boost's Gauss-Kronrod integrator
  result = boost::math::quadrature::gauss_kronrod<double, GKNPOINTS>::integrate(
      [&](double x) { return func_integrand(x, &params); }, a, b, 15, epsrel, abserr);

#endif
#endif
  return result;
}

#endif  // RATECOEFF_H
