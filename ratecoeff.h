#ifndef RATECOEFF_H
#define RATECOEFF_H

#include "constants.h"
#include "mpi_logging.h"

#ifdef USE_SIMPSON_INTEGRATOR
#include <algorithm>

#include "globals.h"

#else
#include <cstdlib>

#pragma clang unsafe_buffer_usage begin
#include <boost/math/quadrature/gauss_kronrod.hpp>
#pragma clang unsafe_buffer_usage end

#endif

void ratecoefficients_init();

void setup_photoion_luts();

[[nodiscard]] DEVICE_FUNC auto select_continuum_nu(int element, int lowerion, int lower, int upperionlevel, float T_e)
    -> double;
[[nodiscard]] DEVICE_FUNC auto get_spontrecombcoeff(int uniquelevelindex, int phixstargetindex, float T_e) -> double;

[[nodiscard]] DEVICE_FUNC auto get_bfcoolingcoeff(int element, int lowerion, int lowerionlevel, int phixstargetindex,
                                                  float T_e) -> double;

[[nodiscard]] DEVICE_FUNC auto get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex,
                                                     int nonemptymgi, bool use_cellcache) -> double;
[[nodiscard]] auto get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, int nonemptymgi)
    -> double;

[[nodiscard]] auto iongamma_is_zero(int nonemptymgi, int element, int ion) -> bool;

[[nodiscard]] auto calculate_iongamma_per_gspop(int nonemptymgi, int element, int ion) -> double;
[[nodiscard]] auto calculate_iongamma_per_ionpop(int nonemptymgi, int element, int lowerion,
                                                 bool collisional_not_radiative, bool force_bfintegral) -> double;

[[nodiscard]] auto calculate_ionrecombcoeff(int nonemptymgi, float T_e, int element, int upperion, bool assume_lte,
                                            bool collisional_not_radiative, bool lower_superlevel_only,
                                            bool per_groundmultipletpop) -> double;

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_ion_spontrecombcoeff(int uniqueionindex, float T_e) -> double;

template <class F>
constexpr auto simpson_integrator(F func_integrand, const double a, const double b, const int samplecount) -> double {
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

    integral += weight * func_integrand(x) * deltax;
  }
  integral /= 3.;

  return integral;
}

template <int GKNPOINTS = 61, class F>
auto integrator(F func_integrand, const double a, const double b, [[maybe_unused]] const double epsrel, double* abserr)
    -> double {
  static_assert(GKNPOINTS == 15 || GKNPOINTS == 31 || GKNPOINTS == 41 || GKNPOINTS == 51 || GKNPOINTS == 61,
                "Unsupported GKNPOINTS value");
  double result{0.};
#ifdef USE_SIMPSON_INTEGRATOR
  // need an odd number for Simpson rule
  const int samplecount = std::max(3, (globals::NPHIXSPOINTS * 3) + 1);

  result = simpson_integrator(func_integrand, a, b, samplecount);
  *abserr = 0.;

#else

  // Boost's Gauss-Kronrod integrator
  result =
      boost::math::quadrature::gauss_kronrod<double, GKNPOINTS>::integrate(func_integrand, a, b, 15, epsrel, abserr);

#endif
  return result;
}

#endif  // RATECOEFF_H
