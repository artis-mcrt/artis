// The numerical integrator used for the rate coefficient integrals: Boost's Gauss-Kronrod
// quadrature by default, or Simpson's rule when USE_SIMPSON_INTEGRATOR is defined.

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#ifdef USE_SIMPSON_INTEGRATOR
#include <algorithm>

#include "globals.h"
#include "mpi_logging.h"

template <class F>
constexpr auto simpson_integrator(const F& func_integrand, const double a, const double b, const int samplecount)
    -> double {
  assert_testmodeonly(samplecount % 2 == 1);

  const double deltax = (b - a) / (samplecount - 1);

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
#else
#include <cstdlib>

#pragma clang unsafe_buffer_usage begin
#include <boost/math/quadrature/gauss_kronrod.hpp>
#pragma clang unsafe_buffer_usage end

#endif

template <int GKNPOINTS = 61, class F>
auto integrator(const F& func_integrand, const double a, const double b, [[maybe_unused]] const double epsrel,
                double* abserr) -> double {
  static_assert(GKNPOINTS == 15 || GKNPOINTS == 31 || GKNPOINTS == 41 || GKNPOINTS == 51 || GKNPOINTS == 61,
                "Unsupported GKNPOINTS value");
#ifdef USE_SIMPSON_INTEGRATOR
  // need an odd number for Simpson rule
  const int samplecount = std::max(3, (globals::NPHIXSPOINTS * 2) + 1);

  *abserr = 0.;
  return simpson_integrator(func_integrand, a, b, samplecount);

#else

  // Boost's Gauss-Kronrod integrator
  return boost::math::quadrature::gauss_kronrod<double, GKNPOINTS>::integrate(func_integrand, a, b, 15, epsrel, abserr);

#endif
}

#endif  // INTEGRATOR_H
