// TOMS Algorithm 748 bracketing root finder (Alefeld, Potra & Shi 1995, ACM Trans. Math. Softw. 21, 327,
// doi:10.1145/210089.210111), extracted from Boost.Math (boost/math/tools/toms748_solve.hpp) so that ARTIS does not
// depend on the full library. The floating-point operation sequence is copied exactly, so results are bit-identical
// to boost::math::tools::toms748_solve<double>.
//
// (C) Copyright John Maddock 2006.
// (C) Copyright Matt Borland 2024.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef TOMS748_H
#define TOMS748_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace toms748_detail {

[[nodiscard]] constexpr auto sign(const double z) -> int {
  if (z == 0) {
    return 0;
  }
  return std::signbit(z) ? -1 : 1;
}

// Given a point c inside the existing enclosing interval [a, b] sets a = c if f(c) == 0, otherwise finds the new
// enclosing interval: either [a, c] or [c, b] and sets d and fd to the point that has just been removed from the
// interval. In other words d is the third best guess to the root.
template <class F>
void bracket(F f, double& a, double& b, double c, double& fa, double& fb, double& d, double& fd) {
  const double tol = std::numeric_limits<double>::epsilon() * 2;

  // If the interval [a,b] is very small, or if c is too close to one end of the interval then we need to adjust the
  // location of c accordingly:
  if ((b - a) < 2 * tol * a) {
    c = a + ((b - a) / 2);
  } else if (c <= a + (std::fabs(a) * tol)) {
    c = a + (std::fabs(a) * tol);
  } else if (c >= b - (std::fabs(b) * tol)) {
    c = b - (std::fabs(b) * tol);
  }

  const double fc = f(c);

  // if we have a zero then we have an exact solution to the root:
  if (fc == 0) {
    a = c;
    fa = 0;
    d = 0;
    fd = 0;
    return;
  }

  // Non-zero fc, update the interval:
  if (sign(fa) * sign(fc) < 0) {
    d = b;
    fd = fb;
    b = c;
    fb = fc;
  } else {
    d = a;
    fd = fa;
    a = c;
    fa = fc;
  }
}

// return num / denom without overflow, return r if overflow would occur
[[nodiscard]] inline auto safe_div(const double num, const double denom, const double r) -> double {
  if (std::fabs(denom) < 1 && std::fabs(denom * std::numeric_limits<double>::max()) <= std::fabs(num)) {
    return r;
  }
  return num / denom;
}

// Performs standard secant interpolation of [a,b] given function evaluations f(a) and f(b). Performs a bisection if
// secant interpolation would leave us very close to either a or b. Rationale: we only call this function when at
// least one other form of interpolation has already failed, so we know that the function is unlikely to be smooth
// with a root very close to a or b.
[[nodiscard]] inline auto secant_interpolate(const double a, const double b, const double fa, const double fb)
    -> double {
  const double tol = std::numeric_limits<double>::epsilon() * 5;
  const double c = a - ((fa / (fb - fa)) * (b - a));
  if ((c <= a + (std::fabs(a) * tol)) || (c >= b - (std::fabs(b) * tol))) {
    return (a + b) / 2;
  }
  return c;
}

// Performs quadratic interpolation to determine the next point, takes count Newton steps to find the location of the
// quadratic polynomial. Point d must lie outside of the interval [a,b], it is the third best approximation to the
// root, after a and b. Note: this does not guarantee to find a root inside [a, b], so we fall back to a secant step
// should the result be out of range.
[[nodiscard]] inline auto quadratic_interpolate(const double a, const double b, const double d, const double fa,
                                                const double fb, const double fd, const unsigned count) -> double {
  // Start by obtaining the coefficients of the quadratic polynomial:
  const double B = safe_div(fb - fa, b - a, std::numeric_limits<double>::max());
  double A = safe_div(fd - fb, d - b, std::numeric_limits<double>::max());
  A = safe_div(A - B, d - a, 0.);

  if (A == 0) {
    // failure to determine coefficients, try a secant step:
    return secant_interpolate(a, b, fa, fb);
  }

  // Determine the starting point of the Newton steps:
  double c = (sign(A) * sign(fa) > 0) ? a : b;

  // Take the Newton steps:
  for (unsigned i = 1; i <= count; ++i) {
    c -= safe_div(fa + ((B + (A * (c - b))) * (c - a)), B + (A * ((2 * c) - a - b)), 1 + c - a);
  }
  if ((c <= a) || (c >= b)) {
    // Oops, failure, try a secant step:
    c = secant_interpolate(a, b, fa, fb);
  }
  return c;
}

// Uses inverse cubic interpolation of f(x) at points [a,b,d,e] to obtain an approximate root of f(x). Points d and e
// lie outside the interval [a,b] and are the third and fourth best approximations to the root that we have found so
// far. Note: this does not guarantee to find a root inside [a, b], so we fall back to quadratic interpolation in
// case of an erroneous result.
[[nodiscard]] inline auto cubic_interpolate(const double a, const double b, const double d, const double e,
                                            const double fa, const double fb, const double fd, const double fe)
    -> double {
  const double q11 = (d - e) * fd / (fe - fd);
  const double q21 = (b - d) * fb / (fd - fb);
  const double q31 = (a - b) * fa / (fb - fa);
  const double d21 = (b - d) * fd / (fd - fb);
  const double d31 = (a - b) * fb / (fb - fa);

  const double q22 = (d21 - q11) * fb / (fe - fb);
  const double q32 = (d31 - q21) * fa / (fd - fa);
  const double d32 = (d31 - q21) * fd / (fd - fa);
  const double q33 = (d32 - q22) * fa / (fe - fa);
  double c = q31 + q32 + q33 + a;

  if ((c <= a) || (c >= b)) {
    // Out of bounds step, fall back to quadratic interpolation:
    c = quadratic_interpolate(a, b, d, fa, fb, fd, 3);
  }

  return c;
}

}  // namespace toms748_detail

// Main entry point and logic for Toms Algorithm 748 root finder. Given a function f and an interval [ax, bx] with
// fax = f(ax) and fbx = f(bx) of opposite signs, iteratively narrows the bracket until the tolerance functor
// tol(a, b) is satisfied or max_iter function evaluations have been used. On return max_iter is set to the number of
// function evaluations performed and the current bracket is returned.
template <class F, class Tol>
auto toms748_solve(F f, const double ax, const double bx, const double fax, const double fbx, Tol tol,
                   std::uintmax_t& max_iter) -> std::pair<double, double> {
  // Sanity check - are we allowed to iterate at all?
  if (max_iter == 0) {
    return {ax, bx};
  }

  std::uintmax_t count = max_iter;
  constexpr double mu = 0.5;

  // initialise a, b and fa, fb:
  double a = ax;
  double b = bx;
  if (a >= b) {
    // GPU builds cannot throw; return a NaN bracket like Boost.Math's BOOST_MATH_NO_EXCEPTIONS error handler
#ifdef GPU_ON
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
#else
    throw std::domain_error("toms748_solve: parameters a and b out of order");
#endif
  }
  double fa = fax;
  double fb = fbx;

  if (tol(a, b) || (fa == 0) || (fb == 0)) {
    max_iter = 0;
    if (fa == 0) {
      b = a;
    } else if (fb == 0) {
      a = b;
    }
    return {a, b};
  }

  if (toms748_detail::sign(fa) * toms748_detail::sign(fb) > 0) {
#ifdef GPU_ON
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
#else
    throw std::domain_error("toms748_solve: parameters a and b do not bracket the root");
#endif
  }

  // dummy value for fd, e and fe:
  double d = 0.;
  double fd = 1e5;
  double e = 1e5;
  double fe = 1e5;
  // cppcheck-suppress unreadVariable
  double c{0.};

  if (fa != 0) {
    // On the first step we take a secant step:
    c = toms748_detail::secant_interpolate(a, b, fa, fb);
    toms748_detail::bracket(f, a, b, c, fa, fb, d, fd);
    --count;

    if (count && (fa != 0) && !tol(a, b)) {
      // On the second step we take a quadratic interpolation:
      c = toms748_detail::quadratic_interpolate(a, b, d, fa, fb, fd, 2);
      e = d;
      fe = fd;
      toms748_detail::bracket(f, a, b, c, fa, fb, d, fd);
      --count;
    }
  }

  while (count && (fa != 0) && !tol(a, b)) {
    // save our brackets:
    const double a0 = a;
    const double b0 = b;

    // Starting with the third step taken we can use either quadratic or cubic interpolation. Cubic interpolation
    // requires that all four function values fa, fb, fd, and fe are distinct, should that not be the case then
    // variable prof will get set to true, and we'll end up taking a quadratic step instead.
    const double min_diff = std::numeric_limits<double>::min() * 32;
    bool prof = (std::fabs(fa - fb) < min_diff) || (std::fabs(fa - fd) < min_diff) || (std::fabs(fa - fe) < min_diff) ||
                (std::fabs(fb - fd) < min_diff) || (std::fabs(fb - fe) < min_diff) || (std::fabs(fd - fe) < min_diff);
    if (prof) {
      c = toms748_detail::quadratic_interpolate(a, b, d, fa, fb, fd, 2);
    } else {
      c = toms748_detail::cubic_interpolate(a, b, d, e, fa, fb, fd, fe);
    }

    // re-bracket, and check for termination:
    e = d;
    fe = fd;
    toms748_detail::bracket(f, a, b, c, fa, fb, d, fd);
    if ((0 == --count) || (fa == 0) || tol(a, b)) {
      break;
    }

    // Now another interpolated step:
    prof = (std::fabs(fa - fb) < min_diff) || (std::fabs(fa - fd) < min_diff) || (std::fabs(fa - fe) < min_diff) ||
           (std::fabs(fb - fd) < min_diff) || (std::fabs(fb - fe) < min_diff) || (std::fabs(fd - fe) < min_diff);
    if (prof) {
      c = toms748_detail::quadratic_interpolate(a, b, d, fa, fb, fd, 3);
    } else {
      c = toms748_detail::cubic_interpolate(a, b, d, e, fa, fb, fd, fe);
    }

    // Bracket again, and check termination condition, update e:
    toms748_detail::bracket(f, a, b, c, fa, fb, d, fd);
    if ((0 == --count) || (fa == 0) || tol(a, b)) {
      break;
    }

    // Now we take a double-length secant step:
    double u{0.};
    double fu{0.};
    if (std::fabs(fa) < std::fabs(fb)) {
      u = a;
      fu = fa;
    } else {
      u = b;
      fu = fb;
    }
    c = u - (2 * (fu / (fb - fa)) * (b - a));
    if (std::fabs(c - u) > (b - a) / 2) {
      c = a + ((b - a) / 2);
    }

    // Bracket again, and check termination condition:
    e = d;
    fe = fd;
    toms748_detail::bracket(f, a, b, c, fa, fb, d, fd);
    if ((0 == --count) || (fa == 0) || tol(a, b)) {
      break;
    }

    // And finally... check to see if an additional bisection step is to be taken, we do this if we're not converging
    // fast enough:
    if ((b - a) < mu * (b0 - a0)) {
      continue;
    }

    // bracket again on a bisection:
    e = d;
    fe = fd;
    toms748_detail::bracket(f, a, b, a + ((b - a) / 2), fa, fb, d, fd);
    --count;
  }  // while loop

  max_iter -= count;
  if (fa == 0) {
    b = a;
  } else if (fb == 0) {
    a = b;
  }
  return {a, b};
}

#endif  // TOMS748_H
