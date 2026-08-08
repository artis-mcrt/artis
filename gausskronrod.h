// Adaptive Gauss-Kronrod quadrature extracted from Boost.Math (boost/math/quadrature/gauss_kronrod.hpp and
// gauss.hpp) so that ARTIS does not depend on the full library. The abscissae/weight tables and the floating-point
// operation sequence are copied exactly, so results are bit-identical to those from
// boost::math::quadrature::gauss_kronrod<double, NPOINTS>::integrate(). Tables for other NPOINTS values (15, 21, 41,
// 51) can be copied from Boost.Math if they are ever needed.
//
// (C) Copyright John Maddock 2017.
// (C) Copyright Nick Thompson 2017.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef GAUSSKRONROD_H
#define GAUSSKRONROD_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace gausskronrod_detail {

// For each NPOINTS-point Kronrod rule: the non-negative Kronrod abscissae (zero first, then interleaved Gauss and
// Kronrod-only points in increasing order) with their Kronrod weights, plus the weights of the embedded
// (NPOINTS - 1) / 2 point Gauss rule. The values are the same 20-significant-digit decimal literals used by
// Boost.Math (which stores them cast from long double, reproduced here so the rounding matches on every platform).
template <int NPOINTS>
struct GaussKronrodTables;

// The long double literals below are Boost.Math's own storage format for the double-precision tables; keeping the
// long-double-to-double conversion reproduces Boost's exact values on every platform.
// NOLINTBEGIN(google-runtime-float)
template <>
struct GaussKronrodTables<31> {
  static constexpr std::array<double, 16> abscissa{
      static_cast<double>(0.00000000000000000000000000000000000e+00L),
      static_cast<double>(1.01142066918717499027074231447392339e-01L),
      static_cast<double>(2.01194093997434522300628303394596208e-01L),
      static_cast<double>(2.99180007153168812166780024266388963e-01L),
      static_cast<double>(3.94151347077563369897207370981045468e-01L),
      static_cast<double>(4.85081863640239680693655740232350613e-01L),
      static_cast<double>(5.70972172608538847537226737253910641e-01L),
      static_cast<double>(6.50996741297416970533735895313274693e-01L),
      static_cast<double>(7.24417731360170047416186054613938010e-01L),
      static_cast<double>(7.90418501442465932967649294817947347e-01L),
      static_cast<double>(8.48206583410427216200648320774216851e-01L),
      static_cast<double>(8.97264532344081900882509656454495883e-01L),
      static_cast<double>(9.37273392400705904307758947710209471e-01L),
      static_cast<double>(9.67739075679139134257347978784337225e-01L),
      static_cast<double>(9.87992518020485428489565718586612581e-01L),
      static_cast<double>(9.98002298693397060285172840152271209e-01L),
  };
  static constexpr std::array<double, 16> weights{
      static_cast<double>(1.01330007014791549017374792767492547e-01L),
      static_cast<double>(1.00769845523875595044946662617569722e-01L),
      static_cast<double>(9.91735987217919593323931734846031311e-02L),
      static_cast<double>(9.66427269836236785051799076275893351e-02L),
      static_cast<double>(9.31265981708253212254868727473457186e-02L),
      static_cast<double>(8.85644430562117706472754436937743032e-02L),
      static_cast<double>(8.30805028231330210382892472861037896e-02L),
      static_cast<double>(7.68496807577203788944327774826590067e-02L),
      static_cast<double>(6.98541213187282587095200770991474758e-02L),
      static_cast<double>(6.20095678006706402851392309608029322e-02L),
      static_cast<double>(5.34815246909280872653431472394302968e-02L),
      static_cast<double>(4.45897513247648766082272993732796902e-02L),
      static_cast<double>(3.53463607913758462220379484783600481e-02L),
      static_cast<double>(2.54608473267153201868740010196533594e-02L),
      static_cast<double>(1.50079473293161225383747630758072681e-02L),
      static_cast<double>(5.37747987292334898779205143012764982e-03L),
  };
  static constexpr std::array<double, 8> gauss_weights{
      static_cast<double>(2.02578241925561272880620199967519315e-01L),
      static_cast<double>(1.98431485327111576456118326443839325e-01L),
      static_cast<double>(1.86161000015562211026800561866422825e-01L),
      static_cast<double>(1.66269205816993933553200860481208811e-01L),
      static_cast<double>(1.39570677926154314447804794511028323e-01L),
      static_cast<double>(1.07159220467171935011869546685869303e-01L),
      static_cast<double>(7.03660474881081247092674164506673385e-02L),
      static_cast<double>(3.07532419961172683546283935772044177e-02L),
  };
};

template <>
struct GaussKronrodTables<61> {
  static constexpr std::array<double, 31> abscissa{
      static_cast<double>(0.00000000000000000000000000000000000e+00L),
      static_cast<double>(5.14718425553176958330252131667225737e-02L),
      static_cast<double>(1.02806937966737030147096751318000592e-01L),
      static_cast<double>(1.53869913608583546963794672743255920e-01L),
      static_cast<double>(2.04525116682309891438957671002024710e-01L),
      static_cast<double>(2.54636926167889846439805129817805108e-01L),
      static_cast<double>(3.04073202273625077372677107199256554e-01L),
      static_cast<double>(3.52704725530878113471037207089373861e-01L),
      static_cast<double>(4.00401254830394392535476211542660634e-01L),
      static_cast<double>(4.47033769538089176780609900322854000e-01L),
      static_cast<double>(4.92480467861778574993693061207708796e-01L),
      static_cast<double>(5.36624148142019899264169793311072794e-01L),
      static_cast<double>(5.79345235826361691756024932172540496e-01L),
      static_cast<double>(6.20526182989242861140477556431189299e-01L),
      static_cast<double>(6.60061064126626961370053668149270753e-01L),
      static_cast<double>(6.97850494793315796932292388026640068e-01L),
      static_cast<double>(7.33790062453226804726171131369527646e-01L),
      static_cast<double>(7.67777432104826194917977340974503132e-01L),
      static_cast<double>(7.99727835821839083013668942322683241e-01L),
      static_cast<double>(8.29565762382768397442898119732501916e-01L),
      static_cast<double>(8.57205233546061098958658510658943857e-01L),
      static_cast<double>(8.82560535792052681543116462530225590e-01L),
      static_cast<double>(9.05573307699907798546522558925958320e-01L),
      static_cast<double>(9.26200047429274325879324277080474004e-01L),
      static_cast<double>(9.44374444748559979415831324037439122e-01L),
      static_cast<double>(9.60021864968307512216871025581797663e-01L),
      static_cast<double>(9.73116322501126268374693868423706885e-01L),
      static_cast<double>(9.83668123279747209970032581605662802e-01L),
      static_cast<double>(9.91630996870404594858628366109485725e-01L),
      static_cast<double>(9.96893484074649540271630050918695283e-01L),
      static_cast<double>(9.99484410050490637571325895705810819e-01L),
  };
  static constexpr std::array<double, 31> weights{
      static_cast<double>(5.14947294294515675583404336470993075e-02L),
      static_cast<double>(5.14261285374590259338628792157812598e-02L),
      static_cast<double>(5.12215478492587721706562826049442083e-02L),
      static_cast<double>(5.08817958987496064922974730498046919e-02L),
      static_cast<double>(5.04059214027823468408930856535850289e-02L),
      static_cast<double>(4.97956834270742063578115693799423285e-02L),
      static_cast<double>(4.90554345550297788875281653672381736e-02L),
      static_cast<double>(4.81858617570871291407794922983045926e-02L),
      static_cast<double>(4.71855465692991539452614781810994865e-02L),
      static_cast<double>(4.60592382710069881162717355593735806e-02L),
      static_cast<double>(4.48148001331626631923555516167232438e-02L),
      static_cast<double>(4.34525397013560693168317281170732581e-02L),
      static_cast<double>(4.19698102151642461471475412859697578e-02L),
      static_cast<double>(4.03745389515359591119952797524681142e-02L),
      static_cast<double>(3.86789456247275929503486515322810503e-02L),
      static_cast<double>(3.68823646518212292239110656171359677e-02L),
      static_cast<double>(3.49793380280600241374996707314678751e-02L),
      static_cast<double>(3.29814470574837260318141910168539275e-02L),
      static_cast<double>(3.09072575623877624728842529430922726e-02L),
      static_cast<double>(2.87540487650412928439787853543342111e-02L),
      static_cast<double>(2.65099548823331016106017093350754144e-02L),
      static_cast<double>(2.41911620780806013656863707252320268e-02L),
      static_cast<double>(2.18280358216091922971674857383389934e-02L),
      static_cast<double>(1.94141411939423811734089510501284559e-02L),
      static_cast<double>(1.69208891890532726275722894203220924e-02L),
      static_cast<double>(1.43697295070458048124514324435800102e-02L),
      static_cast<double>(1.18230152534963417422328988532505929e-02L),
      static_cast<double>(9.27327965951776342844114689202436042e-03L),
      static_cast<double>(6.63070391593129217331982636975016813e-03L),
      static_cast<double>(3.89046112709988405126720184451550328e-03L),
      static_cast<double>(1.38901369867700762455159122675969968e-03L),
  };
  static constexpr std::array<double, 15> gauss_weights{
      static_cast<double>(1.02852652893558840341285636705415044e-01L),
      static_cast<double>(1.01762389748405504596428952168554045e-01L),
      static_cast<double>(9.95934205867952670627802821035694765e-02L),
      static_cast<double>(9.63687371746442596394686263518098651e-02L),
      static_cast<double>(9.21225222377861287176327070876187672e-02L),
      static_cast<double>(8.68997872010829798023875307151257026e-02L),
      static_cast<double>(8.07558952294202153546949384605297309e-02L),
      static_cast<double>(7.37559747377052062682438500221907342e-02L),
      static_cast<double>(6.59742298821804951281285151159623612e-02L),
      static_cast<double>(5.74931562176190664817216894020561288e-02L),
      static_cast<double>(4.84026728305940529029381404228075178e-02L),
      static_cast<double>(3.87991925696270495968019364463476920e-02L),
      static_cast<double>(2.87847078833233693497191796112920436e-02L),
      static_cast<double>(1.84664683110909591423021319120472691e-02L),
      static_cast<double>(7.96819249616660561546588347467362245e-03L),
  };
};
// NOLINTEND(google-runtime-float)

// Evaluate the Kronrod sum over [-1, 1] and, from the difference to the embedded Gauss sum, the error estimate.
template <int NPOINTS, class F>
auto integrate_non_adaptive_m1_1(F f, double* error) -> double {
  using Tables = GaussKronrodTables<NPOINTS>;
  constexpr unsigned gauss_order = (NPOINTS - 1) / 2;

  unsigned gauss_start = 2;
  unsigned kronrod_start = 1;
  const double f_centre = f(0.);
  double kronrod_result = f_centre * Tables::weights[0];
  double gauss_result = 0.;
  if constexpr (gauss_order & 1U) {
    gauss_result += f_centre * Tables::gauss_weights[0];
  } else {
    gauss_start = 1;
    kronrod_start = 2;
  }
  for (unsigned i = gauss_start; i < Tables::abscissa.size(); i += 2) {
    const double fp = f(Tables::abscissa[i]);
    const double fm = f(-Tables::abscissa[i]);
    kronrod_result += (fp + fm) * Tables::weights[i];
    gauss_result += (fp + fm) * Tables::gauss_weights[i / 2];
  }
  for (unsigned i = kronrod_start; i < Tables::abscissa.size(); i += 2) {
    const double fp = f(Tables::abscissa[i]);
    const double fm = f(-Tables::abscissa[i]);
    kronrod_result += (fp + fm) * Tables::weights[i];
  }
  if (error != nullptr) {
    *error = std::max(std::fabs(kronrod_result - gauss_result),
                      std::fabs(kronrod_result * std::numeric_limits<double>::epsilon() * 2));
  }
  return kronrod_result;
}

// Apply the rule to [a, b] and recursively bisect while the error estimate exceeds both the local and the inherited
// absolute tolerance (each derived from the relative tolerance tol), up to max_levels of bisection.
template <int NPOINTS, class F>
auto recursive_adaptive_integrate(const F& f, const double tol, const double a, const double b,
                                  const unsigned max_levels, double abs_tol, double* error) -> double {
  double error_local = 0.;
  const double mean = (b + a) / 2;
  const double scale = (b - a) / 2;
  const auto ff = [&](const double x) -> double { return f((scale * x) + mean); };
  const double r1 = integrate_non_adaptive_m1_1<NPOINTS>(ff, &error_local);
  double estimate = scale * r1;

  const double abs_tol1 = std::fabs(estimate * tol);
  if (abs_tol == 0) {
    abs_tol = abs_tol1;
  }

  if ((max_levels != 0) && (abs_tol1 < error_local) && (abs_tol < error_local)) {
    const double mid = (a + b) / 2;
    estimate = recursive_adaptive_integrate<NPOINTS>(f, tol, a, mid, max_levels - 1, abs_tol / 2, error);
    estimate += recursive_adaptive_integrate<NPOINTS>(f, tol, mid, b, max_levels - 1, abs_tol / 2, &error_local);
    if (error != nullptr) {
      *error += error_local;
    }
    return estimate;
  }
  if (error != nullptr) {
    *error = error_local;
  }
  return estimate;
}

}  // namespace gausskronrod_detail

// Integrate f over the finite interval [a, b] to relative tolerance tol with an adaptive NPOINTS-point Gauss-Kronrod
// rule, subdividing at most max_depth times. If error is non-null it receives an estimate of the absolute error.
template <int NPOINTS, class F>
auto gauss_kronrod_integrate(F f, const double a, const double b, const unsigned max_depth, const double tol,
                             double* error) -> double {
  static_assert(NPOINTS == 31 || NPOINTS == 61, "no Gauss-Kronrod table for this NPOINTS value");
  if (!std::isfinite(a) || !std::isfinite(b)) {
    throw std::domain_error("gauss_kronrod_integrate: integration bounds must be finite");
  }
  if (a == b) {
    return 0.;
  }
  if (b < a) {
    return -gausskronrod_detail::recursive_adaptive_integrate<NPOINTS>(f, tol, b, a, max_depth, 0., error);
  }
  return gausskronrod_detail::recursive_adaptive_integrate<NPOINTS>(f, tol, a, b, max_depth, 0., error);
}

#endif  // GAUSSKRONROD_H
