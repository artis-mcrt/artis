// Declarations for gamma-ray packet transport (gammapkt.cc).

#ifndef GAMMAPKT_H
#define GAMMAPKT_H

#include <array>
#include <cmath>

#include "constants.h"
#include "mpi_logging.h"
#include "packet.h"
#include "random.h"

namespace gammapkt {
constexpr double pair_rest_mass_energy = 1.022 * MEV;
constexpr double pair_production_threshold_nu = pair_rest_mass_energy / H;

// Below the threshold, pair production transfers no energy to kinetic energy.
[[nodiscard]] constexpr auto get_pair_annihilation_probability(const double nu_cmf) -> double {
  return nu_cmf <= pair_production_threshold_nu ? 1. : pair_production_threshold_nu / nu_cmf;
}

// energy-dependent factor of the pair-production cross section, from equation 2 of Ambwani & Sutherland (1988),
// ApJ, 325, 820-827, doi:10.1086/166052, after Hubbell (1969), NSRDS-NBS 29, National Bureau of Standards,
// doi:10.6028/NBS.NSRDS.29. Multiply by Z^2 * 1e-27 for a cross section in cm^2. Only valid above the 1.022 MeV
// threshold.
[[nodiscard]] constexpr auto get_sigma_pair_prod_factor(const double nu_cmf) -> double {
  if (nu_cmf <= pair_production_threshold_nu) {
    return 0.;
  }
  const double hnu_over_1MeV = nu_cmf * H / MEV;
  if (nu_cmf > (1.5 * MEV / H)) {
    return 0.0481 + (0.301 * (hnu_over_1MeV - 1.5));
  }
  // the coefficient below 1.5 MeV is 0.10063, not the 1.0063 printed in the paper: the latter makes
  // the fit jump by a factor of ten at the junction between the two branches
  return 0.10063 * (hnu_over_1MeV - 1.022);
}

// The branches agree at 1.5 MeV to the precision of the fit coefficients.
static_assert((get_sigma_pair_prod_factor((1.5 * MEV / H) * (1. + 1e-12)) -
               get_sigma_pair_prod_factor((1.5 * MEV / H))) < 2e-6);
static_assert((get_sigma_pair_prod_factor((1.5 * MEV / H)) -
               get_sigma_pair_prod_factor((1.5 * MEV / H) * (1. + 1e-12))) < 2e-6);

void init_gamma_data();
DEVICE_FUNC void pellet_gamma_decay(Packet& pkt);
DEVICE_FUNC void do_gamma(Packet& pkt, int nts, double t2);

// Frequency [Hz] of a random gamma ray from a nuclide's decay spectrum. A nuclide without a gamma-ray table has a
// spectrum of one line at its mean gamma energy per decay.
auto choose_gamma_ray(int nucindex, rngstate_type& rngstate) -> double;

// The Compton cross section of Klein & Nishina (1929), Z. Phys., 52, 853-868, doi:10.1007/BF01366453, integrated
// over the energy loss factor from 1 up to f_max.
// Passing f_max = 1 + 2x (the maximum possible loss, i.e. backscattering) gives the total cross section.
// - x: the photon energy in units of the electron rest mass
// - f_max: the energy loss factor (nu_before / nu_after) up to which we integrate
[[nodiscard]] constexpr auto sigma_compton_partial(const double x, const double f_max) -> double {
  const double term1 = ((x * x) - (2 * x) - 2) * std::log(f_max) / x / x;
  const double term2 = (((f_max * f_max) - 1) / (f_max * f_max)) / 2;
  const double term3 = ((f_max - 1) / x) * ((1 / x) + (2 / f_max) + (1 / (x * f_max)));

  return (3 * SIGMA_T * (term1 + term2 + term3) / (8 * x));
}

// Sample an energy loss factor f from the Klein-Nishina distribution, by bisecting for the f at which the
// partial cross section reaches the fraction zrand of the total.
[[nodiscard]] inline auto choose_f(const double xx, const double zrand) -> double {
  double f_max = 1 + (2 * xx);
  double f_min = 1;

  const double norm = zrand * sigma_compton_partial(xx, f_max);

  int count = 0;
  double err = 1e20;

  double ftry = (f_max + f_min) / 2;
  while ((err > 1.e-4) && (count < 1000)) {
    ftry = (f_max + f_min) / 2;
    const double sigma_try = sigma_compton_partial(xx, ftry);
    if (sigma_try > norm) {
      f_max = ftry;
      err = (sigma_try - norm) / norm;
    } else {
      f_min = ftry;
      err = (norm - sigma_try) / norm;
    }

    count++;
  }

  assert_testmodeonly(ftry >= 1.);
  assert_testmodeonly(ftry <= ((2 * xx) + 1.));
  return ftry;
}

// Compute the mean energy converted to non-thermal electrons times the Klein-Nishina cross section.
constexpr auto meanf_sigma(const double x) -> double {
  if (x < THOMSON_LIMIT) {
    // The closed form below sums five terms of order 1/x^2 that have to cancel down to O(x), so it
    // loses about 8 decimal digits per decade in x: the relative error is 2e-9 at x = 1e-2, 2e-6 at
    // 1e-3, 25% at 1e-4, and the result turns negative below ~1e-6. Use the exact Taylor series of
    // the same expression instead, which is accurate to 6e-11 at the x = THOMSON_LIMIT crossover
    // and better below it. The leading term is the Thomson-limit mean fractional energy transfer x.
    constexpr std::array taylor_coeffs{
        1., -21. / 5., 147. / 10., -1616. / 35., 940. / 7., -2584. / 7., 14588. / 15., -409088. / 165.,
    };
    // Horner evaluation, starting from the highest-order coefficient
    double series = taylor_coeffs.back();
    for (int i = static_cast<int>(taylor_coeffs.size()) - 2; i >= 0; i--) {
      series = taylor_coeffs[i] + (x * series);
    }
    return SIGMA_T * x * series;
  }

  const double f = 1 + (2 * x);

  const double term0 = 2 / x;
  const double term1 = (1 - (2 / x) - (3 / (x * x))) * log(f);
  const double term2 = ((4 / x) + (3 / (x * x)) - 1) * 2 * x / f;
  const double term3 = (1 - (2 / x) - (1 / (x * x))) * 2 * x * (1 + x) / f / f;
  const double term4 = -2. * x * ((4 * x * x) + (6 * x) + 3) / 3 / f / f / f;

  const double tot = 3 * SIGMA_T * (term0 + term1 + term2 + term3 + term4) / (8 * x);

  return tot;
}

}  // namespace gammapkt

#endif  // GAMMAPKT_H
