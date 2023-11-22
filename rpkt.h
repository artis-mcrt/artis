#ifndef RPKT_H
#define RPKT_H

#include <ctime>

#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "sn3d.h"

void do_rpkt(struct packet *pkt_ptr, double t2);
void emit_rpkt(struct packet *pkt_ptr);
[[nodiscard]] auto closest_transition(double nu_cmf, int next_trans) -> int;
[[nodiscard]] auto calculate_chi_bf_gammacontr(int modelgridindex, double nu) -> double;
void calculate_chi_rpkt_cont(double nu_cmf, struct rpkt_continuum_absorptioncoeffs &chi_rpkt_cont, int modelgridindex,
                             bool usecellhistupdatephixslist);

[[nodiscard]] constexpr auto get_linedistance(const double prop_time, const double nu_cmf, const double nu_trans,
                                              const double d_nu_on_d_l) -> double {
  // distance from packet position to redshifting into line at frequency nu_trans

  if (nu_cmf <= nu_trans) {
    return 0;  /// photon was propagated too far, make sure that we don't miss a line
  }

  if constexpr (USE_RELATIVISTIC_DOPPLER_SHIFT) {
    // With special relativity, the Doppler shift formula has an extra factor of 1/gamma in it,
    // which changes the distance reach a line resonance and creates a dependence
    // on packet position and direction

    // use linear interpolation of frequency along the path
    return (nu_trans - nu_cmf) / d_nu_on_d_l;
  }

  return CLIGHT * prop_time * (nu_cmf / nu_trans - 1);
}

// wavelength bins are ordered by ascending wavelength (descending frequency)

constexpr auto get_wavelengthbin_nu_upper(const size_t binindex) -> double {
  const auto lambda_lower = expopac_lambdamin + binindex * expopac_deltalambda;
  return 1e8 * CLIGHT / lambda_lower;
}

constexpr auto get_wavelengthbin_nu_lower(const size_t binindex) -> double {
  const auto lambda_upper = expopac_lambdamin + (binindex + 1) * expopac_deltalambda;
  return 1e8 * CLIGHT / lambda_upper;
}

constexpr void calculate_binned_opacities(auto &expansionopacities, const int modelgridindex) {
  const time_t sys_time_start_calc_kpkt_rates = time(nullptr);

  printout("calculating binned expansion opacities for cell %d...", modelgridindex);

  const auto t_mid = globals::timesteps[globals::timestep].mid;

  // find the first line with nu below the upper limit of the first bin
  const auto *matchline =
      std::lower_bound(&globals::linelist[0], &globals::linelist[globals::nlines], get_wavelengthbin_nu_upper(0),
                       [](const auto &line, const double nu_cmf) -> bool { return line.nu > nu_cmf; });
  int lineindex = std::distance(globals::linelist, matchline);

  for (size_t wlbin = 0; wlbin < expopac_nbins; wlbin++) {
    double bin_linesum = 0.;

    const auto bin_nu_lower = get_wavelengthbin_nu_lower(wlbin);
    while (lineindex < globals::nlines && globals::linelist[lineindex].nu >= bin_nu_lower) {
      const float tau_line = get_tau_sobolev(modelgridindex, lineindex, t_mid);
      const auto linelambda = 1e8 * CLIGHT / globals::linelist[lineindex].nu;
      bin_linesum += (linelambda / expopac_deltalambda) * -std::expm1(-tau_line);
      lineindex++;
    }

    const float bin_kappa = 1. / (CLIGHT * t_mid * grid::get_rho(modelgridindex)) * bin_linesum;
    assert_always(std::isfinite(bin_kappa));
    expansionopacities[wlbin] = bin_kappa;
    // printout("bin %d: lambda %g to %g kappa %g kappa_grey %g\n", wlbin, get_wavelengthbin_lambda_lower(wlbin),
    //          get_wavelengthbin_lambda_upper(wlbin), bin_kappa, grid::modelgrid[modelgridindex].kappagrey);
  }
  printout("took %ld seconds\n", time(nullptr) - sys_time_start_calc_kpkt_rates);
}
#endif  // RPKT_H
