#ifndef CHARGETRANSFER_H
#define CHARGETRANSFER_H

#include <span>

namespace chargetransfer {

void init();

// total charge transfer recombination rate [1/s] per ion of (element, upperion), for the transition
// upperion -> (upperion - 1). Includes the partner densities and the clumping factor. The NLTE
// matrix of the element under assembly covers the ions [first_ion_used, first_ion_used + nions_used).
// Call these two functions only with ENABLE_CHARGE_TRANSFER_REACTIONS, after init().
[[nodiscard]] auto ct_recombination_rate(int nonemptymgi, int element, int upperion, int first_ion_used, int nions_used)
    -> double;

// total charge transfer ionisation rate [1/s] per ion of (element, ion), for the transition
// ion -> (ion + 1). Includes the partner densities and the clumping factor. The NLTE matrix of the
// element under assembly covers the ions [first_ion_used, first_ion_used + nions_used).
[[nodiscard]] auto ct_ionisation_rate(int nonemptymgi, int element, int ion, int first_ion_used, int nions_used)
    -> double;

// the flat estimate [cm3/s] for an exothermic capture by a singly charged ion from a neutral donor,
// from the energy release deltae_erg. The unit tests call this directly.
[[nodiscard]] auto singly_charged_ratecoeff(double deltae_erg) -> double;

// evaluate the Kingdon & Ferland (1996) fit form. The unit tests call this directly.
// k = a * (T/1e4 K)^b * (1 + c * exp(d * T/1e4 K)) * exp(-eexp/T) [cm3/s], with a in [cm3/s].
// The fit part uses T clamped into [tmin, tmax]. The Boltzmann factor exp(-eexp/T) uses the given T.
[[nodiscard]] auto evaluate_ctfit(double a, double b, double c, double d, double eexp, double tmin, double tmax,
                                  double T) -> double;

// Landau-Zener cross section [cm2] for electron capture by an ion of the given charge from a neutral donor
// (method of Butler & Dalgarno 1980, ApJ, 241, 838), for the capture channels with the energy releases
// in deltae_erg_list. The trajectory crosses the channels in sequence, so the total transfer probability
// stays below one. ip_donor_erg is the ionisation energy of the neutral donor, which sets the exponent of
// the coupling. The charge must be two or more, and every energy release must be positive.
// The unit tests call this directly.
[[nodiscard]] auto sigma_lz_channels(int ioncharge, std::span<const double> deltae_erg_list, double ip_donor_erg,
                                     double v_cms) -> double;

// the cross section of a single capture channel, see sigma_lz_channels()
[[nodiscard]] auto sigma_lz_channel(int ioncharge, double deltae_erg, double ip_donor_erg, double v_cms) -> double;

}  // namespace chargetransfer

#endif  // CHARGETRANSFER_H
