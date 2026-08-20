#ifndef CHARGETRANSFER_H
#define CHARGETRANSFER_H

namespace chargetransfer {

void init();

[[nodiscard]] auto get_reaction_count() -> int;

// total charge transfer recombination rate [1/s] per ion of (element, upperion),
// for the transition upperion -> (upperion - 1). Includes the partner densities and the clumping factor.
[[nodiscard]] auto ct_recombination_rate(int nonemptymgi, int element, int upperion) -> double;

// total charge transfer ionisation rate [1/s] per ion of (element, ion),
// for the transition ion -> (ion + 1). Includes the partner densities and the clumping factor.
[[nodiscard]] auto ct_ionisation_rate(int nonemptymgi, int element, int ion) -> double;

// evaluate the Kingdon & Ferland (1996) fit form. The unit tests call this directly.
// k = a * (T/1e4 K)^b * (1 + c * exp(d * T/1e4 K)) * exp(-eexp/T) [cm3/s], with a in [cm3/s].
// The fit part uses T clamped into [tmin, tmax]. The Boltzmann factor exp(-eexp/T) uses the given T.
[[nodiscard]] auto evaluate_ctfit(double a, double b, double c, double d, double eexp, double tmin, double tmax,
                                  double T) -> double;

// Landau-Zener cross section [cm2] for electron capture by an ion of the given charge from a neutral donor,
// for a single final channel with energy defect deltae_erg (method of Butler & Dalgarno 1980, ApJ, 241, 838).
// ip_donor_erg is the ionisation energy of the neutral donor, which sets the exponent of the coupling.
// The unit tests call this directly.
[[nodiscard]] auto sigma_lz_channel(int ioncharge, double deltae_erg, double ip_donor_erg, double v_cms) -> double;

}  // namespace chargetransfer

#endif  // CHARGETRANSFER_H
