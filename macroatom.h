#ifndef MACROATOM_H
#define MACROATOM_H

#include <cmath>

#include "constants.h"
#include "packet.h"

void macroatom_open_file(int my_rank);
void macroatom_close_file();

DEVICE_FUNC void do_macroatom(Packet& pkt, const MacroAtomState& pktmastate);

[[gnu::pure]] [[nodiscard]] auto rad_excitation_ratecoeff(int nonemptymgi, double upper_statweight, double einstein_A,
                                                          double epsilon_trans, double nnlevel_lower,
                                                          double nnlevel_upper, double statweight_lower,
                                                          int alltransindex, double t_current) -> double;

[[gnu::pure]] [[nodiscard]] auto rad_recombination_ratecoeff(float T_e, float nne, int element, int upperion,
                                                             int upperionlevel, int lowerionlevel, int nonemptymgi)
    -> double;

[[gnu::pure]] [[nodiscard]] auto stim_recombination_ratecoeff(float nne, int element, int upperion, int upper,
                                                              int lower, int nonemptymgi) -> double;

[[gnu::pure]] [[nodiscard]] auto col_recombination_ratecoeff(float T_e, float nne, int element, int upperion, int upper,
                                                             int lower, double epsilon_trans) -> double;

[[gnu::pure]] [[nodiscard]] auto col_ionisation_ratecoeff(float T_e, float nne, int element, int ion, int lower,
                                                          int phixstargetindex, double epsilon_trans) -> double;

[[gnu::pure]] [[nodiscard]] auto col_deexcitation_ratecoeff(float T_e, float nne, double epsilon_trans,
                                                            double upperstatweight, double lowerstatweight,
                                                            int alltransindex) -> double;

[[gnu::pure]] [[nodiscard]] auto col_excitation_ratecoeff(float T_e, float nne, double upperstatweight,
                                                          int alltransindex, double epsilon_trans,
                                                          double lowerstatweight) -> double;

// radiative deexcitation rate: paperII 3.5.2
// multiply by upper level population to get a rate per second
[[gnu::const]] [[nodiscard]] constexpr auto rad_deexcitation_ratecoeff(
    const double epsilon_trans, const float A_ul, const double upperstatweight, const double lowerstatweight,
    const double nnlevelupper, const double nnlevellower, const double t_current) -> double {
  const double nu_trans = epsilon_trans / H;

  const double B_ul = CLIGHTSQUAREDOVERTWOH / pow3(nu_trans) * A_ul;
  const double B_lu = upperstatweight / lowerstatweight * B_ul;

  const double tau_sobolev = ((B_lu * nnlevellower) - (B_ul * nnlevelupper)) * HCLIGHTOVERFOURPI * t_current;

  if (tau_sobolev > 1e-100) {
    const double beta = 1.0 / tau_sobolev * (-std::expm1(-tau_sobolev));
    const auto R = A_ul * beta;
    return R;
  }
  return 0.;
}

#endif  // MACROATOM_H
