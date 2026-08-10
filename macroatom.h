// Declarations for the macro-atom machinery (macroatom.cc).

#ifndef MACROATOM_H
#define MACROATOM_H

#include <cmath>

#include "constants.h"
#include "packet.h"

void macroatom_open_file();

// Follow an activated macro-atom through stochastic internal energy-flow transitions until it deactivates as an
// r-packet or k-packet. Actions are selected in proportion to the local radiative, collisional, and non-thermal rates.
// Lucy (2002), doi:10.1051/0004-6361:20011756; Lucy (2003), arXiv:astro-ph/0303202.
DEVICE_FUNC void do_macroatom(Packet& pkt, const MacroAtomState& pktmastate);

// prepopulate one unique level's macroatom transition rates into the cellcache. Used in GPU mode, where
// the lazy mutex-guarded calculation in do_macroatom() cannot run safely on the device.
DEVICE_FUNC void calculate_cellcache_macroatom_transitionrates(int nonemptymgi, int uniquelevelindex, double t_mid);

[[gnu::pure]] [[nodiscard]] auto rad_excitation_ratecoeff(int nonemptymgi, double upper_statweight, double einstein_A,
                                                          double epsilon_trans, double nnlevel_lower,
                                                          double nnlevel_upper, double statweight_lower,
                                                          int alltransindex, double t_current) -> double;

[[gnu::pure]] [[nodiscard]] auto rad_recombination_ratecoeff(float T_e, float clumpednne, int element, int upperion,
                                                             int upperionlevel, int lowerionlevel) -> double;

[[gnu::pure]] [[nodiscard]] auto col_recombination_ratecoeff(float T_e, float clumpednne, int element, int upperion,
                                                             int upper, int lower, double epsilon_trans) -> double;

[[gnu::pure]] [[nodiscard]] auto col_ionisation_ratecoeff(float T_e, float clumpednne, int element, int ion, int lower,
                                                          int phixstargetindex, double epsilon_trans) -> double;

[[gnu::pure]] [[nodiscard]] auto col_deexcitation_ratecoeff(float T_e, float clumpednne, double epsilon_trans,
                                                            double upperstatweight, double lowerstatweight,
                                                            int alltransindex) -> double;

[[gnu::pure]] [[nodiscard]] auto col_excitation_ratecoeff(float T_e, float clumpednne, double upperstatweight,
                                                          int alltransindex, double epsilon_trans,
                                                          double lowerstatweight) -> double;

// all four bound-bound rate coefficients of a single transition, as computed together by
// calculate_boundbound_ratecoeffs()
struct BoundBoundRatecoeffs {
  double rad_deexc;
  double rad_exc;
  double col_deexc;
  double col_exc;
};

[[gnu::pure]] [[nodiscard]] auto calculate_boundbound_ratecoeffs(int nonemptymgi, int alltransindex,
                                                                 double epsilon_trans, double upperstatweight,
                                                                 double lowerstatweight, double nnlevel_upper,
                                                                 double nnlevel_lower, float T_e, float clumpednne,
                                                                 double t_current) -> BoundBoundRatecoeffs;

// Sobolev-escape radiative deexcitation rate; multiply by the upper-level population to obtain a rate per second.
// Kromer & Sim (2009), Section 3.5.2, doi:10.1111/j.1365-2966.2009.15256.x.
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
  // tau_sobolev ~ zero (optically thin, escape probability beta -> 1) or negative (inverted level
  // populations, where the Sobolev beta would exceed one and is clamped to the thin-limit value here,
  // consistent with get_tau_sobolev() clamping tau >= 0 so that r-packets see a transparent line)
  return A_ul;
}

#endif  // MACROATOM_H
