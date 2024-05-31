#pragma once
#ifndef RPKT_H
#define RPKT_H

#include <ctime>
#include <vector>

struct Rpkt_continuum_absorptioncoeffs {
  double nu{-1.};  // frequency at which opacity was calculated
  double total{0.};
  double ffescat{0.};
  double ffheat{0.};
  double bf{0.};
};

struct Phixslist {
  std::vector<double> groundcont_gamma_contr;  // for either USE_LUT_PHOTOION = true or !USE_LUT_BFHEATING = false
  std::vector<double> chi_bf_sum;
  std::vector<double> gamma_contr;  // needed for DETAILED_BF_ESTIMATORS_ON
  int allcontend{-1};
  int allcontbegin{0};
  int bfestimend{-1};
  int bfestimbegin{0};
};

#include "artisoptions.h"
#include "atomic.h"
#include "sn3d.h"

void do_rpkt(Packet &pkt, double t2);
void emit_rpkt(Packet &pkt);
[[nodiscard]] auto closest_transition(double nu_cmf, int next_trans) -> int;
void calculate_chi_rpkt_cont(double nu_cmf, Rpkt_continuum_absorptioncoeffs &chi_rpkt_cont, Phixslist *phixslist,
                             int modelgridindex);
[[nodiscard]] auto sample_planck_times_expansion_opacity(int nonemptymgi) -> double;
void allocate_expansionopacities();
void calculate_expansion_opacities(int modelgridindex);
void MPI_Bcast_binned_opacities(int modelgridindex, int root_node_id);

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

[[nodiscard]] inline auto get_ionestimindex_nonemptymgi(const int nonemptymgi, const int element,
                                                        const int ion) -> int {
  assert_testmodeonly(ion >= 0);
  assert_testmodeonly(ion < get_nions(element) - 1);
  const int groundcontindex = globals::elements[element].ions[ion].groundcontindex;
  assert_always(groundcontindex >= 0);
  return nonemptymgi * globals::nbfcontinua_ground + groundcontindex;
}

#endif  // RPKT_H
