#pragma once
#ifndef VPKT_H
#define VPKT_H

#include <cmath>
#include <span>

#include "artisoptions.h"
#include "constants.h"
#include "packet.h"
#include "vectors.h"

void read_parameterfile_vpkt();
void vpkt_init(int nts, int my_rank, bool continued_from_saved);
void vpkt_call_estimators(Packet &pkt, enum packet_type type_before_rpkt);
void vpkt_write_timestep(int nts, int my_rank, bool is_final);

void vpkt_remove_temp_file(int nts, int my_rank);

constexpr int VGRID_NY = 50;
constexpr int VGRID_NZ = 50;

// FREQUENCY
// dlognu = (log(numax) - log(numin)) / VMNUBINS ~ 3.9e-4 (10'000 over 1e14-5e15 Hz)
constexpr double VSPEC_NUMIN = CLIGHT / 10000 * 1e8;
constexpr double VSPEC_NUMAX = CLIGHT / 3500 * 1e8;
constexpr int VMNUBINS = 2500;

// TIME
// dlogt = (log(globals::tmin) - log(globals::tmax)) / VMTBINS ~ 3.69e-2 (111 over 2-120 d)
constexpr double VSPEC_TIMEMIN = 10 * DAY;
constexpr double VSPEC_TIMEMAX = 30 * DAY;
constexpr int VMTBINS = 30;

// number of virtual packets in a given timestep
inline int nvpkt{0};
inline int nvpkt_esc1{0};  // electron scattering event
inline int nvpkt_esc2{0};  // kpkt deactivation
inline int nvpkt_esc3{0};  // macroatom deactivation

inline double cell_is_optically_thick_vpkt;

[[nodiscard]] constexpr auto rot_angle(std::span<const double, 3> n1, std::span<const double, 3> n2,
                                       std::span<double, 3> ref1, std::span<double, 3> ref2) -> double {
  // Rotation angle from the scattering plane
  // We need to rotate Stokes Parameters to (or from) the scattering plane from (or to)
  // the meridian frame such that Q=1 is in the scattering plane and along ref1

  // ref1_sc is the ref1 axis in the scattering plane ref1 = n1 x ( n1 x n2 )
  const double n1_dot_n2 = dot(n1, n2);
  auto ref1_sc = std::array<double, 3>{n1[0] * n1_dot_n2 - n2[0], n1[1] * n1_dot_n2 - n2[1], n1[2] * n1_dot_n2 - n2[2]};
  ref1_sc = vec_norm(ref1_sc);

  double cos_stokes_rot_1 = dot(ref1_sc, ref1);
  const double cos_stokes_rot_2 = dot(ref1_sc, ref2);

  if (cos_stokes_rot_1 < -1) {
    cos_stokes_rot_1 = -1;
  }
  if (cos_stokes_rot_1 > 1) {
    cos_stokes_rot_1 = 1;
  }

  double i = 0;
  if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 > 0)) {
    i = acos(cos_stokes_rot_1);
  } else if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 > 0)) {
    i = PI - acos(fabs(cos_stokes_rot_1));
  } else if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 < 0)) {
    i = 2 * PI - acos(cos_stokes_rot_1);
  } else if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 < 0)) {
    i = PI + acos(fabs(cos_stokes_rot_1));
  }
  if (cos_stokes_rot_1 == 0) {
    i = PI / 2.;
  }
  if (cos_stokes_rot_2 == 0) {
    i = 0.;
  }

  return i;
}

// Routine to compute the meridian frame axes ref1 and ref2
[[nodiscard]] constexpr auto meridian(std::span<const double, 3> n)
    -> std::tuple<std::array<double, 3>, std::array<double, 3>> {
  // for ref_1 use (from triple product rule)
  const double n_xylen = std::sqrt(n[0] * n[0] + n[1] * n[1]);
  const auto ref1 =
      std::array<double, 3>{-1. * n[0] * n[2] / n_xylen, -1. * n[1] * n[2] / n_xylen, (1 - (n[2] * n[2])) / n_xylen};

  // for ref_2 use vector product of n_cmf with ref1
  const auto ref2 = cross_prod(ref1, n);
  return {ref1, ref2};
}

[[nodiscard]] constexpr auto lorentz(std::span<const double, 3> e_rf, std::span<const double, 3> n_rf,
                                     std::span<const double, 3> v) -> std::array<double, 3> {
  // Use Lorentz transformations to get e_cmf from e_rf

  const auto beta = std::array<const double, 3>{v[0] / CLIGHT, v[1] / CLIGHT, v[2] / CLIGHT};
  const double vsqr = dot(beta, beta);

  const double gamma_rel = 1. / (sqrt(1 - vsqr));

  const auto e_par = std::array<const double, 3>{dot(e_rf, beta) * beta[0] / (vsqr), dot(e_rf, beta) * beta[1] / (vsqr),
                                                 dot(e_rf, beta) * beta[2] / (vsqr)};

  const auto e_perp = std::array<const double, 3>{e_rf[0] - e_par[0], e_rf[1] - e_par[1], e_rf[2] - e_par[2]};

  const auto b_rf = cross_prod(n_rf, e_rf);

  // const double b_par[3] = {dot(b_rf, beta) * beta[0] / (vsqr), dot(b_rf, beta) * beta[1] / (vsqr),
  //                          dot(b_rf, beta) * beta[2] / (vsqr)};

  // const double b_perp[3] = {b_rf[0] - b_par[0], b_rf[1] - b_par[1], b_rf[2] - b_par[2]};

  const auto v_cr_b = cross_prod(beta, b_rf);

  // const double v_cr_e[3] = {beta[1] * e_rf[2] - beta[2] * e_rf[1], beta[2] * e_rf[0] - beta[0] * e_rf[2],
  //                           beta[0] * e_rf[1] - beta[1] * e_rf[0]};

  auto e_cmf = std::array<double, 3>{e_par[0] + gamma_rel * (e_perp[0] + v_cr_b[0]),
                                     e_par[1] + gamma_rel * (e_perp[1] + v_cr_b[1]),
                                     e_par[2] + gamma_rel * (e_perp[2] + v_cr_b[2])};
  return vec_norm(e_cmf);
}

// Routine to transform the Stokes Parameters from RF to CMF
constexpr auto frame_transform(std::span<const double, 3> n_rf, double *Q, double *U,
                               std::span<const double, 3> v) -> std::array<double, 3> {
  // Meridian frame in the RF
  const auto [ref1_rf, ref2_rf] = meridian(n_rf);

  const double Q0 = *Q;
  const double U0 = *U;

  // Compute polarisation (which is invariant)
  const double p = sqrt(Q0 * Q0 + U0 * U0);

  // We want to compute the angle between ref1 and the electric field
  double rot_angle = 0;

  if (p > 0) {
    const double cos2rot_angle = Q0 / p;
    const double sin2rot_angle = U0 / p;

    if ((cos2rot_angle > 0) && (sin2rot_angle > 0)) {
      rot_angle = acos(Q0 / p) / 2.;
    } else if ((cos2rot_angle < 0) && (sin2rot_angle > 0)) {
      rot_angle = (PI - acos(fabs(cos2rot_angle))) / 2.;
    } else if ((cos2rot_angle < 0) && (sin2rot_angle < 0)) {
      rot_angle = (PI + acos(fabs(cos2rot_angle))) / 2.;
    } else if ((cos2rot_angle > 0) && (sin2rot_angle < 0)) {
      rot_angle = (2. * PI - acos(fabs(cos2rot_angle))) / 2.;
    } else if (cos2rot_angle == 0) {
      rot_angle = 0.25 * PI;
      if (U0 < 0) {
        rot_angle = 0.75 * PI;
      }
    }
    if (sin2rot_angle == 0) {
      rot_angle = 0.;
      if (Q0 < 0) {
        rot_angle = 0.5 * PI;
      }
    }
  }

  // Define electric field by linear combination of ref1 and ref2 (using the angle just computed)

  const auto elec_rf = std::array<double, 3>{cos(rot_angle) * ref1_rf[0] - sin(rot_angle) * ref2_rf[0],
                                             cos(rot_angle) * ref1_rf[1] - sin(rot_angle) * ref2_rf[1],
                                             cos(rot_angle) * ref1_rf[2] - sin(rot_angle) * ref2_rf[2]};

  // Aberration
  const auto n_cmf = angle_ab(n_rf, v);

  // Lorentz transformation of E
  const auto elec_cmf = lorentz(elec_rf, n_rf, v);

  // Meridian frame in the CMF
  const auto [ref1_cmf, ref2_cmf] = meridian(n_cmf);

  // Projection of E onto ref1 and ref2
  const double cosine_elec_ref1 = dot(elec_cmf, ref1_cmf);
  const double cosine_elec_ref2 = dot(elec_cmf, ref2_cmf);

  // Compute the angle between ref1 and the electric field
  double theta_rot = 0.;
  if ((cosine_elec_ref1 > 0) && (cosine_elec_ref2 < 0)) {
    theta_rot = acos(cosine_elec_ref1);
  } else if ((cosine_elec_ref1 < 0) && (cosine_elec_ref2 > 0)) {
    theta_rot = PI + acos(fabs(cosine_elec_ref1));
  } else if ((cosine_elec_ref1 < 0) && (cosine_elec_ref2 < 0)) {
    theta_rot = PI - acos(fabs(cosine_elec_ref1));
  } else if ((cosine_elec_ref1 > 0) && (cosine_elec_ref2 > 0)) {
    theta_rot = 2 * PI - acos(cosine_elec_ref1);
  }
  if (cosine_elec_ref1 == 0) {
    theta_rot = PI / 2.;
  }
  if (cosine_elec_ref2 == 0) {
    theta_rot = 0.;
  }
  if (cosine_elec_ref1 > 1) {
    theta_rot = 0.;
  }
  if (cosine_elec_ref1 < -1) {
    theta_rot = PI;
  }

  // Compute Stokes Parameters in the CMF
  *Q = cos(2 * theta_rot) * p;
  *U = sin(2 * theta_rot) * p;

  return n_cmf;
}

#endif  // VPKT_H
