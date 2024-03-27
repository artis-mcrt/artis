#pragma once
#ifndef VPKT_H
#define VPKT_H

#include <cmath>
#include <span>

#include "artisoptions.h"
#include "constants.h"
#include "packet.h"
#include "sn3d.h"
#include "vectors.h"

auto frame_transform(std::span<const double, 3> n_rf, double *Q, double *U,
                     std::span<const double, 3> v) -> std::array<double, 3>;

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

  assert_always(std::isfinite(i));

  return i;
}

// Routine to compute the meridian frame axes ref1 and ref2
[[nodiscard]] constexpr auto meridian(std::span<const double, 3> n,
                                      std::span<double, 3> ref1) -> std::array<double, 3> {
  // for ref_1 use (from triple product rule)
  const double n_xylen = std::sqrt(n[0] * n[0] + n[1] * n[1]);
  ref1[0] = -1. * n[0] * n[2] / n_xylen;
  ref1[1] = -1. * n[1] * n[2] / n_xylen;
  ref1[2] = (1 - (n[2] * n[2])) / n_xylen;

  // for ref_2 use vector product of n_cmf with ref1
  const auto ref2 = cross_prod(ref1, n);
  return ref2;
}

#endif  // VPKT_H
