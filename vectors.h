#pragma once
#ifndef VECTORS_H
#define VECTORS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <tuple>

#include "constants.h"
#include "exspec.h"
#include "packet.h"
#include "sn3d.h"

// return the the magnitude of a vector
template <size_t VECDIM>
[[nodiscard]] constexpr auto vec_len(const std::array<double, VECDIM> &vec) -> double {
  const double squaredlen = std::accumulate(vec.begin(), vec.end(), 0., [](auto a, auto b) { return a + b * b; });

  return std::sqrt(squaredlen);
}

// get a normalized copy of vec_in
[[nodiscard]] constexpr auto vec_norm(const std::array<double, 3> &vec_in) {
  const double magnitude = vec_len(vec_in);
  const std::array<double, 3> vec_out{vec_in[0] / magnitude, vec_in[1] / magnitude, vec_in[2] / magnitude};

  assert_testmodeonly(fabs(vec_len(vec_out) - 1.) < 1.e-10);
  return vec_out;
}

// vector dot product
template <size_t S1, size_t S2>
[[nodiscard]] constexpr auto dot(const std::array<double, S1> &x, const std::array<double, S2> &y) -> double {
  // if len(x) < len(y), the extra elements of y are ignored
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.);
}

// Get velocity vector of the flow at a position with homologous expansion.
[[nodiscard]] constexpr auto get_velocity(const std::array<double, 3> &x, const double t) -> std::array<double, 3> {
  return std::array<double, 3>{x[0] / t, x[1] / t, x[2] / t};
}

[[nodiscard]] constexpr auto cross_prod(const std::array<double, 3> &vec_a, const std::array<double, 3> &vec_b) {
  return std::array<double, 3>{(vec_a[1] * vec_b[2]) - (vec_b[1] * vec_a[2]),
                               (vec_a[2] * vec_b[0]) - (vec_b[2] * vec_a[0]),
                               (vec_a[0] * vec_b[1]) - (vec_b[0] * vec_a[1])};
}

[[nodiscard]] constexpr auto vec_scale(const std::array<double, 3> vec, const double scalefactor) {
  return std::array<double, 3>{vec[0] * scalefactor, vec[1] * scalefactor, vec[2] * scalefactor};
}

// aberation of angles in special relativity
//   dir1: direction unit vector in frame1
//   vel: velocity of frame2 relative to frame1
//   dir2: direction vector in frame2
[[nodiscard]] constexpr auto angle_ab(const std::array<double, 3> &dir1,
                                      const std::array<double, 3> &vel) -> std::array<double, 3> {
  const double vsqr = dot(vel, vel) / CLIGHTSQUARED;
  const double gamma_rel = 1. / std::sqrt(1 - vsqr);

  const double ndotv = dot(dir1, vel);
  const double fact1 = gamma_rel * (1 - (ndotv / CLIGHT));
  const double fact2 = (gamma_rel - (gamma_rel * gamma_rel * ndotv / (gamma_rel + 1) / CLIGHT)) / CLIGHT;

  const auto dir2 = std::array<double, 3>{(dir1[0] - (vel[0] * fact2)) / fact1, (dir1[1] - (vel[1] * fact2)) / fact1,
                                          (dir1[2] - (vel[2] * fact2)) / fact1};

  return vec_norm(dir2);
}

// Doppler factor
// arguments:
//   dir_rf: the rest frame direction (unit vector) of light propagation
//   vel_rf: velocity of the comoving frame relative to the rest frame
// returns: the ratio f = nu_cmf / nu_rf
[[nodiscard]] constexpr auto doppler_nucmf_on_nurf(const std::array<double, 3> &dir_rf,
                                                   const std::array<double, 3> &vel_rf) -> double {
  assert_testmodeonly(dot(vel_rf, vel_rf) / CLIGHTSQUARED >= 0.);
  assert_testmodeonly(dot(vel_rf, vel_rf) / CLIGHTSQUARED < 1.);

  const double ndotv = dot(dir_rf, vel_rf);
  double dopplerfactor = 1. - (ndotv / CLIGHT);

  if (USE_RELATIVISTIC_DOPPLER_SHIFT) {
    const double betasq = dot(vel_rf, vel_rf) / CLIGHTSQUARED;
    assert_testmodeonly(betasq >= 0.);  // v < c
    assert_testmodeonly(betasq < 1.);   // v < c
    dopplerfactor = dopplerfactor / std::sqrt(1 - betasq);
  }

  assert_testmodeonly(std::isfinite(dopplerfactor));
  assert_testmodeonly(dopplerfactor > 0);

  return dopplerfactor;
}

[[nodiscard]] constexpr auto doppler_squared_nucmf_on_nurf(const std::array<double, 3> &pos_rf,
                                                           const std::array<double, 3> &dir_rf,
                                                           const double prop_time) -> double
// Doppler factor squared, either to first order v/c or fully relativisitic
// depending on USE_RELATIVISTIC_DOPPLER_SHIFT
//
// arguments:
//   pos_rf: the rest frame position of the packet
//   dir_rf: the rest frame direction (unit vector) of light propagation
//   prop_time: the propagation time of the packet
// returns: the ratio f = (nu_cmf / nu_rf) ^ 2
{
  // velocity of the comoving frame relative to the rest frame
  const auto vel_rf = get_velocity(pos_rf, prop_time);

  assert_testmodeonly(dot(vel_rf, vel_rf) / CLIGHTSQUARED >= 0.);
  assert_testmodeonly(dot(vel_rf, vel_rf) / CLIGHTSQUARED < 1.);

  const double ndotv_on_c = dot(dir_rf, vel_rf) / CLIGHT;
  const double dopplerfactorsq = USE_RELATIVISTIC_DOPPLER_SHIFT
                                     ? std::pow(1. - ndotv_on_c, 2) / (1 - (dot(vel_rf, vel_rf) / CLIGHTSQUARED))
                                     : (1. - 2 * ndotv_on_c);

  assert_testmodeonly(std::isfinite(dopplerfactorsq));
  assert_testmodeonly(dopplerfactorsq > 0);

  return dopplerfactorsq;
}

[[nodiscard]] constexpr auto doppler_packet_nucmf_on_nurf(const std::array<double, 3> &pos_rf,
                                                          const std::array<double, 3> &dir_rf,
                                                          const double prop_time) -> double {
  return doppler_nucmf_on_nurf(dir_rf, get_velocity(pos_rf, prop_time));
}

// Move a packet along a straight line (specified by current dir vector). The distance moved is in the rest frame.
constexpr auto move_pkt_withtime(std::array<double, 3> &pos_rf, const std::array<double, 3> &dir_rf, double &prop_time,
                                 const double nu_rf, double &nu_cmf, const double e_rf, double &e_cmf,
                                 const double distance) -> double {
  assert_always(distance >= 0);

  const double nu_cmf_old = nu_cmf;
  prop_time += distance / CLIGHT_PROP;

  pos_rf[0] += (dir_rf[0] * distance);
  pos_rf[1] += (dir_rf[1] * distance);
  pos_rf[2] += (dir_rf[2] * distance);

  // During motion, rest frame energy and frequency are conserved.
  // But need to update the co-moving ones.
  const double dopplerfactor = doppler_packet_nucmf_on_nurf(pos_rf, dir_rf, prop_time);

  nu_cmf = nu_rf * dopplerfactor;
  e_cmf = e_rf * dopplerfactor;

  // frequency should only over decrease due to packet movement
  // enforce this to overcome numerical error
  nu_cmf = std::min(nu_cmf, nu_cmf_old);

  return dopplerfactor;
}

constexpr auto move_pkt_withtime(Packet &pkt, const double distance) -> double {
  return move_pkt_withtime(pkt.pos, pkt.dir, pkt.prop_time, pkt.nu_rf, pkt.nu_cmf, pkt.e_rf, pkt.e_cmf, distance);
}

[[nodiscard]] constexpr auto get_arrive_time(const Packet &pkt) -> double
// We know that a packet escaped at "escape_time". However, we have
// to allow for travel time. Use the formula in Leon's paper. The extra
// distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
{
  return pkt.escape_time - (dot(pkt.pos, pkt.dir) / CLIGHT_PROP);
}

[[nodiscard]] constexpr auto get_escapedirectionbin(const std::array<double, 3> &dir_in,
                                                    const std::array<double, 3> &syn_dir) -> int {
  constexpr auto xhat = std::array<double, 3>{1.0, 0.0, 0.0};

  // sometimes dir vectors aren't accurately normalised
  const double dirmag = vec_len(dir_in);
  const auto dir = std::array<double, 3>{dir_in[0] / dirmag, dir_in[1] / dirmag, dir_in[2] / dirmag};

  // Angle resolved case: need to work out the correct angle bin
  const double costheta = dot(dir, syn_dir);
  const int costhetabin = static_cast<int>((costheta + 1.0) * NPHIBINS / 2.0);
  assert_testmodeonly(costhetabin < NCOSTHETABINS);

  const auto vec1 = cross_prod(dir, syn_dir);

  const auto vec2 = cross_prod(xhat, syn_dir);
  const double cosphi = dot(vec1, vec2) / vec_len(vec1) / vec_len(vec2);

  const auto vec3 = cross_prod(vec2, syn_dir);
  const double testphi = dot(vec1, vec3);

  // with phi defined according to y = cos(theta) * sin(phi), the
  // phibins are in decreasing phi order (i.e. the upper side of bin zero 0 is 2pi)
  const int phibin = static_cast<int>((testphi >= 0 ? acos(cosphi) : acos(cosphi) + PI) / 2. / PI * NPHIBINS);

  assert_testmodeonly(phibin >= 0);
  assert_testmodeonly(phibin < NPHIBINS);
  const int na = static_cast<int>((costhetabin * NPHIBINS) + phibin);
  assert_always(na < MABINS);

  return na;
}

// Assuming isotropic distribution, get a random direction vector
[[nodiscard]] inline auto get_rand_isotropic_unitvec() -> std::array<double, 3> {
  const double costheta = -1 + (2. * rng_uniform());

  const double phi = rng_uniform() * 2 * PI;

  const double sintheta = std::sqrt(1. - (costheta * costheta));

  return std::array<double, 3>{sintheta * std::cos(phi), sintheta * std::sin(phi), costheta};
}

// Rotation angle from the scattering plane
[[nodiscard]] constexpr auto get_rot_angle(const std::array<double, 3> &n1, const std::array<double, 3> &n2,
                                           const std::array<double, 3> &ref1,
                                           const std::array<double, 3> &ref2) -> double {
  // We need to rotate Stokes Parameters to (or from) the scattering plane from (or to)
  // the meridian frame such that Q=1 is in the scattering plane and along ref1

  // ref1_sc is the ref1 axis in the scattering plane ref1 = n1 x ( n1 x n2 )
  const double n1_dot_n2 = dot(n1, n2);
  auto ref1_sc = std::array<double, 3>{n1[0] * n1_dot_n2 - n2[0], n1[1] * n1_dot_n2 - n2[1], n1[2] * n1_dot_n2 - n2[2]};
  ref1_sc = vec_norm(ref1_sc);

  const double cos_stokes_rot_1 = std::clamp(dot(ref1_sc, ref1), -1., 1.);
  const double cos_stokes_rot_2 = dot(ref1_sc, ref2);

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
[[nodiscard]] constexpr auto meridian(const std::array<double, 3> &n)
    -> std::tuple<std::array<double, 3>, std::array<double, 3>> {
  // for ref_1 use (from triple product rule)
  const double n_xylen = std::sqrt(n[0] * n[0] + n[1] * n[1]);
  const auto ref1 =
      std::array<double, 3>{-1. * n[0] * n[2] / n_xylen, -1. * n[1] * n[2] / n_xylen, (1 - (n[2] * n[2])) / n_xylen};

  // for ref_2 use vector product of n_cmf with ref1
  const auto ref2 = cross_prod(ref1, n);
  return {ref1, ref2};
}

[[nodiscard]] constexpr auto lorentz(const std::array<double, 3> &e_rf, const std::array<double, 3> &n_rf,
                                     const std::array<double, 3> &v) -> std::array<double, 3> {
  // Use Lorentz transformations to get e_cmf from e_rf

  const auto beta = std::array<double, 3>{v[0] / CLIGHT, v[1] / CLIGHT, v[2] / CLIGHT};
  const double vsqr = dot(beta, beta);

  const double gamma_rel = 1. / (sqrt(1 - vsqr));

  const auto e_par = std::array<double, 3>{dot(e_rf, beta) * beta[0] / (vsqr), dot(e_rf, beta) * beta[1] / (vsqr),
                                           dot(e_rf, beta) * beta[2] / (vsqr)};

  const auto e_perp = std::array<double, 3>{e_rf[0] - e_par[0], e_rf[1] - e_par[1], e_rf[2] - e_par[2]};

  const auto b_rf = cross_prod(n_rf, e_rf);

  // const double b_par[3] = {dot(b_rf, beta) * beta[0] / (vsqr), dot(b_rf, beta) * beta[1] / (vsqr),
  //                          dot(b_rf, beta) * beta[2] / (vsqr)};

  // const double b_perp[3] = {b_rf[0] - b_par[0], b_rf[1] - b_par[1], b_rf[2] - b_par[2]};

  const auto v_cr_b = cross_prod(beta, b_rf);

  // const double v_cr_e[3] = {beta[1] * e_rf[2] - beta[2] * e_rf[1], beta[2] * e_rf[0] - beta[0] * e_rf[2],
  //                           beta[0] * e_rf[1] - beta[1] * e_rf[0]};

  const auto e_cmf = std::array<double, 3>{e_par[0] + gamma_rel * (e_perp[0] + v_cr_b[0]),
                                           e_par[1] + gamma_rel * (e_perp[1] + v_cr_b[1]),
                                           e_par[2] + gamma_rel * (e_perp[2] + v_cr_b[2])};
  return vec_norm(e_cmf);
}

// Routine to transform the Stokes Parameters from RF to CMF
constexpr auto frame_transform(const std::array<double, 3> n_rf, double *Q, double *U,
                               const std::array<double, 3> v) -> std::array<double, 3> {
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

#endif  // VECTORS_H
