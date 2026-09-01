// Inline 3-vector geometry and special relativity helpers: dot/cross products, Doppler
// factors, angle aberration, and observer-frame/comoving-frame transformations.

#ifndef VECTORS_H
#define VECTORS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>

#include "artisoptions.h"
#include "constants.h"
#include "exspec.h"
#include "mpi_logging.h"
#include "packet.h"
#include "random.h"

// return the magnitude of a vector
template <size_t VECDIM>
[[gnu::pure]] [[nodiscard]] constexpr auto vec_len(const std::array<double, VECDIM>& vec) -> double {
  double squaredlen = 0.;
  for (auto i = 0ZU; i < VECDIM; i++) {
    squaredlen += pow2(vec[i]);
  }
  return std::sqrt(squaredlen);
}

// get a normalized copy of vec_in
[[gnu::pure]] [[nodiscard]] constexpr auto vec_norm(const Vec3d& vec_in) -> Vec3d {
  const double magnitude = vec_len(vec_in);
  assert_testmodeonly(magnitude > 0.);

  return Vec3d{vec_in[0] / magnitude, vec_in[1] / magnitude, vec_in[2] / magnitude};
}

// vector dot product
template <size_t VECDIM>
[[gnu::pure]] [[nodiscard]] constexpr auto dot(const std::array<double, VECDIM>& x, const std::array<double, VECDIM>& y)
    -> double {
  double sum = 0.;
  for (auto i = 0ZU; i < VECDIM; i++) {
    sum += x[i] * y[i];
  }
  return sum;
}

// Get velocity vector of the flow at a position with homologous expansion.
[[gnu::pure]] [[nodiscard]] constexpr auto get_velocity(const Vec3d& x, const double t) -> Vec3d {
  return Vec3d{x[0] / t, x[1] / t, x[2] / t};
}

[[gnu::pure]] [[nodiscard]] constexpr auto cross_prod(const Vec3d& vec_a, const Vec3d& vec_b) -> Vec3d {
  return Vec3d{
      (vec_a[1] * vec_b[2]) - (vec_b[1] * vec_a[2]),
      (vec_a[2] * vec_b[0]) - (vec_b[2] * vec_a[0]),
      (vec_a[0] * vec_b[1]) - (vec_b[0] * vec_a[1]),
  };
}

[[gnu::pure]] [[nodiscard]] constexpr auto vec_scale(const Vec3d& vec, const double scalefactor) -> Vec3d {
  return Vec3d{vec[0] * scalefactor, vec[1] * scalefactor, vec[2] * scalefactor};
}

// aberration of angles in special relativity
//   dir1: direction unit vector in frame1
//   vel: velocity of frame2 relative to frame1
//   dir2: direction vector in frame2
[[gnu::pure]] [[nodiscard]] constexpr auto angle_ab(const Vec3d& dir1, const Vec3d& vel) -> Vec3d {
  const double vsqr = dot(vel, vel) / CLIGHTSQUARED;
  const double gamma_rel = 1. / std::sqrt(1 - vsqr);

  const double ndotv = dot(dir1, vel);
  const double fact1 = gamma_rel * (1 - (ndotv / CLIGHT));
  const double fact2 = (gamma_rel - (pow2(gamma_rel) * ndotv / (gamma_rel + 1) / CLIGHT)) / CLIGHT;

  return vec_norm({
      (dir1[0] - (vel[0] * fact2)) / fact1,
      (dir1[1] - (vel[1] * fact2)) / fact1,
      (dir1[2] - (vel[2] * fact2)) / fact1,
  });
}

// Doppler factor either to first order v/c or fully relativisitic depending on USE_RELATIVISTIC_DOPPLER_SHIFT
// Arguments:
//   pos_rf: the rest frame position of the packet
//   dir_rf: the rest frame direction (unit vector) of light propagation
//   prop_time: the propagation time of the packet
// Returns: the ratio f = nu_cmf / nu_rf
[[nodiscard]] DEVICE_FUNC constexpr auto calculate_doppler_nucmf_on_nurf(const Vec3d& pos_rf, const Vec3d& dir_rf,
                                                                         const double prop_time) -> double {
  // velocity of the comoving frame relative to the rest frame
  const auto vel_rf = get_velocity(pos_rf, prop_time);

  assert_testmodeonly(dot(vel_rf, vel_rf) / CLIGHTSQUARED >= 0.);
  assert_testmodeonly(dot(vel_rf, vel_rf) / CLIGHTSQUARED < 1.);

  const double ndotv = dot(dir_rf, vel_rf);
  double dopplerfactor = 1. - (ndotv / CLIGHT);

  if (USE_RELATIVISTIC_DOPPLER_SHIFT) {
    const double betasq = dot(vel_rf, vel_rf) / CLIGHTSQUARED;
    assert_testmodeonly(betasq >= 0.);  // v < c
    assert_testmodeonly(betasq < 1.);  // v < c
    dopplerfactor = dopplerfactor / std::sqrt(1 - betasq);
  }

  assert_testmodeonly(std::isfinite(dopplerfactor));
  assert_testmodeonly(dopplerfactor > 0);

  return dopplerfactor;
}

// Move a packet along a straight line (specified by current dir vector). The distance moved is in the rest frame.
constexpr void move_pkt_withtime(Vec3d& pos_rf, const Vec3d& dir_rf, double& prop_time, const double nu_rf,
                                 double& nu_cmf, const double e_rf, double& e_cmf, const double distance) {
  assert_always(distance >= 0);

  const double nu_cmf_old = nu_cmf;
  prop_time += distance / CLIGHT_PROP;

  pos_rf = {pos_rf[0] + (dir_rf[0] * distance), pos_rf[1] + (dir_rf[1] * distance), pos_rf[2] + (dir_rf[2] * distance)};

  // During motion, rest frame energy and frequency are conserved. But need to update the co-moving ones.
  const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pos_rf, dir_rf, prop_time);

  // frequency should only ever decrease during packet movement with homologous expansion
  // enforce this to overcome numerical error
  nu_cmf = std::min(nu_rf * dopplerfactor, nu_cmf_old);

  e_cmf = e_rf * dopplerfactor;
}

constexpr void move_pkt_withtime(Packet& pkt, const double distance) {
  move_pkt_withtime(pkt.pos, pkt.dir, pkt.prop_time, pkt.nu_rf, pkt.nu_cmf, pkt.e_rf, pkt.e_cmf, distance);
}

// Set the packet's rest-frame frequency and energy from its co-moving frame values using the
// Doppler factor for its current position, direction, and propagation time.
DEVICE_FUNC constexpr void set_pkt_restframe_from_cmf(Packet& pkt) {
  const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
  pkt.nu_rf = pkt.nu_cmf / dopplerfactor;
  pkt.e_rf = pkt.e_cmf / dopplerfactor;
}

[[gnu::pure]] [[nodiscard]] constexpr auto get_escapedirectionbin(const Vec3d& dir_in) -> int {
  constexpr auto xhat = Vec3d{1.0, 0.0, 0.0};

  // sometimes dir vectors aren't accurately normalised
  const double dirmag = vec_len(dir_in);
  const auto dir = Vec3d{dir_in[0] / dirmag, dir_in[1] / dirmag, dir_in[2] / dirmag};

  // Angle resolved case: need to work out the correct angle bin
  const double costheta = dot(dir, syn_dir);
  const int costhetabin = std::clamp(static_cast<int>((costheta + 1.0) * NCOSTHETABINS / 2.0), 0, NCOSTHETABINS - 1);

  const auto vec1 = cross_prod(dir, syn_dir);
  constexpr auto vec2 = cross_prod(xhat, syn_dir);
  const double vec1_len = vec_len(vec1);
  const double cosphi = vec1_len > 1e-12 ? std::clamp(dot(vec1, vec2) / vec1_len, -1.0, 1.0) : 1.0;

  constexpr auto vec3 = cross_prod(vec2, syn_dir);
  const double testphi = dot(vec1, vec3);

  // With syn_dir = z and phi defined by dir = (sin(theta) cos(phi), sin(theta) sin(phi), cos(theta)),
  // the mapped angle is phi + pi for phi in (0, pi) and 2 pi - phi for phi in (pi, 2 pi).
  // The bins 0..(NPHIBINS/2 - 1) therefore cover phi from 2 pi down to pi in decreasing order, and
  // the bins NPHIBINS/2..(NPHIBINS - 1) cover phi from 0 up to pi in increasing order.
  // artistools contains the same mapping, so a change here breaks the analysis of the output files.
  // test_escapedirectionbin() in unittests.cc pins this mapping.
  const int phibin =
      std::clamp(static_cast<int>((testphi > 0 ? std::acos(cosphi) : std::acos(cosphi) + PI) / 2. / PI * NPHIBINS), 0,
                 NPHIBINS - 1);

  const int na = ((costhetabin * NPHIBINS) + phibin);

  return na;
}

// Assuming isotropic distribution, get a random direction vector
[[nodiscard]] DEVICE_FUNC inline auto get_rand_isotropic_unitvec(rngstate_type& rngstate) -> Vec3d {
  const double u = rng_uniform(rngstate);  // [0, 1)
  const double costheta = (2. * u) - 1.;
  const double sintheta = 2. * std::sqrt(u * (1. - u));

  const double phi = rng_uniform(rngstate) * 2 * PI;

  return {sintheta * std::cos(phi), sintheta * std::sin(phi), costheta};
}

// Rotation angle from the scattering plane
[[nodiscard]] constexpr auto get_rot_angle(const Vec3d& n1, const Vec3d& n2, const Vec3d& ref1, const Vec3d& ref2)
    -> double {
  // We need to rotate Stokes Parameters to (or from) the scattering plane from (or to)
  // the meridian frame such that Q=1 is in the scattering plane and along ref1

  // ref1_sc is the ref1 axis in the scattering plane ref1 = n1 x ( n1 x n2 )
  const double n1_dot_n2 = dot(n1, n2);
  const Vec3d ref1_sc_unnorm{(n1[0] * n1_dot_n2) - n2[0], (n1[1] * n1_dot_n2) - n2[1], (n1[2] * n1_dot_n2) - n2[2]};
  const double len = vec_len(ref1_sc_unnorm);
  if (len < 1e-12) {
    return 0.0;
  }
  const auto ref1_sc = Vec3d{ref1_sc_unnorm[0] / len, ref1_sc_unnorm[1] / len, ref1_sc_unnorm[2] / len};

  const double cos_stokes_rot_1 = std::clamp(dot(ref1_sc, ref1), -1., 1.);
  const double cos_stokes_rot_2 = dot(ref1_sc, ref2);

  const double rot_angle = std::atan2(cos_stokes_rot_2, cos_stokes_rot_1);
  return rot_angle < 0 ? rot_angle + (2 * PI) : rot_angle;
}

// Compute the meridian frame axes ref1 and ref2
[[gnu::pure]] [[nodiscard]] constexpr auto meridian(const Vec3d& dir) -> std::tuple<Vec3d, Vec3d> {
  // for ref_1 use (from triple product rule)
  const double n_xylen = std::sqrt(pow2(dir[0]) + pow2(dir[1]));
  if (n_xylen == 0.) {
    // if n is along z axis, we can just use x and y as the meridian frame axes.
    // Each vector has a name here. nvc++ 26.5 writes an invalid constant for a vector that goes
    // directly into the return value. The device compiler then stops.
    const Vec3d ref1_zaxis{1., 0., 0.};
    const Vec3d ref2_zaxis{0., 1., 0.};
    return {ref1_zaxis, ref2_zaxis};
  }
  const auto ref1 = Vec3d{-dir[0] * dir[2] / n_xylen, -dir[1] * dir[2] / n_xylen, (1 - pow2(dir[2])) / n_xylen};

  // for ref_2 use vector product of n_cmf with ref1
  const auto ref2 = cross_prod(ref1, dir);
  return {ref1, ref2};
}

[[gnu::pure]] [[nodiscard]] constexpr auto lorentz(const Vec3d& elec_rf, const Vec3d& n_rf, const Vec3d& v) -> Vec3d {
  // Use Lorentz transformations to get elec_cmf from elec_rf

  const Vec3d beta{v[0] / CLIGHT, v[1] / CLIGHT, v[2] / CLIGHT};
  const double betasquared = dot(beta, beta);
  if (betasquared == 0.) {
    return elec_rf;
  }

  const double gamma_rel = 1. / sqrt(1 - betasquared);
  const double elec_rf_dot_beta = dot(elec_rf, beta);

  const Vec3d elec_par{
      elec_rf_dot_beta * beta[0] / betasquared,
      elec_rf_dot_beta * beta[1] / betasquared,
      elec_rf_dot_beta * beta[2] / betasquared,
  };

  const Vec3d elec_perp{elec_rf[0] - elec_par[0], elec_rf[1] - elec_par[1], elec_rf[2] - elec_par[2]};

  const auto b_rf = cross_prod(n_rf, elec_rf);

  const auto v_cross_b = cross_prod(beta, b_rf);

  const auto elec_cmf = vec_norm({
      elec_par[0] + (gamma_rel * (elec_perp[0] + v_cross_b[0])),
      elec_par[1] + (gamma_rel * (elec_perp[1] + v_cross_b[1])),
      elec_par[2] + (gamma_rel * (elec_perp[2] + v_cross_b[2])),
  });
  return elec_cmf;
}

// Transform a direction and Stokes Parameters from RF to CMF
constexpr auto frame_transform(const Vec3d& n_rf, const double q0, const double u0, const Vec3d& v)
    -> std::tuple<Vec3d, double, double> {
  // Meridian frame in the RF
  const auto [ref1_rf, ref2_rf] = meridian(n_rf);

  // Compute polarisation (which is invariant)
  const double p = std::sqrt(pow2(q0) + pow2(u0));

  // We want to compute the angle between ref1 and the electric field
  double rot_angle = 0;

  if (p > 0) {
    const double pol_angle = std::atan2(u0, q0);
    rot_angle = (pol_angle < 0 ? pol_angle + (2. * PI) : pol_angle) / 2.;
  }

  const double cos_rot_angle = std::cos(rot_angle);
  const double sin_rot_angle = std::sin(rot_angle);

  // Define electric field by linear combination of ref1 and ref2 (using the angle just computed)
  const auto elec_rf = Vec3d{
      (cos_rot_angle * ref1_rf[0]) - (sin_rot_angle * ref2_rf[0]),
      (cos_rot_angle * ref1_rf[1]) - (sin_rot_angle * ref2_rf[1]),
      (cos_rot_angle * ref1_rf[2]) - (sin_rot_angle * ref2_rf[2]),
  };

  // Aberration
  const auto n_cmf = angle_ab(n_rf, v);

  // Lorentz transformation of E
  const auto elec_cmf = lorentz(elec_rf, n_rf, v);

  // Meridian frame in the CMF
  const auto [ref1_cmf, ref2_cmf] = meridian(n_cmf);

  // Projection of E onto ref1 and ref2
  const double cosine_elec_ref1 = dot(elec_cmf, ref1_cmf);
  const double cosine_elec_ref2 = dot(elec_cmf, ref2_cmf);

  double theta_rot = std::atan2(-cosine_elec_ref2, cosine_elec_ref1);
  if (theta_rot < 0) {
    theta_rot += 2 * PI;
  }

  // Compute Stokes Parameters in the CMF
  const auto q_cmf = cos(2 * theta_rot) * p;
  const auto u_cmf = sin(2 * theta_rot) * p;

  return {n_cmf, q_cmf, u_cmf};
}

// Compute the new Stokes Parameters after scattering and transform them back to the RF.
// Return a tuple of the new direction in the RF, the new q and u in the RF and the scattering phase-function
// probability pn
constexpr auto scatter_polarisation_to_rf(const Vec3d& old_dir_cmf, const Vec3d& new_dir_cmf, const double q_i_cmf,
                                          const double u_i_cmf, const Vec3d& vel_vec)
    -> std::tuple<Vec3d, double, double, double> {
  const auto [ref1_olddir, ref2_olddir] = meridian(old_dir_cmf);

  // This is the i1 angle of Bulla+2015, obtained by computing the angle between the
  // reference axes ref1 and ref2 in the meridian frame and the corresponding axes
  // ref1_sc and ref2_sc in the scattering plane.
  const double i1 = get_rot_angle(old_dir_cmf, new_dir_cmf, ref1_olddir, ref2_olddir);
  const double cos2i1 = cos(2 * i1);
  const double sin2i1 = sin(2 * i1);

  const double q_old = (q_i_cmf * cos2i1) - (u_i_cmf * sin2i1);
  const double u_old = (q_i_cmf * sin2i1) + (u_i_cmf * cos2i1);

  // Scattering

  const double mu = dot(old_dir_cmf, new_dir_cmf);
  const double musquared = pow2(mu);

  const double I_new = 0.75 * ((musquared + 1.) + (q_old * (musquared - 1.)));
  const double q_new = (0.75 * ((musquared - 1.) + (q_old * (musquared + 1.)))) / I_new;
  const double u_new = (1.5 * mu * u_old) / I_new;

  // Need to rotate Stokes Parameters out of the scattering plane to the meridian frame (Clockwise rotation of PI-i2)

  const auto [ref1, ref2] = meridian(new_dir_cmf);

  // This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
  // reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
  // meridian frame. NB: we need to add PI to transform THETA to i2
  const double i2 = PI + get_rot_angle(new_dir_cmf, old_dir_cmf, ref1, ref2);
  const double cos2i2 = cos(2 * i2);
  const double sin2i2 = sin(2 * i2);

  const double q_cmf = (q_new * cos2i2) + (u_new * sin2i2);
  const double u_cmf = (-q_new * sin2i2) + (u_new * cos2i2);

  const auto [new_dir_rf, q_rf, u_rf] =
      frame_transform(new_dir_cmf, q_cmf, u_cmf, Vec3d{-vel_vec[0], -vel_vec[1], -vel_vec[2]});

  const double pn = 3. / (16. * PI) * (1. + musquared + ((musquared - 1.) * q_old));
  return {new_dir_rf, q_rf, u_rf, pn};
}

#endif  // VECTORS_H
