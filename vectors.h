#ifndef VECTORS_H
#define VECTORS_H

#include <array>
#include <cmath>
#include <numeric>
#include <span>

#include "constants.h"
#include "exspec.h"
#include "packet.h"
#include "sn3d.h"

void scatter_dir(std::span<const double, 3> dir_in, double cos_theta, std::span<double, 3> dir_out);
void get_rand_isotropic_unitvec(std::span<double, 3> vecout);

[[nodiscard]] [[gnu::pure]] constexpr auto vec_len(std::span<const double> vec) -> double
// return the the magnitude of a vector
{
  const double squaredlen = std::accumulate(vec.begin(), vec.end(), 0., [](auto a, auto b) { return a + b * b; });

  return std::sqrt(squaredlen);
}

constexpr void vec_norm(std::span<const double, 3> vec_in, std::span<double, 3> vec_out)
// normalizing a copy of vec_in and save it to vec_out
{
  const double magnitude = vec_len(vec_in);

  vec_out[0] = vec_in[0] / magnitude;
  vec_out[1] = vec_in[1] / magnitude;
  vec_out[2] = vec_in[2] / magnitude;

  assert_testmodeonly(fabs(vec_len(vec_out) - 1.) < 1.e-10);
}

[[nodiscard]] [[gnu::pure]] constexpr auto dot(std::span<const double> x, std::span<const double> y) -> double
// vector dot product
{
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.);
}

constexpr void get_velocity(std::span<const double, 3> x, std::span<double, 3> y, const double t)
// Routine for getting velocity vector of the flow at a position with homologous expansion.
{
  y[0] = x[0] / t;
  y[1] = x[1] / t;
  y[2] = x[2] / t;
}

constexpr std::array<double, 3> cross_prod(std::span<const double, 3> vec1, std::span<const double, 3> vec2) {
  std::array<double, 3> vecout = {(vec1[1] * vec2[2]) - (vec2[1] * vec1[2]), (vec1[2] * vec2[0]) - (vec2[2] * vec1[0]),
                                  (vec1[0] * vec2[1]) - (vec2[0] * vec1[1])};
  return vecout;
}

constexpr void vec_scale(std::span<double, 3> vec, const double scalefactor) {
  vec[0] *= scalefactor;
  vec[1] *= scalefactor;
  vec[2] *= scalefactor;
}

constexpr void vec_copy(std::span<double, 3> destination, std::span<const double, 3> source) {
  destination[0] = source[0];
  destination[1] = source[1];
  destination[2] = source[2];
}

constexpr void angle_ab(std::span<const double, 3> dir1, std::span<const double, 3> vel, std::span<double, 3> dir2)
// aberation of angles in special relativity
//   dir1: direction unit vector in frame1
//   vel: velocity of frame2 relative to frame1
//   dir2: direction vector in frame2
{
  const double vsqr = dot(vel, vel) / CLIGHTSQUARED;
  const double gamma_rel = 1. / std::sqrt(1 - vsqr);

  const double ndotv = dot(dir1, vel);
  const double fact1 = gamma_rel * (1 - (ndotv / CLIGHT));
  const double fact2 = (gamma_rel - (gamma_rel * gamma_rel * ndotv / (gamma_rel + 1) / CLIGHT)) / CLIGHT;

  for (int d = 0; d < 3; d++) {
    dir2[d] = (dir1[d] - (vel[d] * fact2)) / fact1;
  }

  vec_norm(dir2, dir2);
}

[[gnu::pure]] [[nodiscard]] constexpr auto doppler_nucmf_on_nurf(std::span<const double, 3> dir_rf,
                                                                 std::span<const double, 3> vel_rf) -> double
// Doppler factor
// arguments:
//   dir_rf: the rest frame direction (unit vector) of light propagation
//   vel_rf: velocity of the comoving frame relative to the rest frame
// returns: the ratio f = nu_cmf / nu_rf
{
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

[[gnu::pure]] [[nodiscard]] constexpr auto doppler_squared_nucmf_on_nurf(std::span<const double, 3> pos_rf,
                                                                         std::span<const double, 3> dir_rf,
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
  std::array<double, 3> vel_rf = {0, 0, 0};  // homologous flow velocity
  get_velocity(pos_rf, vel_rf, prop_time);

  assert_testmodeonly(dot(vel_rf, vel_rf) / CLIGHTSQUARED >= 0.);
  assert_testmodeonly(dot(vel_rf, vel_rf) / CLIGHTSQUARED < 1.);

  const double ndotv_on_c = dot(dir_rf, vel_rf) / CLIGHT;
  double dopplerfactorsq = USE_RELATIVISTIC_DOPPLER_SHIFT
                               ? std::pow(1. - ndotv_on_c, 2) / (1 - (dot(vel_rf, vel_rf) / CLIGHTSQUARED))
                               : (1. - 2 * ndotv_on_c);

  assert_testmodeonly(std::isfinite(dopplerfactorsq));
  assert_testmodeonly(dopplerfactorsq > 0);

  return dopplerfactorsq;
}

[[gnu::pure]] [[nodiscard]] constexpr auto doppler_packet_nucmf_on_nurf(std::span<const double, 3> pos_rf,
                                                                        std::span<const double, 3> dir_rf,
                                                                        const double prop_time) -> double {
  double flow_velocity[3] = {0, 0, 0};  // homologous flow velocity
  get_velocity(pos_rf, flow_velocity, prop_time);
  return doppler_nucmf_on_nurf(dir_rf, flow_velocity);
}

constexpr auto move_pkt_withtime(struct packet *pkt_ptr, const double distance) -> double
/// Subroutine to move a packet along a straight line (specified by current
/// dir vector). The distance moved is in the rest frame.
{
  assert_always(distance >= 0);

  const double nu_cmf_old = pkt_ptr->nu_cmf;
  pkt_ptr->prop_time += distance / CLIGHT_PROP;

  pkt_ptr->pos[0] += (pkt_ptr->dir[0] * distance);
  pkt_ptr->pos[1] += (pkt_ptr->dir[1] * distance);
  pkt_ptr->pos[2] += (pkt_ptr->dir[2] * distance);

  /// During motion, rest frame energy and frequency are conserved.
  /// But need to update the co-moving ones.
  const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr->pos, pkt_ptr->dir, pkt_ptr->prop_time);

  pkt_ptr->nu_cmf = pkt_ptr->nu_rf * dopplerfactor;
  pkt_ptr->e_cmf = pkt_ptr->e_rf * dopplerfactor;

  // frequency should only over decrease due to packet movement
  // enforce this to overcome numerical error
  pkt_ptr->nu_cmf = std::min(pkt_ptr->nu_cmf, nu_cmf_old);

  return dopplerfactor;
}

[[nodiscard]] [[gnu::pure]] constexpr auto get_arrive_time(const struct packet *const pkt_ptr) -> double
/// We know that a packet escaped at "escape_time". However, we have
/// to allow for travel time. Use the formula in Leon's paper. The extra
/// distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
{
  return pkt_ptr->escape_time - (dot(pkt_ptr->pos, pkt_ptr->dir) / CLIGHT_PROP);
}

constexpr auto get_escapedirectionbin(std::span<const double, 3> dir_in, std::span<const double, 3> syn_dir) -> int {
  constexpr std::array<double, 3> xhat = {1.0, 0.0, 0.0};

  // sometimes dir vectors aren't accurately normalised
  const double dirmag = vec_len(dir_in);
  std::array<const double, 3> dir = {dir_in[0] / dirmag, dir_in[1] / dirmag, dir_in[2] / dirmag};

  /// Angle resolved case: need to work out the correct angle bin
  const double costheta = dot(dir, syn_dir);
  const int costhetabin = static_cast<int>((costheta + 1.0) * NPHIBINS / 2.0);
  assert_testmodeonly(costhetabin < NCOSTHETABINS);

  const std::array<double, 3> vec1 = cross_prod(dir, syn_dir);

  const std::array<double, 3> vec2 = cross_prod(xhat, syn_dir);
  const double cosphi = dot(vec1, vec2) / vec_len(vec1) / vec_len(vec2);

  const std::array<double, 3> vec3 = cross_prod(vec2, syn_dir);
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

#endif  // VECTORS_H