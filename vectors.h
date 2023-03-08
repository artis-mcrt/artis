#ifndef VECTORS_H
#define VECTORS_H

#include <cmath>

#include "constants.h"
#include "cuda.h"
#include "exspec.h"
#include "packet.h"
#include "sn3d.h"

__host__ __device__ void scatter_dir(const double dir_in[3], double cos_theta, double dir_out[3]);
__host__ __device__ void get_rand_isotropic_unitvec(double vecout[3]);

constexpr double vec_len(const double x[3])
// return the the magnitude of a vector
{
  return std::sqrt((x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]));
}

constexpr void vec_norm(const double vec_in[3], double vec_out[3])
// normalizing a copy of vec_in and save it to vec_out
{
  const double magnitude = vec_len(vec_in);

  vec_out[0] = vec_in[0] / magnitude;
  vec_out[1] = vec_in[1] / magnitude;
  vec_out[2] = vec_in[2] / magnitude;

  assert_testmodeonly(fabs(vec_len(vec_out) - 1.) < 1.e-10);
}

constexpr double dot(const double x[3], const double y[3])
// vector dot product
{
  return (x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]);
}

constexpr void get_velocity(const double x[3], double y[3], const double t)
// Routine for getting velocity vector of the flow at a position with homologous expansion.
{
  y[0] = x[0] / t;
  y[1] = x[1] / t;
  y[2] = x[2] / t;
}

constexpr void cross_prod(const double vec1[3], const double vec2[3], double vecout[3]) {
  vecout[0] = (vec1[1] * vec2[2]) - (vec2[1] * vec1[2]);
  vecout[1] = (vec1[2] * vec2[0]) - (vec2[2] * vec1[0]);
  vecout[2] = (vec1[0] * vec2[1]) - (vec2[0] * vec1[1]);
}

constexpr void vec_scale(double vec[3], const double scalefactor) {
  vec[0] *= scalefactor;
  vec[1] *= scalefactor;
  vec[2] *= scalefactor;
}

constexpr void vec_copy(double destination[3], const double source[3]) {
  destination[0] = source[0];
  destination[1] = source[1];
  destination[2] = source[2];
}

constexpr void angle_ab(const double dir1[3], const double vel[3], double dir2[3])
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

constexpr double doppler_nucmf_on_nurf(const double dir_rf[3], const double vel_rf[3])
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

constexpr double doppler_squared_nucmf_on_nurf(const double dir_rf[3], const double vel_rf[3])
// Doppler factor squared, either to first order v/c or fully relativisitic
// depending on USE_RELATIVISTIC_DOPPLER_SHIFT
//
// arguments:
//   dir_rf: the rest frame direction (unit vector) of light propagation
//   vel_rf: velocity of the comoving frame relative to the rest frame
// returns: the ratio f = (nu_cmf / nu_rf) ^ 2
{
  assert_testmodeonly(dot(vel_rf, vel_rf) / CLIGHTSQUARED >= 0.);
  assert_testmodeonly(dot(vel_rf, vel_rf) / CLIGHTSQUARED < 1.);

  const double ndotv = dot(dir_rf, vel_rf);
  double dopplerfactorsq;

  if (USE_RELATIVISTIC_DOPPLER_SHIFT) {
    const double betasq = dot(vel_rf, vel_rf) / CLIGHTSQUARED;
    assert_testmodeonly(betasq >= 0.);  // v < c
    assert_testmodeonly(betasq < 1.);   // v < c
    dopplerfactorsq = std::pow(1. - (ndotv / CLIGHT), 2) / (1 - betasq);
  } else {
    dopplerfactorsq = 1. - 2 * (ndotv / CLIGHT);
  }

  assert_testmodeonly(std::isfinite(dopplerfactorsq));
  assert_testmodeonly(dopplerfactorsq > 0);

  return dopplerfactorsq;
}

constexpr double doppler_packet_nucmf_on_nurf(const struct packet *const pkt_ptr) {
  double flow_velocity[3] = {0, 0, 0};  // homologous flow velocity
  get_velocity(pkt_ptr->pos, flow_velocity, pkt_ptr->prop_time);
  return doppler_nucmf_on_nurf(pkt_ptr->dir, flow_velocity);
}

constexpr void move_pkt(struct packet *pkt_ptr, const double distance)
/// Subroutine to move a packet along a straight line (specified by current
/// dir vector). The distance moved is in the rest frame.
{
  /// First update pos.
  assert_always(distance >= 0);

  pkt_ptr->pos[0] += (pkt_ptr->dir[0] * distance);
  pkt_ptr->pos[1] += (pkt_ptr->dir[1] * distance);
  pkt_ptr->pos[2] += (pkt_ptr->dir[2] * distance);

  /// During motion, rest frame energy and frequency are conserved.
  /// But need to update the co-moving ones.
  const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr);
  pkt_ptr->nu_cmf = pkt_ptr->nu_rf * dopplerfactor;
  pkt_ptr->e_cmf = pkt_ptr->e_rf * dopplerfactor;
}

constexpr void move_pkt_withtime(struct packet *pkt_ptr, const double distance)
/// Subroutine to move a packet along a straight line (specified by current
/// dir vector). The distance moved is in the rest frame.
{
  const double nu_cmf_old = pkt_ptr->nu_cmf;
  pkt_ptr->prop_time += distance / CLIGHT_PROP;
  move_pkt(pkt_ptr, distance);

  // frequency should only over decrease due to packet movement
  // enforce this to overcome numerical error
  if (pkt_ptr->nu_cmf > nu_cmf_old) {
    pkt_ptr->nu_cmf = nu_cmf_old;
  }
}

constexpr double get_arrive_time(const struct packet *pkt_ptr)
/// We know that a packet escaped at "escape_time". However, we have
/// to allow for travel time. Use the formula in Leon's paper. The extra
/// distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
{
  return pkt_ptr->escape_time - (dot(pkt_ptr->pos, pkt_ptr->dir) / CLIGHT_PROP);
}

inline double get_arrive_time_cmf(const struct packet *pkt_ptr) {
  return pkt_ptr->escape_time * std::sqrt(1. - (globals::vmax * globals::vmax / CLIGHTSQUARED));
}

constexpr int get_escapedirectionbin(const double dir_in[3], const double syn_dir[3]) {
  constexpr double xhat[3] = {1.0, 0.0, 0.0};

  // sometimes dir vectors aren't accurately normalised
  const double dirmag = vec_len(dir_in);
  const double dir[3] = {dir_in[0] / dirmag, dir_in[1] / dirmag, dir_in[2] / dirmag};

  /// Angle resolved case: need to work out the correct angle bin
  const double costheta = dot(dir, syn_dir);
  const int costhetabin = ((costheta + 1.0) * NPHIBINS / 2.0);
  assert_testmodeonly(costhetabin < NCOSTHETABINS);

  double vec1[3] = {0};
  cross_prod(dir, syn_dir, vec1);

  double vec2[3] = {0};
  cross_prod(xhat, syn_dir, vec2);
  const double cosphi = dot(vec1, vec2) / vec_len(vec1) / vec_len(vec2);

  double vec3[3] = {0};
  cross_prod(vec2, syn_dir, vec3);
  const double testphi = dot(vec1, vec3);

  int phibin = 0;
  if (testphi > 0) {
    phibin = (acos(cosphi) / 2. / PI * NPHIBINS);
  } else {
    phibin = ((acos(cosphi) + PI) / 2. / PI * NPHIBINS);
  }
  assert_testmodeonly(phibin < NPHIBINS);
  const int na = (costhetabin * NPHIBINS) + phibin;
  assert_always(na < MABINS);

  return na;
}

#endif  // VECTORS_H