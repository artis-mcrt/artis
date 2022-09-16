#ifndef VECTORS_H
#define VECTORS_H

#include <gsl/gsl_blas.h>

#include <cmath>

#include "artisoptions.h"
#include "constants.h"
#include "cuda.h"
#include "packet.h"

__host__ __device__ void scatter_dir(const double dir_in[3], double cos_theta, double dir_out[3]);
__host__ __device__ void get_rand_isotropic_unitvec(double vecout[3]);
__host__ __device__ void move_pkt_withtime(struct packet *pkt_ptr, double distance);

__host__ __device__ constexpr double vec_len(const double x[3])
// return the the magnitude of a vector
{
  return sqrt((x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]));
}

__host__ __device__ constexpr void vec_norm(const double vec_in[3], double vec_out[3])
// Routine for normalizing a vector.
{
  const double magnitude = vec_len(vec_in);

  vec_out[0] = vec_in[0] / magnitude;
  vec_out[1] = vec_in[1] / magnitude;
  vec_out[2] = vec_in[2] / magnitude;

  // vec_copy(vec_out, vec_in);
}

__host__ __device__ constexpr double dot(const double x[3], const double y[3])
// Routine for taking dot product.
{
  return (x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]);
}

__host__ __device__ constexpr void get_velocity(const double x[3], double y[3], const double t)
// Routine for getting velocity vector of the flow at a position with homologous expansion.
{
  y[0] = x[0] / t;
  y[1] = x[1] / t;
  y[2] = x[2] / t;
}

__host__ __device__ constexpr void cross_prod(const double vec1[3], const double vec2[3], double vecout[3]) {
  vecout[0] = (vec1[1] * vec2[2]) - (vec2[1] * vec1[2]);
  vecout[1] = (vec1[2] * vec2[0]) - (vec2[2] * vec1[0]);
  vecout[2] = (vec1[0] * vec2[1]) - (vec2[0] * vec1[1]);
}

__host__ __device__ constexpr void vec_scale(double vec[3], const double scalefactor) {
  vec[0] *= scalefactor;
  vec[1] *= scalefactor;
  vec[2] *= scalefactor;
}

__host__ __device__ constexpr void vec_copy(double destination[3], const double source[3]) {
  destination[0] = source[0];
  destination[1] = source[1];
  destination[2] = source[2];
}

__host__ __device__ constexpr double doppler_nucmf_on_nurf(const double dir_rf[3], const double vel_rf[3])
// Doppler factor
// arguments:
//   dir_rf: the rest frame direction (unit vector) of light propagation
//   vel_rf: velocity of the comoving frame relative to the rest frame
// returns: the ratio f = nu_cmf / nu_rf
{
  const double ndotv = dot(dir_rf, vel_rf);
  double dopplerfactor = 1. - (ndotv / CLIGHT);

  if (USE_RELATIVISTIC_DOPPLER_SHIFT) {
    const double betasq = dot(vel_rf, vel_rf) / CLIGHTSQUARED;
    dopplerfactor = dopplerfactor / sqrt(1 - betasq);
  }

  return dopplerfactor;
}

__host__ __device__ inline double doppler_packet_nucmf_on_nurf(const struct packet *const pkt_ptr) {
  double flow_velocity[3];  // homologous flow velocity
  get_velocity(pkt_ptr->pos, flow_velocity, pkt_ptr->prop_time);
  return doppler_nucmf_on_nurf(pkt_ptr->dir, flow_velocity);
}

__host__ __device__ constexpr void angle_ab(const double dir1[3], const double vel[3], double dir2[3])
// aberation of angles in special relativity
//   dir1: direction unit vector in frame1
//   vel: velocity of frame2 relative to frame1
//   dir2: direction vector in frame2
{
  const double vsqr = dot(vel, vel) / CLIGHTSQUARED;
  const double gamma_rel = 1. / (sqrt(1 - vsqr));

  const double ndotv = dot(dir1, vel);
  const double fact1 = gamma_rel * (1 - (ndotv / CLIGHT));
  const double fact2 = (gamma_rel - (gamma_rel * gamma_rel * ndotv / (gamma_rel + 1) / CLIGHT)) / CLIGHT;

  for (int d = 0; d < 3; d++) {
    dir2[d] = (dir1[d] - (vel[d] * fact2)) / fact1;
  }
}

__host__ __device__ constexpr void move_pkt(struct packet *pkt_ptr, const double distance)
/// Subroutine to move a packet along a straight line (specified by current
/// dir vector). The distance moved is in the rest frame.
{
  /// First update pos.

  pkt_ptr->pos[0] += (pkt_ptr->dir[0] * distance);
  pkt_ptr->pos[1] += (pkt_ptr->dir[1] * distance);
  pkt_ptr->pos[2] += (pkt_ptr->dir[2] * distance);

  /// During motion, rest frame energy and frequency are conserved.
  /// But need to update the co-moving ones.
  const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr);
  pkt_ptr->nu_cmf = pkt_ptr->nu_rf * dopplerfactor;
  pkt_ptr->e_cmf = pkt_ptr->e_rf * dopplerfactor;
}

#endif  // VECTORS_H
