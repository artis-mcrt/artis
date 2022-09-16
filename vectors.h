#ifndef VECTORS_H
#define VECTORS_H

#include <cmath>

#include "cuda.h"
#include "packet.h"

__host__ __device__ void angle_ab(const double dir1[3], const double vel[3], double dir2[3]);
__host__ __device__ double doppler_nucmf_on_nurf(const double dir_rf[3], const double vel_rf[3]);
__host__ __device__ void scatter_dir(const double dir_in[3], double cos_theta, double dir_out[3]);
__host__ __device__ void get_rand_isotropic_unitvec(double vecout[3]);
__host__ __device__ void move_pkt(struct packet *pkt_ptr, double distance);
__host__ __device__ void move_pkt_withtime(struct packet *pkt_ptr, double distance);

__host__ __device__ static inline double vec_len(const double x[3])
// return the the magnitude of a vector
{
  return std::sqrt((x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]));
}

__host__ __device__ inline void vec_norm(const double vec_in[3], double vec_out[3])
// normalizing a copy of vec_in and save it to vec_out
{
  const double magnitude = vec_len(vec_in);

  vec_out[0] = vec_in[0] / magnitude;
  vec_out[1] = vec_in[1] / magnitude;
  vec_out[2] = vec_in[2] / magnitude;
}

__host__ __device__ constexpr double dot(const double x[3], const double y[3])
// vector dot product
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

__host__ __device__ inline double doppler_packet_nucmf_on_nurf(const struct packet *const pkt_ptr) {
  double flow_velocity[3];  // homologous flow velocity
  get_velocity(pkt_ptr->pos, flow_velocity, pkt_ptr->prop_time);
  return doppler_nucmf_on_nurf(pkt_ptr->dir, flow_velocity);
}

#endif  // VECTORS_H
