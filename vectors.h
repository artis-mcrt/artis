#ifndef VECTORS_H
#define VECTORS_H

#include <math.h>
#include <gsl/gsl_blas.h>
#include "types.h"

__host__ __device__ void angle_ab(const double dir1[3], const double vel[3], double dir2[3]);
__host__ __device__ double doppler(const double dir1[3], const double vel[3]);
__host__ __device__ void scatter_dir(const double dir_in[3], double cos_theta, double dir_out[3]);
__host__ __device__ void get_rand_isotropic_unitvec(double vecout[3]);
__host__ __device__ void move_pkt(PKT *pkt_ptr, double distance, double time);

// #define vec_len(x)   (cblas_dnrm2(3, x, 1))

__host__ __device__
inline double vec_len(const double x[3])
// return the the magnitude of a vector
{
#ifdef __CUDA_ARCH__
  return sqrt((x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]));
#else
  return cblas_dnrm2(3, x, 1);
#endif
}


__host__ __device__
inline void vec_norm(const double vec_in[3], double vec_out[3])
// Routine for normalizing a vector.
{
  const double magnitude = vec_len(vec_in);

  vec_out[0] = vec_in[0] / magnitude;
  vec_out[1] = vec_in[1] / magnitude;
  vec_out[2] = vec_in[2] / magnitude;

  // vec_copy(vec_out, vec_in);
  // cblas_dscal(3, 1 / magnitude, vec_out, 1)
}


__host__ __device__
inline double dot(const double x[3], const double y[3])
// Routine for taking dot product.
{

#ifdef __CUDA_ARCH__
  return (x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]);
#else
  return cblas_ddot(3, x, 1, y, 1);
#endif
}


__host__ __device__
inline void get_velocity(const double x[3], double y[3], const double t)
// Routine for getting velocity vector of the flow at a position with homologous expansion.
{
  y[0] = x[0] / t;
  y[1] = x[1] / t;
  y[2] = x[2] / t;
}


__host__ __device__
inline void cross_prod(const double vec1[3], const double vec2[3], double vecout[3])
{
  vecout[0] = (vec1[1] * vec2[2]) - (vec2[1] * vec1[2]);
  vecout[1] = (vec1[2] * vec2[0]) - (vec2[2] * vec1[0]);
  vecout[2] = (vec1[0] * vec2[1]) - (vec2[0] * vec1[1]);
}


__host__ __device__
inline void
vec_scale(double vec[3], const double scalefactor)
{
#ifdef __CUDA_ARCH__
  for (int d = 0; d < 3; d++)
  {
    vec[d] *= scalefactor;
  }
#else
  cblas_dscal(3, scalefactor, vec, 1);
#endif
}


__host__ __device__
inline void vec_copy(double destination[3], const double source[3])
{
#ifdef __CUDA_ARCH__
  for (int d = 0; d < 3; d++)
  {
    destination[d] = source[d];
  }
#else
  cblas_dcopy(3, source, 1, destination, 1);
#endif
}


__host__ __device__
inline double
doppler_packetpos(const PKT *const pkt_ptr, const double t)
{
  double vel_vec[3];
  get_velocity(pkt_ptr->pos, vel_vec, t); // homologous flow velocity
  const double dopplerfactor = doppler(pkt_ptr->dir, vel_vec);
  return dopplerfactor;
}

#endif //VECTORS_H
