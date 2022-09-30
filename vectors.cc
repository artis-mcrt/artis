#include "vectors.h"

// #include <gsl/gsl_randist.h>

#include <cmath>

#include "artisoptions.h"
#include "sn3d.h"

__host__ __device__ void angle_ab(const double dir1[3], const double vel[3], double dir2[3])
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
}

__host__ __device__ void scatter_dir(const double dir_in[3], const double cos_theta, double dir_out[3])
// Routine for scattering a direction through angle theta.
{
  // begin with setting the direction in coordinates where original direction
  // is parallel to z-hat.

  const double zrand = gsl_rng_uniform(rng);
  const double phi = zrand * 2 * PI;

  const double sin_theta_sq = 1. - (cos_theta * cos_theta);
  const double sin_theta = std::sqrt(sin_theta_sq);
  const double zprime = cos_theta;
  const double xprime = sin_theta * cos(phi);
  const double yprime = sin_theta * sin(phi);

  // Now need to derotate the coordinates back to real x,y,z.
  // Rotation matrix is determined by dir_in.

  const double norm1 = 1. / std::sqrt((dir_in[0] * dir_in[0]) + (dir_in[1] * dir_in[1]));
  const double norm2 = 1. / std::sqrt((dir_in[0] * dir_in[0]) + (dir_in[1] * dir_in[1]) + (dir_in[2] * dir_in[2]));

  const double r11 = dir_in[1] * norm1;
  const double r12 = -1 * dir_in[0] * norm1;
  const double r13 = 0.0;
  const double r21 = dir_in[0] * dir_in[2] * norm1 * norm2;
  const double r22 = dir_in[1] * dir_in[2] * norm1 * norm2;
  const double r23 = -1 * norm2 / norm1;
  const double r31 = dir_in[0] * norm2;
  const double r32 = dir_in[1] * norm2;
  const double r33 = dir_in[2] * norm2;

  dir_out[0] = (r11 * xprime) + (r21 * yprime) + (r31 * zprime);
  dir_out[1] = (r12 * xprime) + (r22 * yprime) + (r32 * zprime);
  dir_out[2] = (r13 * xprime) + (r23 * yprime) + (r33 * zprime);
}

__host__ __device__ void get_rand_isotropic_unitvec(double vecout[3])
// Assuming isotropic distribution, get a random direction vector
{
  // alternatively, use GSL's functions:
  // gsl_ran_dir_3d(rng, &vecout[0], &vecout[1], &vecout[2]);
  // or
  // gsl_ran_dir_nd(rng, 3, vecout);
  // but check validity first

  const double zrand = gsl_rng_uniform(rng);
  const double zrand2 = gsl_rng_uniform(rng);

  const double mu = -1 + (2. * zrand);
  const double phi = zrand2 * 2 * PI;
  const double sintheta = std::sqrt(1. - (mu * mu));

  vecout[0] = sintheta * std::cos(phi);
  vecout[1] = sintheta * std::sin(phi);
  vecout[2] = mu;
}

__host__ __device__ void move_pkt(struct packet *pkt_ptr, const double distance)
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

__host__ __device__ void move_pkt_withtime(struct packet *pkt_ptr, const double distance)
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