#include "vectors.h"

#include <cmath>

#include "sn3d.h"

void scatter_dir(std::span<const double, 3> dir_in, const double cos_theta, std::span<double, 3> dir_out)
// Routine for scattering a direction through angle theta.
{
  // begin with setting the direction in coordinates where original direction
  // is parallel to z-hat.

  const double zrand = rng_uniform();
  const double phi = zrand * 2 * PI;

  const double sin_theta_sq = 1. - (cos_theta * cos_theta);
  const double sin_theta = std::sqrt(sin_theta_sq);
  const double zprime = cos_theta;
  const double xprime = sin_theta * cos(phi);
  const double yprime = sin_theta * sin(phi);

  // Now need to derotate the coordinates back to real x,y,z.
  // Rotation matrix is determined by dir_in.

  const double norm1 = 1. / std::sqrt((dir_in[0] * dir_in[0]) + (dir_in[1] * dir_in[1]));
  const double norm2 = 1. / vec_len(dir_in);

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

  assert_testmodeonly(std::fabs(vec_len(dir_out) - 1.) < 1e-10);
}

void get_rand_isotropic_unitvec(std::span<double, 3> vecout)
// Assuming isotropic distribution, get a random direction vector
{
  // alternatively, use GSL's functions:
  // gsl_ran_dir_3d(rng, &vecout[0], &vecout[1], &vecout[2]);
  // or
  // gsl_ran_dir_nd(rng, 3, vecout);
  // but check validity first

  const double zrand = rng_uniform();
  const double zrand2 = rng_uniform();

  const double mu = -1 + (2. * zrand);
  const double phi = zrand2 * 2 * PI;
  const double sintheta = std::sqrt(1. - (mu * mu));

  vecout[0] = sintheta * std::cos(phi);
  vecout[1] = sintheta * std::sin(phi);
  vecout[2] = mu;

  assert_testmodeonly(std::fabs(vec_len(vecout) - 1.) < 1e-10);
}
