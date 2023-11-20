#include "vectors.h"

#include <cmath>
#include <span>

#include "constants.h"
#include "sn3d.h"

auto get_rand_isotropic_unitvec() -> std::array<double, 3>
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

  std::array<double, 3> vecout = {sintheta * std::cos(phi), sintheta * std::sin(phi), mu};

  assert_testmodeonly(std::fabs(vec_len(vecout) - 1.) < 1e-10);
  return vecout;
}
