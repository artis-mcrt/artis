#include "sn3d.h"
#include "vectors.h"


/************************************************************/
/*Routine for aberation of angles in SR. Takes one direction and velocity
 as input and gives back another direction.*/
void angle_ab(const double *restrict const dir1, const double *restrict vel, double *dir2)
{
  const double vsqr = dot(vel,vel)/CLIGHTSQUARED;
  const double gamma_rel = 1./(sqrt(1 - vsqr));

  const double ndotv = dot(dir1,vel);
  const double fact1 = gamma_rel * (1 - (ndotv/CLIGHT));
  const double fact2 = (gamma_rel - (gamma_rel*gamma_rel*ndotv/(gamma_rel + 1)/CLIGHT))/CLIGHT;

  dir2[0] = (dir1[0] - (vel[0] * fact2))/fact1;
  dir2[1] = (dir1[1] - (vel[1] * fact2))/fact1;
  dir2[2] = (dir1[2] - (vel[2] * fact2))/fact1;
}


/************************************************************/
/*Routine for Doppler shift in SR. Takes one direction and velocity
 as input and gives back double.*/
double doppler (const double *restrict dir1, const double *restrict vel)
{
  //double vsqr = dot(vel,vel)/CLIGHTSQUARED;
  //double gamma_rel = 1./(sqrt(1 - vsqr));
  const double gamma_rel = 1.;
  const double ndotv = dot(dir1,vel);
  const double fact1 = gamma_rel * (1. - (ndotv/CLIGHT));

  #ifdef DEBUG_ON
  if (fabs(fact1 - 1) > 0.5)
  {
    printout("Dopper factor > 1.05?? Abort.\n");
    exit(0);
  }
  #endif

  return fact1;
}


/************************************************************/
/*Routine for scattering a direction through angle theta.*/
void scatter_dir(const double *dir_in, double cos_theta, double *dir_out)
{
  /*begin with setting the direction in coordinates where original direction
    is parallel to z-hat.*/

  const double zrand = gsl_rng_uniform(rng);
  const double phi = zrand * 2 * PI;

  const double sin_theta_sq = 1. - (cos_theta * cos_theta);
  const double sin_theta = sqrt(sin_theta_sq);
  const double zprime = cos_theta;
  const double xprime = sin_theta * cos(phi);
  const double yprime = sin_theta * sin(phi);

  /* Now need to derotate the coordinates back to real x,y,z. */
  /* Rotation matrix is determined by dir_in. */

  const double norm1 = 1./( sqrt( (dir_in[0]*dir_in[0]) + (dir_in[1]*dir_in[1]))) ;
  const double norm2 = 1./( sqrt( (dir_in[0]*dir_in[0]) + (dir_in[1]*dir_in[1]) + (dir_in[2]*dir_in[2])) );

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
