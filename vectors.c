#include "sn3d.h"
/***********************************************************/

/*Routine for getting the magnitude of a vector.*/

double
vec_len (x)
     double x[3];
{
  double y;

  y = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);

  y = sqrt(y);

  if (y < 0)
    {
      printout("Error: vec_len. Abort.\n");
      exit(0);
    }

  return(y);
}

/************************************************************/
/*Routine for normalizing a vector.*/

void
vec_norm (double x[3], double z[3])
{
  double y;

  y = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);

  y = sqrt(y);

  if (y < 0)
    {
      printout("Error: vec_len. Abort.\n");
      exit(0);
    }

  z[0]=x[0]/y;
  z[1]=x[1]/y;
  z[2]=x[2]/y;
}

/************************************************************/
/*Routine for taking dot product.*/

double
dot (x,y)
     double *x;
     double *y;
{
  double result;

  result = (x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]);

  return(result);
}

/************************************************************/
/*Routine for getting velocity vector of the flow at a position.*/

int
get_velocity (x,y,t)
     double *x, *y;
     double t;
{
  /* For homologous expansion. */

  y[0] = x[0] / t;
  y[1] = x[1] / t;
  y[2] = x[2] / t;

  return(0);
}
/************************************************************/

/*Routine for aberation of angles in SR. Takes one direction and velocity
 as input and gives back another direction.*/

int
angle_ab (dir1,vel,dir2)
     double *dir1,*vel,*dir2;
{
  double gamma_rel;
  double dot();
  void vec_norm();
  double ndotv, fact2, fact1;
  double vsqr;

  vsqr = dot(vel,vel)/CLIGHT2;
  gamma_rel = 1./(sqrt(1 - vsqr));
  
  ndotv = dot(dir1,vel);
  fact1 = gamma_rel * (1 - (ndotv/CLIGHT));
  fact2 = (gamma_rel - (gamma_rel*gamma_rel*ndotv/(gamma_rel + 1)/CLIGHT))/CLIGHT;

  dir2[0] = (dir1[0] - (vel[0] * fact2))/fact1;
  dir2[1] = (dir1[1] - (vel[1] * fact2))/fact1;
  dir2[2] = (dir1[2] - (vel[2] * fact2))/fact1;
  
  //dir2[0]=dir1[0];
  //dir2[1]=dir1[1];
  //dir2[2]=dir1[2];
    
  vec_norm(dir2,dir2);
    
  return(0);
}


int
angle_ab2 (dir1,vel,dir2)
double *dir1,*vel,*dir2;
{
    double gamma_rel;
    double dot();
    double ndotv, fact2, fact1;
    double vsqr;
    
    vsqr = dot(vel,vel)/CLIGHT2;
    gamma_rel = 1./(sqrt(1 - vsqr));
    
    ndotv = dot(dir1,vel);
    fact1 = gamma_rel * (1 - (ndotv/CLIGHT));
    fact2 = (gamma_rel - (gamma_rel*gamma_rel*ndotv/(gamma_rel + 1)/CLIGHT))/CLIGHT;
  
    dir2[0]=dir1[0];
    dir2[1]=dir1[1];
    dir2[2]=dir1[2];
    
    return(0);
}


/************************************************************/

/*Routine for Doppler shift in SR. Takes one direction and velocity
 as input and gives back double.*/

double
doppler (dir1,vel)
     double *dir1,*vel;
{
  double gamma_rel;
  double dot();
  double ndotv, fact1;
  double vsqr;

  //vsqr = dot(vel,vel)/CLIGHT2;
  //gamma_rel = 1./(sqrt(1 - vsqr));
  gamma_rel = 1.;

  ndotv = dot(dir1,vel);
 
  fact1 = gamma_rel * (1. - (ndotv/CLIGHT));

  if (fabs(fact1-1) > 0.5)
    {
      printout("Dopper factor > 1.05?? Abort.\n");
      exit(0);
    }

  return(fact1);
}
/************************************************************/

/*Routine for scattering a direction through angle theta.*/

int
scatter_dir(dir_in, cos_theta, dir_out)
     double *dir_in,*dir_out;
     double cos_theta;
{
  double xprime,yprime,zprime;
  double zrand, phi;
  double norm1, norm2;
  double r11,r12,r13,r21,r22,r23,r31,r32,r33;
  double sin_theta;
  
  /*begin with setting the direction in coordinates where original direction 
    is parallel to z-hat.*/

  zrand = gsl_rng_uniform(rng);
  phi = zrand * 2 * PI;

  sin_theta = 1. - (cos_theta*cos_theta);
  sin_theta = sqrt(sin_theta);
  zprime = cos_theta;
  xprime = sin_theta * cos(phi);
  yprime = sin_theta * sin(phi);

  /* Now need to derotate the coordinates back to real x,y,z. */
  /* Rotation matrix is determined by dir_in. */

  norm1=1./(sqrt( (dir_in[0]*dir_in[0]) + (dir_in[1]*dir_in[1])));
  norm2=1./(sqrt( (dir_in[0]*dir_in[0]) + (dir_in[1]*dir_in[1]) + (dir_in[2]*dir_in[2])));

  r11 = dir_in[1] * norm1;
  r12 = -1 * dir_in[0] * norm1;
  r13 = 0.0;
  r21 = dir_in[0] * dir_in[2] * norm1 * norm2;
  r22 = dir_in[1] * dir_in[2] * norm1 * norm2;
  r23 = -1 * norm2 / norm1;
  r31 = dir_in[0] * norm2;
  r32 = dir_in[1] * norm2;
  r33 = dir_in[2] * norm2;
  
  dir_out[0] = (r11 * xprime) + (r21 * yprime) + (r31 * zprime);
  dir_out[1] = (r12 * xprime) + (r22 * yprime) + (r32 * zprime);
  dir_out[2] = (r13 * xprime) + (r23 * yprime) + (r33 * zprime);
  
  return(0);


}


/**********************************************************************/
int cross_prod (v1, v2, v3)
     double v1[3], v2[3], v3[3];
{
  v3[0] = (v1[1]*v2[2]) - (v2[1]*v1[2]);
  v3[1] = (v1[2]*v2[0]) - (v2[2]*v1[0]);
  v3[2] = (v1[0]*v2[1]) - (v2[0]*v1[1]);

  return(0);
}

