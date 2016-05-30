#ifndef VECTORS_H
#define VECTORS_H

#include <math.h>

void angle_ab(const double *restrict dir1, const double *restrict vel, double *dir2);
double doppler(const double *restrict dir1, const double *restrict vel);
void scatter_dir(const double *restrictdir_in, double cos_theta, double *dir_out);

/*Routine for getting the magnitude of a vector.*/
static inline
double vec_len(const double x[3])
{
  return sqrt((x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]));
}


/*Routine for normalizing a vector.*/
static inline
void vec_norm(const double vec_in[3], double vec_out[3])
{
  const double magnitude = vec_len(vec_in);

  vec_out[0] = vec_in[0] / magnitude;
  vec_out[1] = vec_in[1] / magnitude;
  vec_out[2] = vec_in[2] / magnitude;
}


/*Routine for taking dot product.*/
static inline
double dot(const double *const restrict x, const double *const restrict y)
{
  return (x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]);
}


/*Routine for getting velocity vector of the flow at a position.*/
static inline
void get_velocity(const double *const restrict x, double *restrict y, const double t)
{
  /* For homologous expansion. */

  y[0] = x[0] / t;
  y[1] = x[1] / t;
  y[2] = x[2] / t;
}


static inline
void cross_prod(const double vec1[3], const double vec2[3], double vecout[3])
{
  const double vec1x = vec1[0];
  const double vec1y = vec1[1];
  const double vec1z = vec1[2];
  const double vec2x = vec2[0];
  const double vec2y = vec2[1];
  const double vec2z = vec2[2];
  vecout[0] = (vec1y*vec2z) - (vec2y*vec1z);
  vecout[1] = (vec1z*vec2x) - (vec2z*vec1x);
  vecout[2] = (vec1x*vec2y) - (vec2x*vec1y);
}

#endif //VECTORS_H
