#ifndef VECTORS_H
#define VECTORS_H

#include <math.h>

int angle_ab(const double *dir1, const double *vel, double *dir2);
double doppler(const double *dir1, const double *vel);
int scatter_dir(const double *dir_in, double cos_theta, double *dir_out);

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
  double magnitude = vec_len(vec_in);

  vec_out[0] = vec_in[0] / magnitude;
  vec_out[1] = vec_in[1] / magnitude;
  vec_out[2] = vec_in[2] / magnitude;
}


/*Routine for taking dot product.*/
static inline
double dot(const double *x, const double *y)
{
  return (x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]);
}


/*Routine for getting velocity vector of the flow at a position.*/
static inline
void get_velocity(const double *x, double *y, const double t)
{
  /* For homologous expansion. */

  y[0] = x[0] / t;
  y[1] = x[1] / t;
  y[2] = x[2] / t;
}


static inline
void cross_prod(const double vec1[3], const double vec2[3], double vecout[3])
{
  vecout[0] = (vec1[1]*vec2[2]) - (vec2[1]*vec1[2]);
  vecout[1] = (vec1[2]*vec2[0]) - (vec2[2]*vec1[0]);
  vecout[2] = (vec1[0]*vec2[1]) - (vec2[0]*vec1[1]);
}

#endif //VECTORS_H
