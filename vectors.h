#ifndef VECTORS_H
#define VECTORS_H

double vec_len(double x[3]);
void vec_norm(double x[3], double z[3]);
double dot(const double *x, const double *y);
int get_velocity(const double *x, double *y, double t);
int angle_ab(const double *dir1, const double *vel, double *dir2);
double doppler(const double *dir1, const double *vel);
int scatter_dir(const double *dir_in, double cos_theta, double *dir_out);
int cross_prod(double v1[3], double v2[3], double v3[3]);

#endif //VECTORS_H
