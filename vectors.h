#ifndef VECTORS_H
#define VECTORS_H

double vec_len(double x[3]);
void vec_norm(double x[3], double z[3]);
double dot(double *x, double *y);
int get_velocity(double *x, double *y, double t);
int angle_ab(double *dir1, double *vel, double *dir2);
double doppler(double *dir1, double *vel);
int scatter_dir(double *dir_in, double cos_theta, double *dir_out);
int cross_prod(double v1[3], double v2[3], double v3[3]);

#endif //VECTORS_H
