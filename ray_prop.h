#ifndef RAY_PROP_H
#define RAY_PROP_H

#include "types.h"

double do_gamma_ray(RAY *ray_ptr, double t1, double t2);
int get_nul(double freq);
double get_gam_freq(const LIST *restrict line_list, int n);

#endif //RAY_PROP_H
