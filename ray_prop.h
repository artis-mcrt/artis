#ifndef RAY_PROP_H
#define RAY_PROP_H

int ray_prop(RAY *ray_ptr, double t1, double t2, int nts);
double do_gamma_ray(RAY *ray_ptr, double t1, double t2);
double boundary_cross_ray(RAY *ray_ptr, double tstart, int *snext);
int change_cell_ray(RAY *ray_ptr, int snext, int *end_packet, double t_current);
int copy_ray (RAY *ray1, RAY *ray2);
int move_ray(RAY *ray_ptr, double dist, double time);
int move_one_ray(RAY *ray_ptr, int nray, double dist, double *single_pos, double single_t);
int get_nul(double freq);
double get_gam_freq(LIST *line_list, int n);

#endif //RAY_PROP_H
