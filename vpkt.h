#ifndef VPKT_H
#define VPKT_H

int rlc_emiss_vpkt(PKT *pkt_ptr, double t_current, int bin, double *obs, int realtype);
double rot_angle(double *n1, double *n2, double *ref1, double *ref2);
void meridian(double *n, double *ref1, double *ref2);
void frame_transform(double *n_rf, double *Q, double *U, double *v, double *n_cmf);
void lorentz(double *e_rf, double *n_rf, double *v, double *e_cmf);
int add_to_vspecpol(PKT *pkt_ptr, int bin, double t_arrive);
void init_vspecpol(void);
int write_vspecpol(FILE *specpol_file);
int read_vspecpol(FILE *specpol_file);
void init_vpkt_grid(void);
int add_to_vpkt_grid(PKT *dummy_ptr, double *vel, int bin_range, int bin, double *obs);
int write_vpkt_grid(FILE *vpkt_grid_file);
int read_vpkt_grid(FILE *vpkt_grid_file);

#endif //VPKT_H
