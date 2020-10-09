#ifndef VPKT_H
#define VPKT_H

#include <stdio.h>
#include "types.h"

double rot_angle(double *n1, double *n2, double *ref1, double *ref2);
void meridian(double *n, double *ref1, double *ref2);
void frame_transform(double *n_rf, double *Q, double *U, double *v, double *n_cmf);
void lorentz(double *e_rf, double *n_rf, double *v, double *e_cmf);

#ifdef VPKT_ON
void rlc_emiss_vpkt(PKT *pkt_ptr, double t_current, int bin, double *obs, int realtype);
int add_to_vspecpol(PKT *pkt_ptr, int bin, int ind, double t_arrive);
void init_vspecpol(void);
int write_vspecpol(FILE *specpol_file);
int read_vspecpol(FILE *specpol_file);
void init_vpkt_grid(void);
int add_to_vpkt_grid(PKT *dummy_ptr, double *vel, int bin_range, int bin, double *obs);
int write_vpkt_grid(FILE *vpkt_grid_file);
int read_vpkt_grid(FILE *vpkt_grid_file);
int check_tau(double *tau, double *tau_max);

extern int vgrid_flag;

extern int nvpkt;
extern int nvpkt_esc1; /* electron scattering event */
extern int nvpkt_esc2; /* kpkt deactivation */
extern int nvpkt_esc3; /* macroatom deactivation */

#endif

#endif //VPKT_H
