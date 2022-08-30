#ifndef VPKT_H
#define VPKT_H

#include <cstdio>

#include "artisoptions.h"
#include "cuda.h"
#include "constants.h"

double rot_angle(double *n1, double *n2, double *ref1, double *ref2);
void meridian(const double *n, double *ref1, double *ref2);
void frame_transform(double *n_rf, double *Q, double *U, double *v, double *n_cmf);
void lorentz(double *e_rf, double *n_rf, double *v, double *e_cmf);

void rlc_emiss_vpkt(struct packet *pkt_ptr, double t_current, int bin, double *obs, int realtype);
void add_to_vspecpol(struct packet *pkt_ptr, int bin, int ind, double t_arrive);
void init_vspecpol(void);
void read_parameterfile_vpkt(void);
void write_vspecpol(FILE *specpol_file);
void read_vspecpol(int my_rank, int nts);
void init_vpkt_grid(void);
void add_to_vpkt_grid(struct packet *dummy_ptr, double *vel, int bin_range, int bin, double *obs);
void write_vpkt_grid(FILE *vpkt_grid_file);
void read_vpkt_grid(FILE *vpkt_grid_file);
int check_tau(double *tau, double *tau_max);
__host__ __device__ int vpkt_call_estimators(struct packet *pkt_ptr, double t_current, int realtype);

// --------------------------------------------------------------------------------
// ---------------------------  VIRTUAL PACKETS -----------------------------------
// --------------------------------------------------------------------------------
#define MRANGE_GRID 5
#define NY_VGRID 50
#define NZ_VGRID 50

// FREQUENCY
// dlognu = (log(numax) - log(numin)) / VMNUBINS ~ 3.9e-4 (10'000 over 1e14-5e15 Hz)
#define numin_vspec (CLIGHT / 10000 * 1e8)
#define numax_vspec (CLIGHT / 3500 * 1e8)
#define VMNUBINS 2500

// TIME
// dlogt = (log(globals::tmin) - log(globals::tmax)) / VMTBINS ~ 3.69e-2 (111 over 2-120 d)
#define tmin_vspec (10 * DAY)
#define tmax_vspec (30 * DAY)
#define VMTBINS 30

// Total number of frequency ranges
#define MRANGE 1

extern int vgrid_flag;

extern int nvpkt;
extern int nvpkt_esc1;  // electron scattering event
extern int nvpkt_esc2;  // kpkt deactivation
extern int nvpkt_esc3;  // macroatom deactivation

extern double cell_is_optically_thick_vpkt;

#endif  // VPKT_H
