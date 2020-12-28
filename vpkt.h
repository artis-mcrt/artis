#ifndef VPKT_H
#define VPKT_H

#include <cstdio>
#include "types.h"

double rot_angle(double *n1, double *n2, double *ref1, double *ref2);
void meridian(double *n, double *ref1, double *ref2);
void frame_transform(double *n_rf, double *Q, double *U, double *v, double *n_cmf);
void lorentz(double *e_rf, double *n_rf, double *v, double *e_cmf);

#ifdef VPKT_ON
void rlc_emiss_vpkt(PKT *pkt_ptr, double t_current, int bin, double *obs, int realtype);
void add_to_vspecpol(PKT *pkt_ptr, int bin, int ind, double t_arrive);
void init_vspecpol(void);
void write_vspecpol(FILE *specpol_file);
void read_vspecpol(FILE *specpol_file);
void init_vpkt_grid(void);
void add_to_vpkt_grid(PKT *dummy_ptr, double *vel, int bin_range, int bin, double *obs);
void write_vpkt_grid(FILE *vpkt_grid_file);
void read_vpkt_grid(FILE *vpkt_grid_file);
int check_tau(double *tau, double *tau_max);

// --------------------------------------------------------------------------------
// ---------------------------  VIRTUAL PACKETS -----------------------------------
// --------------------------------------------------------------------------------
#define MRANGE_GRID 5
#define NY_VGRID 50
#define NZ_VGRID 50

/* FREQUENCY */
/* dlognu = (log(numax) - log(numin)) / VMNUBINS ~ 3.9e-4 (10'000 over 1e14-5e15 Hz) */
#define numin_vspec (CLIGHT / 10000 * 1e8)
#define numax_vspec (CLIGHT / 3500 * 1e8)
#define VMNUBINS 2500

/* TIME */
/* dlogt = (log(globals::tmin) - log(globals::tmax)) / VMTBINS ~ 3.69e-2 (111 over 2-120 d) */
#define tmin_vspec (10 * DAY)
#define tmax_vspec (30 * DAY)
#define VMTBINS 30

/* Number of spectra for each observer (total + elements switched off) */
#define MSPECTRA 12
/* Number of observers */
#define MOBS 5
/* Total number of spectra */
#define MTOT (MSPECTRA * MOBS)
/* Total number of frequency ranges */
#define MRANGE 1

extern double Nrange_grid;
extern double tmin_grid;
extern double tmax_grid;
extern double nu_grid_min[MRANGE_GRID];
extern double nu_grid_max[MRANGE_GRID];
extern int vgrid_flag;
extern double dlogt_vspec;
extern double dlognu_vspec;

extern int nvpkt;
extern int nvpkt_esc1; /* electron scattering event */
extern int nvpkt_esc2; /* kpkt deactivation */
extern int nvpkt_esc3; /* macroatom deactivation */

extern int Nobs;
extern int Nspectra;
extern double nz_obs_vpkt[MOBS];
extern double phiobs[MOBS];
extern double tmin_vspec_input;
extern double tmax_vspec_input;
extern double Nrange;
extern double lmin_vspec_input[MRANGE];
extern double lmax_vspec_input[MRANGE];
extern double numin_vspec_input[MRANGE];
extern double numax_vspec_input[MRANGE];
extern double cell_is_optically_thick_vpkt;
extern double tau_max_vpkt;
extern double exclude[MSPECTRA];
extern double tau_vpkt[MSPECTRA];

#endif

#endif //VPKT_H
