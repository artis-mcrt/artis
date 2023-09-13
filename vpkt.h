#ifndef VPKT_H
#define VPKT_H

#include <cstdio>

#include "artisoptions.h"

double rot_angle(double *n1, double *n2, double *ref1, double *ref2);
void meridian(const double *n, double *ref1, double *ref2);
void frame_transform(const double (&n_rf)[3], double *Q, double *U, const double (&v)[3], double (&n_cmf)[3]);
void lorentz(const double *e_rf, const double *n_rf, const double *v, double *e_cmf);

void rlc_emiss_vpkt(struct packet *pkt_ptr, double t_current, int bin, double (&obs)[3], int realtype);
void add_to_vspecpol(struct packet *pkt_ptr, int bin, int ind, double t_arrive);
void init_vspecpol();
void read_parameterfile_vpkt();
void write_vspecpol(FILE *specpol_file);
void read_vspecpol(int my_rank, int nts);
void init_vpkt_grid();
void add_to_vpkt_grid(struct packet *dummy_ptr, const double *vel, int bin_range, int bin, const double *obs);
void write_vpkt_grid(FILE *vpkt_grid_file);
void read_vpkt_grid(FILE *vpkt_grid_file);
int check_tau(const double *tau, const double *tau_max);
int vpkt_call_estimators(struct packet *pkt_ptr, double t_current, int realtype);

// --------------------------------------------------------------------------------
// ---------------------------  VIRTUAL PACKETS -----------------------------------
// --------------------------------------------------------------------------------
constexpr int MRANGE_GRID = 5;
constexpr int NY_VGRID = 50;
constexpr int NZ_VGRID = 50;

// FREQUENCY
// dlognu = (log(numax) - log(numin)) / VMNUBINS ~ 3.9e-4 (10'000 over 1e14-5e15 Hz)
constexpr double numin_vspec = (CLIGHT / 10000 * 1e8);
constexpr double numax_vspec = (CLIGHT / 3500 * 1e8);
constexpr int VMNUBINS = 2500;

// TIME
// dlogt = (log(globals::tmin) - log(globals::tmax)) / VMTBINS ~ 3.69e-2 (111 over 2-120 d)
constexpr double tmin_vspec = (10 * DAY);
constexpr double tmax_vspec = (30 * DAY);
constexpr int VMTBINS = 30;

// Total number of frequency ranges
constexpr int MRANGE = 2;

extern int vgrid_flag;

extern int nvpkt;
extern int nvpkt_esc1;  // electron scattering event
extern int nvpkt_esc2;  // kpkt deactivation
extern int nvpkt_esc3;  // macroatom deactivation

extern double cell_is_optically_thick_vpkt;

#endif  // VPKT_H
