#ifndef VPKT_H
#define VPKT_H

#include <cstdio>
#include <span>

#include "artisoptions.h"

double rot_angle(std::span<double, 3> n1, std::span<double, 3> n2, std::span<double, 3> ref1,
                 std::span<double, 3> ref2);
void meridian(std::span<const double, 3> n, std::span<double, 3> ref1, std::span<double, 3> ref2);
void frame_transform(std::span<const double, 3> n_rf, double *Q, double *U, std::span<const double, 3> v,
                     std::span<double, 3> n_cmf);
void lorentz(std::span<const double, 3> e_rf, std::span<const double, 3> n_rf, std::span<const double, 3> v,
             std::span<double, 3> e_cmf);

void rlc_emiss_vpkt(const struct packet *const pkt_ptr, double t_current, int bin, std::span<double, 3> obs,
                    int realtype);
void add_to_vspecpol(const struct packet *const pkt_ptr, int bin, int ind, double t_arrive);
void init_vspecpol();
void read_parameterfile_vpkt();
void write_vspecpol(FILE *specpol_file);
void read_vspecpol(int my_rank, int nts);
void init_vpkt_grid();
void add_to_vpkt_grid(const struct packet *const dummy_ptr, std::span<const double, 3> vel, int bin_range, int bin,
                      std::span<const double, 3> obs);
void write_vpkt_grid(FILE *vpkt_grid_file);
void read_vpkt_grid(FILE *vpkt_grid_file);
int check_tau(const double *tau, const double *tau_max);
int vpkt_call_estimators(struct packet *pkt_ptr, int realtype);

// --------------------------------------------------------------------------------
// ---------------------------  VIRTUAL PACKETS -----------------------------------
// --------------------------------------------------------------------------------
constexpr int VGRID_NY = 50;
constexpr int VGRID_NZ = 50;

// FREQUENCY
// dlognu = (log(numax) - log(numin)) / VMNUBINS ~ 3.9e-4 (10'000 over 1e14-5e15 Hz)
constexpr double VSPEC_NUMIN = CLIGHT / 10000 * 1e8;
constexpr double VSPEC_NUMAX = CLIGHT / 3500 * 1e8;
constexpr int VMNUBINS = 2500;

// TIME
// dlogt = (log(globals::tmin) - log(globals::tmax)) / VMTBINS ~ 3.69e-2 (111 over 2-120 d)
constexpr double VSPEC_TIMEMIN = 10 * DAY;
constexpr double VSPEC_TIMEMAX = 30 * DAY;
constexpr int VMTBINS = 30;

extern bool vgrid_on;

extern int nvpkt;
extern int nvpkt_esc1;  // electron scattering event
extern int nvpkt_esc2;  // kpkt deactivation
extern int nvpkt_esc3;  // macroatom deactivation

extern double cell_is_optically_thick_vpkt;

#endif  // VPKT_H
