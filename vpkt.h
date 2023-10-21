#ifndef VPKT_H
#define VPKT_H

#include <span>

#include "artisoptions.h"
#include "packet.h"

double rot_angle(std::span<double, 3> n1, std::span<double, 3> n2, std::span<double, 3> ref1,
                 std::span<double, 3> ref2);
void meridian(std::span<const double, 3> n, std::span<double, 3> ref1, std::span<double, 3> ref2);
void frame_transform(std::span<const double, 3> n_rf, double *Q, double *U, std::span<const double, 3> v,
                     std::span<double, 3> n_cmf);

void read_parameterfile_vpkt();
void vpkt_init(int nts, int my_rank, int tid, bool continued_from_saved);
void vpkt_call_estimators(struct packet *pkt_ptr, const enum packet_type);
void vpkt_write_timestep(int nts, int my_rank, int tid, bool is_final);

void vpkt_remove_temp_file(int nts, int my_rank);

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

extern int nvpkt;
extern int nvpkt_esc1;  // electron scattering event
extern int nvpkt_esc2;  // kpkt deactivation
extern int nvpkt_esc3;  // macroatom deactivation

extern double cell_is_optically_thick_vpkt;

#endif  // VPKT_H
