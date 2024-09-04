#ifndef VPKT_H
#define VPKT_H

#include "constants.h"
#include "packet.h"

void read_parameterfile_vpkt();
void vpkt_init(int nts, int my_rank, bool continued_from_saved);
void vpkt_call_estimators(const Packet &pkt, enum packet_type type_before_rpkt);
void vpkt_write_timestep(int nts, int my_rank, bool is_final);

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
constexpr double VSPEC_TIMEMIN = 3 * DAY;
constexpr double VSPEC_TIMEMAX = 8 * DAY;
constexpr int VMTBINS = 5;

// number of virtual packets in a given timestep
inline int nvpkt{0};
inline int nvpkt_esc1{0};  // electron scattering event
inline int nvpkt_esc2{0};  // kpkt deactivation
inline int nvpkt_esc3{0};  // macroatom deactivation

inline double cell_is_optically_thick_vpkt;

#endif  // VPKT_H
