#ifndef VPKT_H
#define VPKT_H

#include "constants.h"
#include "packet.h"

namespace vpkt {

void read_vpktparameterfile();
void init(int nts, bool continued_from_saved);
void trace_vpkts(const Packet& pkt, enum packet_type type_before_rpkt);
void write_timestep(int nts, bool is_final);

void remove_temp_vpkt_file(int nts, int my_rank);

constexpr int VGRID_NY = 50;
constexpr int VGRID_NZ = 50;

// FREQUENCY
constexpr double VSPEC_NUMIN = CLIGHT / 10000 * 1e8;
constexpr double VSPEC_NUMAX = CLIGHT / 3500 * 1e8;
constexpr int VSPEC_NUBINS = 2500;

// TIME
constexpr double VSPEC_TIMEMIN = 3 * DAY;
constexpr double VSPEC_TIMEMAX = 8 * DAY;
constexpr int VSPEC_TIMEBINS = 5;

// number of virtual packets in a given timestep
inline int nvpkt_created{0};
inline int nvpkt_esc_from_rpkt{0};  // electron scattering event
inline int nvpkt_esc_from_kpkt{0};  // kpkt deactivation
inline int nvpkt_esc_from_macroatom{0};  // macroatom deactivation

inline double optical_depth_is_thick_vpkt;
}  // namespace vpkt

#endif  // VPKT_H
