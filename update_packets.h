// Declarations for the per-timestep packet update (update_packets.cc).

#ifndef UPDATE_PACKETS_H
#define UPDATE_PACKETS_H

#include <span>

#include "packet.h"

void update_packets(int nts, std::span<Packet> packets);

#endif  // UPDATE_PACKETS_H
