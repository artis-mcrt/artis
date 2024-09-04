#ifndef UPDATE_PACKETS_H
#define UPDATE_PACKETS_H

#include <span>

#include "packet.h"

void update_packets(int my_rank, int nts, std::span<Packet> packets);

#endif  // UPDATE_PACKETS_H
