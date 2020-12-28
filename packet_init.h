#ifndef PACKET_INIT_H
#define PACKET_INIT_H

#include <cstdio>
#include "types.h"

void packet_init(int middle_iteration, int my_rank, PKT *pkt);
void write_packets(char filename[], PKT *pkt);
void read_packets(char filename[], PKT *pkt);
void read_temp_packetsfile(const int timestep, const int my_rank, PKT *const pkt);

#endif //PACKET_INIT_H
