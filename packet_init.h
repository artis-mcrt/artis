#ifndef PACKET_INIT_H
#define PACKET_INIT_H

#include <stdio.h>

void packet_init(int middle_iteration, int my_rank, PKT *pkt);
void write_packets(FILE *restrict packets_file, PKT *pkt);
void read_packets(FILE *restrict packets_file, PKT *pkt);

#endif //PACKET_INIT_H
