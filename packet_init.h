#ifndef PACKET_INIT_H
#define PACKET_INIT_H

#include "types.h"

void packet_init(int middle_iteration, int my_rank);
void write_packets(FILE *restrict packets_file);
void read_packets(FILE *restrict packets_file);

#endif //PACKET_INIT_H
