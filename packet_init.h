#ifndef PACKET_INIT_H
#define PACKET_INIT_H

#include "types.h"

void packet_init(int middle_iteration, int my_rank);
double fni(const CELL *restrict grid_ptr);
double f52fe(const CELL *restrict grid_ptr);
double f48cr(const CELL *restrict grid_ptr);
void write_packets(FILE *restrict packets_file);
void read_packets(FILE *restrict packets_file);

#endif //PACKET_INIT_H
