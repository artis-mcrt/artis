#ifndef PACKET_INIT_H
#define PACKET_INIT_H

#include "types.h"

void packet_init(int middle_iteration, int my_rank);
double fni(const CELL *grid_ptr);
double f52fe(const CELL *grid_ptr);
double f48cr(const CELL *grid_ptr);
void write_packets(FILE *packets_file);
void read_packets(FILE *packets_file);

#endif //PACKET_INIT_H
