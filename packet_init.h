#ifndef PACKET_INIT_H
#define PACKET_INIT_H

#include <stdio.h>
#include "types.h"

void packet_init(int middle_iteration, int my_rank, PKT *pkt);
void write_packets(char filename[], PKT *pkt);
void read_packets(char filename[], PKT *pkt);

#endif //PACKET_INIT_H
