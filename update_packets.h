#ifndef UPDATE_PACKETS_H
#define UPDATE_PACKETS_H

#include "sn3d.h"
#include "update_grid.h"

void update_packets_gpu(int nts, PKT *pkt);
void update_packets(int nts, PKT *pkt);

#endif //UPDATE_PACKETS_H
