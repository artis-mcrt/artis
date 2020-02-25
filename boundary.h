#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "types.h"

__host__ __device__ double boundary_cross(PKT *pkt_ptr, double tstart, int *snext);
__host__ __device__ void change_cell(PKT *pkt_ptr, int snext, bool *end_packet, double t_current, int tid);
void change_cell_vpkt(PKT *pkt_ptr, int snext, bool *end_packet, double t_current);

#endif //BOUNDARY_H
