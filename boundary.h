#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "types.h"
#include "cuda.h"

__host__ __device__ double boundary_cross(PKT *pkt_ptr, double tstart, int *snext);
__host__ __device__ void change_cell(PKT *pkt_ptr, int snext, double t_current);

#endif //BOUNDARY_H
