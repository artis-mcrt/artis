#ifndef KPKT_H
#define KPKT_H

#include "sn3d.h"
#include "thermalbalance.h"
#include "types.h"

void setup_coolinglist(void);
__host__ __device__ void calculate_cooling_rates(int modelgridindex, heatingcoolingrates_t *heatingcoolingrates);
__host__ __device__ double do_kpkt_bb(PKT *pkt_ptr);
__host__ __device__ double do_kpkt(PKT *pkt_ptr, double t2, int nts);
__host__ __device__ int get_coolinglistoffset(int element, int ion);

#endif  // KPKT_H
