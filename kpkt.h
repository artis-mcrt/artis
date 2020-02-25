#ifndef KPKT_H
#define KPKT_H

#include "sn3d.h"
#include "types.h"

void calculate_cooling_rates(int modelgridindex, heatingcoolingrates_t *heatingcoolingrates);
__host__ __device__ double do_kpkt_bb(PKT *pkt_ptr, double t1, int tid);
__host__ __device__ double do_kpkt(PKT *pkt_ptr, double t1, double t2, int nts, int tid);

inline __host__ __device__ int get_coolinglistoffset(int element, int ion)
{
  return elements[element].ions[ion].coolingoffset;
}

#endif //KPKT_H
