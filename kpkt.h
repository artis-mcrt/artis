#ifndef KPKT_H
#define KPKT_H

#include "cuda.h"
#include "globals.h"
#include "packet.h"
#include "thermalbalance.h"

namespace kpkt {

void setup_coolinglist(void);
__host__ __device__ void calculate_cooling_rates(int modelgridindex, struct heatingcoolingrates *heatingcoolingrates);
__host__ __device__ double do_kpkt_bb(struct packet *pkt_ptr);
__host__ __device__ double do_kpkt(struct packet *pkt_ptr, double t2, int nts);

__host__ __device__ static inline int get_coolinglistoffset(int element, int ion) {
  return globals::elements[element].ions[ion].coolingoffset;
}

}  // namespace kpkt

#endif  // KPKT_H
