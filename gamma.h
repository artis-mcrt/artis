#ifndef GAMMA_H
#define GAMMA_H

#include "types.h"

void init_gamma_linelist(void);
__host__ __device__ void pellet_decay(int nts, PKT *pkt_ptr);
__host__ __device__ double do_gamma(PKT *pkt_ptr, double t1, double t2, int tid);
__host__ __device__ double get_gam_freq(int n);
__host__ __device__ int get_nul(double freq);

#endif //GAMMA_H
