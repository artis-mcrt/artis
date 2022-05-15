#ifndef RPKT_H
#define RPKT_H

#include "types.h"

__host__ __device__ void do_rpkt(PKT *pkt_ptr, double t2);
__host__ __device__ void emitt_rpkt(PKT *pkt_ptr);
__host__ __device__ int closest_transition(double nu_cmf, int next_trans);
__host__ __device__ double get_rpkt_escape_prob(PKT *pkt_ptr, const double tstart);
__host__ __device__ double calculate_kappa_bf_gammacontr(const int modelgridindex, const double nu);
__host__ __device__ void calculate_kappa_rpkt_cont(const PKT *const pkt_ptr, rpkt_cont_opacity_struct *kappa_rpkt_cont_thisthread);

#endif //RPKT_H
