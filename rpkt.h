#ifndef RPKT_H
#define RPKT_H

#include "types.h"

double do_rpkt(PKT *pkt_ptr, double t1, double t2);
void emitt_rpkt(PKT *pkt_ptr);
int closest_transition(double nu_cmf, int next_trans);
double get_rpkt_escape_prob(PKT *pkt_ptr, const double tstart);
void calculate_kappa_bf_fb_gammacontr(const int modelgridindex, const double nu, double *kappa_bf, double *kappa_fb);
void calculate_kappa_rpkt_cont(const PKT *const pkt_ptr, const int modelgridindex);
void calculate_kappa_vpkt_cont(const PKT *pkt_ptr, double t_current);

#endif //RPKT_H
