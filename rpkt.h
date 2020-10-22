#ifndef RPKT_H
#define RPKT_H

#include "types.h"

void do_rpkt(PKT *pkt_ptr, double t2);
void emitt_rpkt(PKT *pkt_ptr);
int closest_transition(double nu_cmf, int next_trans);
double get_rpkt_escape_prob(PKT *pkt_ptr, const double tstart);
void calculate_kappa_bf_gammacontr(const int modelgridindex, const double nu, double *kappa_bf);
void calculate_kappa_rpkt_cont(const PKT *const pkt_ptr);

#ifdef VPKT_ON
int vpkt_call_estimators(PKT *pkt_ptr, double t_current, int realtype);
#endif

#endif //RPKT_H
