#ifndef RPKT_H
#define RPKT_H

#include "types.h"

double do_rpkt(PKT *restrict pkt_ptr, double t1, double t2);
void emitt_rpkt(PKT *restrict pkt_ptr, double t_current);
void calculate_kappa_rpkt_cont(const PKT *restrict const pkt_ptr, const double t_current);
void calculate_kappa_vpkt_cont(const PKT *pkt_ptr, double t_current);
int compare_phixslistentry_bynuedge(const void *restrict p1, const void *restrict p2);
int compare_groundphixslistentry_bynuedge(const void *restrict p1, const void *restrict p2);

#endif //RPKT_H
