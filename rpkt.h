#ifndef RPKT_H
#define RPKT_H

#include "types.h"

double do_rpkt(PKT *pkt_ptr, double t1, double t2);
double closest_transition(PKT *pkt_ptr);

static inline
double min(double a, double b)
// returns minimum of a and b
{
  if (a >= b)
  {
    return(b);
  }
  else
  {
    return(a);
  }
}

double closest_transition_empty(PKT *pkt_ptr);
void emitt_rpkt(PKT *pkt_ptr, double t_current);
void calculate_kappa_rpkt_cont(const PKT *pkt_ptr, double t_current);
void calculate_kappa_vpkt_cont(const PKT *pkt_ptr, double t_current);
int compare_phixslistentry_bynuedge(const void *p1, const void *p2);
int compare_groundphixslistentry_bynuedge(const void *p1, const void *p2);
double do_rpkt_thickcell(PKT *pkt_ptr, double t1, double t2);

#endif //RPKT_H
