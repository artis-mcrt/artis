#ifndef RPKT_H
#define RPKT_H

#include "types.h"

struct rpkt_cont_opacity_struct
{
  double nu; // frequency at which opacity was calculated
  double total;
  double es;
  double ff;
  double bf;
  double fb;
  double bf_inrest;
  double fb_inrest;
  double ffheating;
  //double bfheating;
  int modelgridindex;
  int timestep;
};

double do_rpkt(PKT *restrict pkt_ptr, double t1, double t2);
void emitt_rpkt(PKT *restrict pkt_ptr, double t_current);
int closest_transition(double nu_cmf, int next_trans);
double get_rpkt_escape_prob(PKT *restrict pkt_ptr, const double tstart);
void calculate_kappa_bf_fb_gammacontr(const int modelgridindex, const double nu, double *kappa_bf, double *kappa_fb);
void calculate_kappa_vpkt_cont(const PKT *pkt_ptr, double t_current, struct rpkt_cont_opacity_struct *kappa_rpkt_cont_thisthread);


#endif //RPKT_H
