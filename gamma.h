#ifndef GAMMA_H
#define GAMMA_H

#include "types.h"

void init_gamma_linelist(void);
double pellet_decay_gamma(PKT *pkt_ptr, const double tdecay, const double t2, const int nts);
double do_gamma(PKT *pkt_ptr, const double t1, const double t2, const int timestep);
double get_gam_freq(int n);
int get_nul(double freq);

#endif //GAMMA_H
