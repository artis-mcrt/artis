#ifndef GAMMA_H
#define GAMMA_H

#include "types.h"

void init_gamma_linelist(void);
void pellet_decay(int nts, PKT *pkt_ptr);
double do_gamma(PKT *pkt_ptr, const double t1, const double t2, const int timestep);
double get_gam_freq(int n);
int get_nul(double freq);

#endif //GAMMA_H
