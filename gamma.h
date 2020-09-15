#ifndef GAMMA_H
#define GAMMA_H

#include "types.h"

void init_gamma_linelist(void);
void pellet_decay(int nts, PKT *pkt_ptr);
double do_gamma(PKT *pkt_ptr, double t2);
double get_gam_freq(int n);
int get_nul(double freq);

#endif //GAMMA_H
