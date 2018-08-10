#ifndef GAMMA_H
#define GAMMA_H

#include "types.h"

void get_gam_ll(void);
void pellet_decay(int nts, PKT *pkt_ptr);
double do_gamma(PKT *pkt_ptr, double t1, double t2);
double get_gam_freq(int n);
int get_nul(double freq);
void read_decaydata(void);

#endif //GAMMA_H
