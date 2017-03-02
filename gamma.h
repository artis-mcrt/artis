#ifndef GAMMA_H
#define GAMMA_H

#include "types.h"

void pellet_decay(int nts, PKT *pkt_ptr);
double do_gamma(PKT *restrict pkt_ptr, double t1, double t2);
double get_gam_freq(const LIST *restrict line_list, int n);
int get_nul(double freq);

#endif //GAMMA_H
