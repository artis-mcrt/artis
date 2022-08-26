#ifndef GAMMA_H
#define GAMMA_H

#include "packet.h"

namespace gamma {
void init_gamma_linelist(void);
void pellet_gamma_decay(int nts, struct packet *pkt_ptr);
void do_gamma(struct packet *pkt_ptr, double t2);
double get_gam_freq(int n);
int get_nul(double freq);
}  // namespace gamma

#endif  // GAMMA_H
