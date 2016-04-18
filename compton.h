#ifndef COMPTON_H
#define COMPTON_H

#include "types.h"

double sig_comp(const PKT *pkt_ptr, double t_current);
int com_sca(PKT *pkt_ptr, double t_current);
double sigma_compton_partial(double x, double f);
double choose_f(double xx, double zrand);
double thomson_angle(void);

#endif //COMPTON_H
