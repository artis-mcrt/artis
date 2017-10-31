#ifndef COMPTON_H
#define COMPTON_H

#include "types.h"

double sig_comp(const PKT *pkt_ptr, double t_current);
void compton_scatter(PKT *pkt_ptr, double t_current);

#endif //COMPTON_H
