#ifndef KPKT_H
#define KPKT_H

#include "types.h"

void calculate_kpkt_rates(int modelgridindex);
double do_kpkt_bb(PKT *pkt_ptr, double t1, double t2);
double planck(double nu, double T);
double do_kpkt(PKT *pkt_ptr, double t1, double t2, int nts);
int get_coolinglistoffset(int element, int ion);
int get_ncoolingterms(int element, int ion);

#endif //KPKT_H
