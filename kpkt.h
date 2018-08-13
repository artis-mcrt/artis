#ifndef KPKT_H
#define KPKT_H

void calculate_kpkt_rates(int modelgridindex);
double do_kpkt_bb(PKT *restrict pkt_ptr, double t1);
double do_kpkt(PKT *restrict pkt_ptr, double t1, double t2, int nts);

inline int get_coolinglistoffset(int element, int ion)
{
  return elements[element].ions[ion].coolingoffset;
}

#endif //KPKT_H
