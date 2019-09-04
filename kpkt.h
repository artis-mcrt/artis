#ifndef KPKT_H
#define KPKT_H

void calculate_cooling_rates(int modelgridindex, heatingcoolingrates_t *heatingcoolingrates);
double do_kpkt_bb(PKT *restrict pkt_ptr, const double t1, const double t2, const int nts);
double do_kpkt(PKT *restrict pkt_ptr, const double t1, const double t2, const int nts);

inline int get_coolinglistoffset(int element, int ion)
{
  return elements[element].ions[ion].coolingoffset;
}

#endif //KPKT_H
