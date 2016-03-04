#ifndef KPKT_H
#define KPKT_H

void calculate_kpkt_rates(int modelgridindex);
void calculate_kpkt_rates_ion(int modelgridindex, int element, int ion, int low, double oldcoolingsum, int high);
double do_kpkt_bb(PKT *pkt_ptr, double t1, double t2);
double sample_planck(double T);

inline
double planck(double nu, double T)
/// returns intensity for frequency nu and temperature T according
/// to the Planck distribution
{
  return TWOHOVERCLIGHTSQUARED * pow(nu,3) / (exp(HOVERKB*nu/T) - 1);
}

double do_kpkt(PKT *pkt_ptr, double t1, double t2, int nts);

inline
int get_coolinglistoffset(int element, int ion)
{
  return elements[element].ions[ion].coolingoffset;
}

inline
int get_ncoolingterms(int element, int ion)
{
  return elements[element].ions[ion].ncoolingterms;
}

#endif //KPKT_H
