#ifndef KPKT_H
  #define KPKT_H

  void calculate_kpkt_rates(int modelgridindex);

  void calculate_kpkt_rates_ion(int modelgridindex, int element, int ion, int low, double oldcoolingsum, int high);

  double do_kpkt_bb(PKT *pkt_ptr, double t1, double t2);
  double sample_planck(double T);
  double planck(double nu, double T);
  double do_kpkt(PKT *pkt_ptr, double t1, double t2, int nts);
  int get_coolinglistoffset(int element, int ion);
  int get_ncoolingterms(int element, int ion);

#endif //KPKT_H
