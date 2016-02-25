#ifndef GREY_EMISSIVITIES_H
  #define GREY_EMISSIVITIES_H

  int rlc_emiss_gamma(PKT *pkt_ptr, double dist, double t_current);
  int rlc_emiss_rpkt(PKT *pkt_ptr, double dist, double t_current);
  int normalise_grey(int nts);
  int write_grey(int nts);
  double meanf_sigma(double x);
  int emiss_rlc_load(int nts);
  int grey_rt(RAY *ray_ptr, int nray, double ldist, double *single_pos, double single_t, int lindex);

#endif //GREY_EMISSIVITIES_H
