#ifndef GAMMA_H
  #define GAMMA_H

  int pellet_decay(int nts, PKT *pkt_ptr);
  int choose_gamma_ray(PKT *pkt_ptr);
  double do_gamma (PKT *pkt_ptr, double t1, double t2);

#endif //GAMMA_H
