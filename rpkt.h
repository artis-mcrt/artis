#ifndef RPKT_H
  #define RPKT_H

  double do_rpkt(PKT *pkt_ptr, double t1, double t2);
  double get_event(PKT *pkt_ptr, int *rpkt_eventtype, double t_current, double tau_rnd, double abort_dist);
  int rpkt_event(PKT *pkt_ptr, int rpkt_eventtype, double t_current);
  double closest_transition(PKT *pkt_ptr);
  double min(double a, double b);
  void emitt_rpkt(PKT *pkt_ptr, double t_current);
  double closest_transition_empty(PKT *pkt_ptr);
  int compare_phixslistentry_bynuedge(const void *p1, const void *p2);
  int compare_groundphixslistentry_bynuedge(const void *p1, const void *p2);
  void rpkt_event_thickcell(PKT *pkt_ptr, double t_current);

#endif //RPKT_H
