#ifndef MOVE_H
  #define MOVE_H

  void update_estimators(PKT *pkt_ptr, double distance);
  int move_pkt(PKT *pkt_ptr, double distance, double time);

#endif //MOVE_H
