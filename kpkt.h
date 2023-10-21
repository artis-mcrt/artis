#ifndef KPKT_H
#define KPKT_H

#include "globals.h"
#include "packet.h"
#include "thermalbalance.h"

constexpr double COOLING_UNDEFINED = -99;

namespace kpkt {

void setup_coolinglist();
void calculate_cooling_rates(int modelgridindex, struct heatingcoolingrates *heatingcoolingrates);
void do_kpkt_blackbody(struct packet *pkt_ptr);
void do_kpkt(struct packet *pkt_ptr, double t2, int nts);

static inline int get_coolinglistoffset(int element, int ion) {
  return globals::elements[element].ions[ion].coolingoffset;
}

}  // namespace kpkt

#endif  // KPKT_H
