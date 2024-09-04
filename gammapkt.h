#ifndef GAMMAPKT_H
#define GAMMAPKT_H

#include "packet.h"

namespace gammapkt {
void init_gamma_data();
void pellet_gamma_decay(Packet &pkt);
void do_gamma(Packet &pkt, int nts, double t2);

}  // namespace gammapkt

#endif  // GAMMAPKT_H
