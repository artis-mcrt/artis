#ifndef GAMMAPKT_H
#define GAMMAPKT_H

#include "constants.h"
#include "packet.h"
#include "random.h"

namespace gammapkt {
void init_gamma_data();
DEVICE_FUNC void pellet_gamma_decay(Packet& pkt);
DEVICE_FUNC void do_gamma(Packet& pkt, int nts, double t2);
auto choose_gamma_ray(int nucindex, rngstate_type& rngstate) -> double;

}  // namespace gammapkt

#endif  // GAMMAPKT_H
