#pragma once
#ifndef GAMMAPKT_H
#define GAMMAPKT_H

#include "packet.h"

namespace gammapkt {
void init_gamma_linelist();
void pellet_gamma_decay(Packet &pkt);
void do_gamma(Packet &pkt, double t2);
void normalise(int nts);

}  // namespace gammapkt

#endif  // GAMMAPKT_H
