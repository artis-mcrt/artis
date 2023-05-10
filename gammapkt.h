#ifndef GAMMAPKT_H
#define GAMMAPKT_H

#include "packet.h"

namespace gammapkt {
void init_gamma_linelist();
void pellet_gamma_decay(struct packet *pkt_ptr);
void do_gamma(struct packet *pkt_ptr, double t2);
void normalise_grey(int nts);

}  // namespace gammapkt

#endif  // GAMMAPKT_H
