#ifndef GREY_EMISSIVITIES_H
#define GREY_EMISSIVITIES_H

#include "types.h"

void rlc_emiss_gamma(const PKT *pkt_ptr, double dist);
void rlc_emiss_rpkt(const PKT *pkt_ptr, double dist);
void normalise_grey(int nts);
void write_grey(int nts);

#endif //GREY_EMISSIVITIES_H
