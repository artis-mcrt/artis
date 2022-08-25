#ifndef EMISSIVITIES_H
#define EMISSIVITIES_H

#include "types.h"

void compton_emiss_cont(const struct packet *pkt_ptr, double dist);
void pp_emiss_cont(const struct packet *pkt_ptr, double dist);
void zero_estimators(void);
void normalise_compton_estimators(const int nts);
void write_compton_estimators(int nts);
bool estim_switch(int nts);

#endif //EMISSIVITIES_H
