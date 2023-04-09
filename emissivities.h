#ifndef EMISSIVITIES_H
#define EMISSIVITIES_H

void compton_emiss_cont(const struct packet *pkt_ptr, double dist);
void pp_emiss_cont(const struct packet *pkt_ptr, double dist);
void zero_estimators();

#endif  // EMISSIVITIES_H
