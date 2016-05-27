#ifndef EMISSIVITIES_H
#define EMISSIVITIES_H

#include "types.h"

void add_gam_line_emissivity(RAY *ray_ptr, int nray, double *single_pos, double single_t, int lindex, double dnuds);
void continuum_rt(RAY *ray_ptr, int nray, double ldist, double *single_pos, double single_t, int lindex);
void compton_emiss_cont(const PKT *pkt_ptr, double dist, double t_current);
void pp_emiss_cont(const PKT *pkt_ptr, double dist, double t_current);
void zero_estimators(void);
void normalise_estimators(int nts);
void write_estimators(int nts);
bool estim_switch(int nts);
void emiss_load(int nts);

#endif //EMISSIVITIES_H
