#ifndef GREY_EMISSIVITIES_H
#define GREY_EMISSIVITIES_H

#include "types.h"

int rlc_emiss_gamma(const PKT *pkt_ptr, double dist, double t_current);
int rlc_emiss_rpkt(const PKT *pkt_ptr, double dist, double t_current);
int normalise_grey(int nts);
int write_grey(int nts);
int grey_rt(RAY *ray_ptr, int nray, double ldist, double *single_pos, double single_t, int lindex);

#endif //GREY_EMISSIVITIES_H
