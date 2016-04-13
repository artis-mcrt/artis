#ifndef EMISSIVITIES_H
#define EMISSIVITIES_H

int add_gam_line_emissivity(RAY *ray_ptr, int nray, double *single_pos, double single_t, int lindex, double dnuds);
int continuum_rt(RAY *ray_ptr, int nray, double ldist, double *single_pos, double single_t, int lindex);
int compton_emiss_cont(const PKT *pkt_ptr, double dist, double t_current);
int pp_emiss_cont(const PKT *pkt_ptr, double dist, double t_current);
int zero_estimators(void);
int normalise_estimators(int nts);
int write_estimators(int nts);
int estim_switch(int nts);
int emiss_load(int nts);

#endif //EMISSIVITIES_H
