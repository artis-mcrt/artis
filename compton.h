#ifndef COMPTON_H
#define COMPTON_H

double sig_comp(PKT *pkt_ptr, double t_current);
int com_sca(PKT *pkt_ptr, double t_current);
double sigma_compton_partial(double x, double f);
double choose_f(double xx, double zrand);
double thomson_angle();

#endif //COMPTON_H
