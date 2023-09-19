#ifndef RPKT_H
#define RPKT_H

void do_rpkt(struct packet *pkt_ptr, double t2);
void emitt_rpkt(struct packet *pkt_ptr);
int closest_transition(double nu_cmf, int next_trans);
double get_rpkt_escape_prob(struct packet *pkt_ptr, double tstart);
double calculate_kappa_bf_gammacontr(int modelgridindex, double nu);
void calculate_kappa_rpkt_cont(const struct packet *pkt_ptr, struct rpkt_cont_opacity *kappa_rpkt_cont_thisthread,
                               const bool usecellhistupdatephixslist);

#endif  // RPKT_H
