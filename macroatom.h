#ifndef MACROATOM_H
#define MACROATOM_H

void macroatom_open_file(const int my_rank);
void macroatom_close_file(void);

double do_macroatom(PKT *restrict pkt_ptr, double t1, double t2, int timestep);

double rad_deexcitation_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower, double epsilon_trans, int lineindex, double t_current);
double rad_excitation_ratecoeff(int modelgridindex, int element, int ion, int lower, int upper, double epsilon_trans, int lineindex, double t_current);
double rad_recombination_ratecoeff(float T_e, float nne, int element, int ion, int upper, int lower, int modelgridindex);
double stim_recombination_ratecoeff(float nne, int element, int upperion, int upper, int lower, int modelgridindex);

double col_deexcitation_ratecoeff(float T_e, float nne, double epsilon_trans, int lineindex);
double col_excitation_ratecoeff(float T_e, float nne, int lineindex, double epsilon_trans);
double col_recombination_ratecoeff(int modelgridindex, int element, int upperion, int upper, int lower, double epsilon_trans);
double col_ionization_ratecoeff(float T_e, float nne, int element, int ion, int lower, int phixstargetindex, double epsilon_trans);


#endif //MACROATOM_H
