#ifndef MACROATOM_H
#define MACROATOM_H

#include "grid_init.h"
#include "sn3d.h"
#include "types.h"

double do_ma(PKT *restrict pkt_ptr, double t1, double t2, int timestep);

double rad_deexcitation_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower, double epsilon_trans, int lineindex, double t_current);
double rad_excitation_ratecoeff(int modelgridindex, int element, int ion, int lower, int upper, double epsilon_trans, int lineindex, double t_current);
double rad_recombination_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower);

double col_deexcitation_ratecoeff(float T_e, float nne, double epsilon_trans, int lineindex);
double col_excitation_ratecoeff(float T_e, float nne, int lineindex, double epsilon_trans);
double col_recombination_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower, double epsilon_trans);
double col_ionization_ratecoeff(float T_e, float nne, int element, int ion, int lower, int phixstargetindex, double epsilon_trans);


#endif //MACROATOM_H
