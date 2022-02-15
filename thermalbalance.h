#ifndef THERMALBALANCE_H
#define THERMALBALANCE_H

void call_T_e_finder(int modelgridindex, int timestep, double t_current, double T_min, double T_max, heatingcoolingrates_t *heatingcoolingrates);
double get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, double T, double W);
void calculate_bfheatingcoeffs(int modelgridindex);

#endif //THERMALBALANCE_H
