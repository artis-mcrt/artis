#ifndef THERMALBALANCE_H
#define THERMALBALANCE_H

void call_T_e_finder(int modelgridindex, int timestep, double t_current, double T_min, double T_max, heatingcoolingrates_t *heatingcoolingrates, int tid);
double get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, float T, double W);
void calculate_bfheatingcoeffs(int modelgridindex, int tid);


#endif //THERMALBALANCE_H
