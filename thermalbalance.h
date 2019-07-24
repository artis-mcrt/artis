#ifndef THERMALBALANCE_H
#define THERMALBALANCE_H

void call_T_e_finder(int modelgridindex, int timestep, double t_current, double T_min, double T_max, heatingcoolingrates_t *heatingcoolingrates);

#endif //THERMALBALANCE_H
