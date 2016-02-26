#ifndef THERMALBALANCE_H
#define THERMALBALANCE_H

double call_T_e_finder(int modelgridindex, double t_current, int tb_info, double T_min, double T_max);
double find_T_e(double T_e, void *paras);
void calculate_heating_rates(int modelgridindex);
void calculate_cooling_rates(int modelgridindex);

#endif //THERMALBALANCE_H
