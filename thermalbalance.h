#ifndef THERMALBALANCE_H
#define THERMALBALANCE_H

double call_T_e_finder(int modelgridindex, double t_current, double T_min, double T_max);
void calculate_heating_rates(int modelgridindex);
void calculate_cooling_rates(int modelgridindex);

#endif //THERMALBALANCE_H
