#ifndef THERMALBALANCE_H
#define THERMALBALANCE_H

struct heatingcoolingrates {
  double cooling_collisional;
  double cooling_fb;
  double cooling_ff;
  double cooling_adiabatic;
  double heating_collisional;
  double heating_bf;
  double heating_ff;
  double heating_dep;
  double nt_frac_heating;
};

void call_T_e_finder(int modelgridindex, int timestep, double t_current, double T_min, double T_max,
                     struct heatingcoolingrates *heatingcoolingrates);
double get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, double T, double W);
void calculate_bfheatingcoeffs(int modelgridindex);

#endif  // THERMALBALANCE_H
