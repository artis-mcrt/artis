#pragma once
#ifndef THERMALBALANCE_H
#define THERMALBALANCE_H

#include <vector>
struct HeatingCoolingRates {
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
                     HeatingCoolingRates *heatingcoolingrates, const std::vector<double> &bfheatingcoeffs);
[[nodiscard]] auto get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, double T_R, double W)
    -> double;
void calculate_bfheatingcoeffs(int modelgridindex, std::vector<double> &bfheatingcoeffs);

#endif  // THERMALBALANCE_H
