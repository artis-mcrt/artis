#ifndef THERMALBALANCE_H
#define THERMALBALANCE_H

#include <vector>

struct HeatingCoolingRates {
  double cooling_collisional{0};
  double cooling_fb{0};
  double cooling_ff{0};
  double cooling_adiabatic{0};
  double heating_collisional{0};
  double heating_bf{0};
  double heating_ff{0};
  double heating_dep{0};
  double dep_frac_heating{0};
  double dep_gamma{0};
  double dep_positron{0};
  double dep_electron{0};
  double dep_alpha{0};
  // analytic rates at the middle of the timestep (t_mid)
  double eps_gamma_ana{0};
  double eps_positron_ana{0};
  double eps_electron_ana{0};
  double eps_alpha_ana{0};
};

void call_T_e_finder(int nonemptymgi, double t_current, double T_min, double T_max,
                     HeatingCoolingRates *heatingcoolingrates, const std::vector<double> &bfheatingcoeffs);
[[nodiscard]] auto get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, double T_R, double W)
    -> double;
void calculate_bfheatingcoeffs(int modelgridindex, std::vector<double> &bfheatingcoeffs);

#endif  // THERMALBALANCE_H
