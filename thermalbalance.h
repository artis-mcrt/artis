// The HeatingCoolingRates struct and declarations for the electron thermal balance solver
// (thermalbalance.cc).

#ifndef THERMALBALANCE_H
#define THERMALBALANCE_H

#include <span>

// Per-cell energy budget of the free electrons. All members are erg/s/cm^3 except dep_frac_heating.
// call_T_e_finder() solves heating_ff + heating_bf + heating_collisional + heating_dep
//                        = cooling_ff + cooling_fb + cooling_collisional + cooling_adiabatic + cooling_heatcapacity
// for T_e. Only the heating_* and cooling_* members are terms in that balance. The dep_* members feed it (via
// heating_dep), dep_frac_heating is computed from them, and the eps_*_ana members are reported to
// estimators.out and used to set some dep_* values, but never enter the balance directly.
struct HeatingCoolingRates {
  double cooling_collisional{0};
  double cooling_fb{0};
  double cooling_ff{0};
  double cooling_adiabatic{0};  // work done by the expanding ejecta, p * (dV/dt) / V
  // the change of the thermal energy density 3/2 k_B n_tot T_e from the previous grid update over the time
  // interval (backward Euler). Zero without NLTE_TIME_DEPENDENT_FIRST_TIMESTEP and in a steady-state timestep.
  double cooling_heatcapacity{0};

  double heating_collisional{0};
  double heating_bf{0};
  double heating_ff{0};
  double heating_dep{0};  // the share of the deposition below that goes to heating

  // heating_dep over the total deposition (leptons + alpha + spontaneous fission). The lepton share comes
  // from the Spencer-Fano solution; alpha and spontaneous fission are assumed to thermalise completely.
  double dep_frac_heating{0};

  // Deposition actually seen in this cell. dep_gamma is always from the Monte Carlo estimators and so noisy,
  // as are the positron/electron/alpha rates unless PARTICLE_THERMALISATION_SCHEME is INSTANTFULLDEPOSITION,
  // which copies them from the analytic rates below. dep_spfission is always eps_spfission_ana.
  double dep_gamma{0};
  double dep_positron{0};
  double dep_electron{0};
  double dep_alpha{0};
  double dep_spfission{0};

  // Analytic rates at t_mid (density times decay power per mass): what the cell would absorb if everything
  // thermalised locally and instantly, so dep_* against eps_*_ana shows how much escapes or moves elsewhere.
  double eps_gamma_ana{0};
  double eps_positron_ana{0};
  double eps_electron_ana{0};
  double eps_alpha_ana{0};
  double eps_spfission_ana{0};

  // the two sides of the thermal balance. The order of the terms is fixed, because the results must not change.
  [[nodiscard]] auto get_total_heating() const -> double {
    return heating_ff + heating_bf + heating_collisional + heating_dep;
  }
  [[nodiscard]] auto get_total_cooling() const -> double {
    return cooling_ff + cooling_fb + cooling_collisional + cooling_adiabatic + cooling_heatcapacity;
  }
};

void call_T_e_finder(int nonemptymgi, double t_current, HeatingCoolingRates& heatingcoolingrates,
                     std::span<const double> bfheatingcoeffs);
auto calculate_bfheatingcoeff(int element, int ion, int level, int phixstargetindex, int nonemptymgi) -> double;
void calculate_bfheatingcoeffs(int nonemptymgi, std::span<double> bfheatingcoeffs);

#endif  // THERMALBALANCE_H
