#ifndef RATECOEFF_H
#define RATECOEFF_H

void ratecoefficients_init(void);

double alpha_sp_integrand_gsl(double nu, void *restrict paras);
double alpha_sp_E_integrand_gsl(double nu, void *restrict paras);

double interpolate_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, double T);

double get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, float T_e);
double get_bfheatingcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_bfcooling(int element, int ion, int level, int phixstargetindex, int modelgridindex);

double calculate_gamma_ion_per_gspop(int modelgridindex, int element, int ion);
double calculate_gamma_ion_per_ionpop(
  const int modelgridindex, const float T_e, const int element, const int lowerion,
  const bool assume_lte, const bool collisional_not_radiative, const bool printdebug);

double calculate_recombcoeff_ion_per_gmpop(
  const int modelgridindex, const float T_e,
  const int element, const int upperion,
  const bool assume_lte, const bool collisional_not_radiative, const bool printdebug,
  const bool lower_superlevel_only);

double calculate_recombcoeff_ion_per_ionpop(
  const int modelgridindex, const float T_e,
  const int element, const int upperion,
  const bool assume_lte, const bool collisional_not_radiative, const bool printdebug,
  const bool lower_superlevel_only);

#endif //RATECOEFF_H