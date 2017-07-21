#ifndef RATECOEFF_H
#define RATECOEFF_H

void ratecoefficients_init(void);

double alpha_sp_integrand_gsl(double nu, void *restrict paras);
double alpha_sp_E_integrand_gsl(double nu, void *restrict paras);

double interpolate_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, double T);
double interpolate_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, double T);

double get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, float T_e);
double get_bfheatingcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_bfcooling(int element, int ion, int level, int phixstargetindex, int modelgridindex);

double calculate_iongamma_per_gspop(int modelgridindex, int element, int ion);
double calculate_iongamma_per_ionpop(
  int modelgridindex, float T_e, int element, int lowerion,
  bool assume_lte, bool collisional_not_radiative, bool printdebug);

double calculate_ionrecombcoeff(
  int modelgridindex, float T_e,
  int element, int upperion,
  bool assume_lte, bool collisional_not_radiative, bool printdebug,
  bool lower_superlevel_only, bool per_groundmultipletpop);

#endif //RATECOEFF_H