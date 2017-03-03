#ifndef RATECOEFF_H
#define RATECOEFF_H

void ratecoefficients_init(void);

double alpha_sp_integrand_gsl(double nu, void *restrict paras);
double alpha_sp_E_integrand_gsl(double nu, void *restrict paras);

// double interpolate_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, double T);
// double interpolate_ions_spontrecombcoeff(int element, int ion, double T);

double get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_bfheatingcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_bfcooling(int element, int ion, int level, int phixstargetindex, int modelgridindex);

double calculate_gamma_ion(int modelgridindex, int element, int ion);

#endif //RATECOEFF_H
