#ifndef RATECOEFF_H
#define RATECOEFF_H

void tabulate_ratecoefficients_gsl();
void calculate_ion_alpha_sp();
double alpha_sp_integrand_gsl(double nu, void *paras);
double alpha_sp_E_integrand_gsl(double nu, void *paras);
double gamma_integrand_gsl(double nu, void *paras);
double gammacorr_integrand_gsl(double nu, void *paras);
double approx_bfheating_integrand_gsl(double nu, void *paras);
double bfcooling_integrand_gsl(double nu, void *paras);
double interpolate_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, double T);

double interpolate_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, double T);
double interpolate_ions_spontrecombcoeff(int element, int ion, double T);

double get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_bfheatingcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_bfheatingcoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_bfcooling(int element, int ion, int level, int phixstargetindex, int modelgridindex);

#endif //RATECOEFF_H
