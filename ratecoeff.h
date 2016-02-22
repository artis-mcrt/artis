#ifndef RATECOEFF_H
  #define RATECOEFF_H

  double alpha_sp_integrand_gsl(double nu, void *paras);
  double alpha_sp_E_integrand_gsl(double nu, void *paras);
  double gamma_integrand_gsl(double nu, void *paras);
  double gammacorr_integrand_gsl(double nu, void *paras);
  double approx_bfheating_integrand_gsl(double nu, void *paras);
  double bfcooling_integrand_gsl(double nu, void *paras);
  double interpolate_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, double T);

#endif //RATECOEFF_H
