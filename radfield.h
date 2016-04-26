#ifndef RADFIELD_H
#define RADFIELD_H

#include "types.h"
#include <stdbool.h>

typedef enum
{
  ONE = 0,
  TIMES_NU = 1,
  TIMES_E = 2
} enum_prefactor;


void zero_radfield_estimators(int modelgridindex);
void update_radfield(int modelgridindex, double distance, double e_cmf,
                     double nu_cmf);
double radfield(double nu, int modelgridindex);
void fit_radfield_parameters(int modelgridindex);
double call_T_R_finder(int modelgridindex, int binindex);
double delta_nu_bar(double T_R, void *paras);
double calculate_planck_integral(double T_R, double nu_lower, double nu_upper,
                                 enum_prefactor prefactor);
double gsl_integrand_planck(double nu, void *paras);
void radfield_set_J_normfactor(int modelgridindex, double normfactor);

#endif //RADFIELD_H
