#ifndef RADFIELD_H
#define RADFIELD_H

#include "sn3d.h"
#include <gsl/gsl_integration.h>

void radfield_zero_estimators(int modelgridindex);
void radfield_init(int my_rank);
void radfield_write_to_file(int modelgridindex, int timestep);
void radfield_close_file(void);
void radfield_update_estimators(int modelgridindex, double distance_e_cmf, double nu_cmf);
double radfield(double nu, int modelgridindex);
void radfield_fit_parameters(int modelgridindex, int timestep);
void radfield_set_J_normfactor(int modelgridindex, double normfactor);
void radfield_normalise_J(const int modelgridindex, const double estimator_normfactor_over4pi);
void radfield_normalise_nuJ(const int modelgridindex, const double estimator_normfactor_over4pi);
double get_T_R_from_J(const int modelgridindex);
void radfield_titer_J(const int modelgridindex);
void radfield_titer_nuJ(const int modelgridindex);
void radfield_reduce_estimators(void);
void radfield_MPI_Bcast(const int my_rank, const int root, const int root_nstart, const int root_ndo);
void radfield_write_restart_data(FILE *gridsave_file);
void radfield_read_restart_data(FILE *gridsave_file);
int radfield_integrate(
  const gsl_function *f, double nu_a, double nu_b, double epsabs, double epsrel,
  size_t limit, int key, gsl_integration_workspace *workspace, double *result, double *abserr);


inline double radfield_dbb(double nu, float T, float W)
// returns J_nu for a diluted black body
{
  return W * TWOHOVERCLIGHTSQUARED * pow(nu,3) / expm1(HOVERKB * nu / T);
}


#endif //RADFIELD_H
