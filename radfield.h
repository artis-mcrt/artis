#ifndef RADFIELD_H
#define RADFIELD_H

#include <stdio.h>
#include <gsl/gsl_integration.h>
#include "types.h"
#include "sn3d.h"

void radfield_zero_estimators(int modelgridindex);
void radfield_jblue_init(void);
void radfield_init(int my_rank);
void initialise_prev_titer_photoionestimators(void);
void radfield_write_to_file(int modelgridindex, int timestep);
void radfield_close_file(void);
__host__ __device__ void radfield_update_estimators(int modelgridindex, double distance_e_cmf, double nu_cmf, const PKT *pkt_ptr, double t_current);
__host__ __device__ void radfield_increment_lineestimator(int modelgridindex, int lineindex, double increment);
__host__ __device__ double radfield(double nu, int modelgridindex);
__host__ __device__ double radfield_dbb_mgi(double nu, int modelgridindex);
void radfield_fit_parameters(int modelgridindex, int timestep);
void radfield_set_J_normfactor(int modelgridindex, double normfactor);
void radfield_normalise_J(int modelgridindex, double estimator_normfactor_over4pi);
void radfield_normalise_nuJ(int modelgridindex, double estimator_normfactor_over4pi);
double get_T_R_from_J(int modelgridindex);
__host__ __device__ int radfield_get_Jblueindex(int lineindex);
__host__ __device__ double radfield_get_Jb_lu(int modelgridindex, int jblueindex);
__host__ __device__ int radfield_get_Jb_lu_contribcount(int modelgridindex, int jblueindex);
void radfield_titer_J(int modelgridindex);
void radfield_titer_nuJ(int modelgridindex);
void radfield_reduce_estimators(void);
void radfield_MPI_Bcast(int modelgridindex, int root);
void radfield_write_restart_data(FILE *gridsave_file);
void radfield_read_restart_data(FILE *gridsave_file);
void radfield_normalise_bf_estimators(int modelgridindex, double estimator_normfactor_over_H);
__host__ __device__ double get_bfrate_estimator(int element, int lowerion, int lower, int phixstargetindex, int modelgridindex);
void print_bfrate_contributions(int element, int lowerion, int lower, int phixstargetindex, int modelgridindex, double nnlowerlevel, double nnlowerion);
void reset_bfrate_contributions(const int modelgridindex);
// int radfield_integrate(
//   const gsl_function *f, double nu_a, double nu_b, double epsabs, double epsrel,
//   size_t limit, int key, gsl_integration_workspace *workspace, double *result, double *abserr);

__host__ __device__
inline double radfield_dbb(double nu, float T, float W)
// returns J_nu for a diluted black body
{
  return W * TWOHOVERCLIGHTSQUARED * pow(nu, 3) / expm1(HOVERKB * nu / T);
}


#endif //RADFIELD_H
