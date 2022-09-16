#ifndef RADFIELD_H
#define RADFIELD_H

#include <gsl/gsl_integration.h>

#include <cstdio>

#include "sn3d.h"

namespace radfield {
void zero_estimators(int modelgridindex);
void init(int my_rank, int ndo, int ndo_nonempty);
void initialise_prev_titer_photoionestimators(void);
void write_to_file(int modelgridindex, int timestep);
void close_file(void);
__host__ __device__ void update_estimators(int modelgridindex, double distance_e_cmf, double nu_cmf,
                                           const struct packet *pkt_ptr, double t_current);
void increment_lineestimator(int modelgridindex, int lineindex, double increment);
__host__ __device__ double radfield(double nu, int modelgridindex);
__host__ __device__ double dbb_mgi(double nu, int modelgridindex);
void fit_parameters(int modelgridindex, int timestep);
__host__ __device__ void set_J_normfactor(int modelgridindex, double normfactor);
__host__ __device__ void normalise_J(int modelgridindex, double estimator_normfactor_over4pi);
__host__ __device__ void normalise_nuJ(int modelgridindex, double estimator_normfactor_over4pi);
__host__ __device__ double get_T_J_from_J(int modelgridindex);
__host__ __device__ int get_Jblueindex(int lineindex);
__host__ __device__ double get_Jb_lu(int modelgridindex, int jblueindex);
__host__ __device__ int get_Jb_lu_contribcount(int modelgridindex, int jblueindex);
void titer_J(int modelgridindex);
void titer_nuJ(int modelgridindex);
void reduce_estimators(void);
void do_MPI_Bcast(int modelgridindex, int root, int root_node_id);
void write_restart_data(FILE *gridsave_file);
void read_restart_data(FILE *gridsave_file);
__host__ __device__ void normalise_bf_estimators(int modelgridindex, double estimator_normfactor_over_H);
double get_bfrate_estimator(int element, int lowerion, int lower, int phixstargetindex, int modelgridindex);
void print_bfrate_contributions(int element, int lowerion, int lower, int phixstargetindex, int modelgridindex,
                                double nnlowerlevel, double nnlowerion);
void reset_bfrate_contributions(const int modelgridindex);
int integrate(const gsl_function *f, double nu_a, double nu_b, double epsabs, double epsrel, size_t limit, int key,
              gsl_integration_workspace *workspace, double *result, double *abserr);

template <typename numtype>
__host__ __device__ static inline double dbb(double nu, numtype T, numtype W)
// returns J_nu for a dilute black body [ergs/s/sr/cm2/Hz]
{
  return W * TWOHOVERCLIGHTSQUARED * pow(nu, 3) / expm1(HOVERKB * nu / T);
}
}  // namespace radfield

#endif  // RADFIELD_H
