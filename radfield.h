#ifndef RADFIELD_H
#define RADFIELD_H

#include <gsl/gsl_integration.h>

#include <cstdio>

#include "sn3d.h"

namespace radfield {

void zero_estimators(int modelgridindex);
void init(int my_rank, int ndo_nonempty);
void initialise_prev_titer_photoionestimators();
void write_to_file(int modelgridindex, int timestep);
void close_file();
void update_estimators(int modelgridindex, double distance_e_cmf, double nu_cmf, const struct packet *pkt_ptr);
void update_lineestimator(int modelgridindex, int lineindex, double increment);
[[nodiscard]] auto radfield(double nu, int modelgridindex) -> double;
void fit_parameters(int modelgridindex, int timestep);
void set_J_normfactor(int modelgridindex, double normfactor);
void normalise_J(int modelgridindex, double estimator_normfactor_over4pi);
void normalise_nuJ(int modelgridindex, double estimator_normfactor_over4pi);
[[nodiscard]] auto get_T_J_from_J(int modelgridindex) -> double;
[[nodiscard]] auto get_Jblueindex(int lineindex) -> int;
[[nodiscard]] auto get_Jb_lu(int modelgridindex, int jblueindex) -> double;
[[nodiscard]] auto get_Jb_lu_contribcount(int modelgridindex, int jblueindex) -> int;
void titer_J(int modelgridindex);
void titer_nuJ(int modelgridindex);
void reduce_estimators();
void do_MPI_Bcast(int modelgridindex, int root, int root_node_id);
void write_restart_data(FILE *gridsave_file);
void read_restart_data(FILE *gridsave_file);
void normalise_bf_estimators(int modelgridindex, double estimator_normfactor_over_H);
auto get_bfrate_estimator(int element, int lowerion, int lower, int phixstargetindex, int modelgridindex) -> double;
void print_bfrate_contributions(int element, int lowerion, int lower, int phixstargetindex, int modelgridindex,
                                double nnlowerlevel, double nnlowerion);
void reset_bfrate_contributions(int modelgridindex);
auto integrate(const gsl_function *f, double nu_a, double nu_b, double epsabs, double epsrel, size_t limit, int key,
               gsl_integration_workspace *workspace, double *result, double *abserr) -> int;

[[nodiscard]] constexpr auto dbb(double nu, auto T, auto W) -> double
// returns J_nu [ergs/s/sr/cm2/Hz] for a dilute black body with temperature T and dilution factor W
{
  return W * TWOHOVERCLIGHTSQUARED * std::pow(nu, 3) / std::expm1(HOVERKB * nu / T);
}

}  // namespace radfield

#endif  // RADFIELD_H
