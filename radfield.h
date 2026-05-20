#ifndef RADFIELD_H
#define RADFIELD_H

#include <cmath>
#include <cstddef>
#include <cstdio>

#include "constants.h"
#include "rpkt.h"

namespace radfield {

void zero_estimators();
void init();
void initialise_prev_titer_photoionestimators();
DEVICE_FUNC void update_estimators(ptrdiff_t nonemptymgi, double distance_e_cmf, double nu_cmf,
                                   double doppler_nucmf_on_nurf, const Phixslist& phixslist, bool thickcell);
DEVICE_FUNC void update_lineestimator(int nonemptymgi, int lineindex, double increment);
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto radfield(double nu, int nonemptymgi) -> double;
void fit_parameters(int nonemptymgi, int timestep);
void set_J_normfactor(int nonemptymgi, double normfactor);
void normalise_J(int nonemptymgi, double estimator_normfactor_over4pi);
void normalise_nuJ(int nonemptymgi, double estimator_normfactor_over4pi);
[[nodiscard]] auto get_T_J_from_J(int nonemptymgi) -> float;
[[nodiscard]] auto get_Jblueindex(int lineindex) -> int;
[[nodiscard]] auto get_Jb_lu(int nonemptymgi, int jblueindex) -> double;
void titer_J(int nonemptymgi);
void titer_nuJ(int nonemptymgi);
void reduce_estimators();
void do_MPI_Bcast(ptrdiff_t nonemptymgi, int root, int root_node_id);
void write_restart_data(FILE* gridsave_file);
void read_restart_data(FILE* gridsave_file);
void normalise_bf_estimators(int nts, int nts_prev, int titer, double deltat);
[[nodiscard]] DEVICE_FUNC auto get_bfrate_estimator(int element, int lowerion, int lower, int phixstargetindex,
                                                    int nonemptymgi) -> double;

// get J_nu [ergs/s/sr/cm2/Hz] for a dilute black body with temperature T
[[gnu::const]] [[nodiscard]] constexpr auto planck(const double nu, const double T) -> double {
  return (2 * H / CLIGHTSQUARED) * pow3(nu) / std::expm1(HOVERKB * nu / T);
}

}  // namespace radfield

#endif  // RADFIELD_H
