#ifndef NONTHERMAL_H
#define NONTHERMAL_H

#include <cstddef>
#include <cstdio>

#include "constants.h"
#include "packet.h"
#include "thermalbalance.h"

namespace nonthermal {
void init();
void close_file();
void solve_spencerfano(int nonemptymgi, int timestep, int iteration);
[[nodiscard]] __host__ __device__ auto nt_ionisation_ratecoeff(int nonemptymgi, int element, int ion) -> double;
[[nodiscard]] __host__ __device__ auto nt_ionisation_upperion_probability(int nonemptymgi, int element, int lowerion,
                                                                          int upperion, bool energyweighted) -> double;
[[nodiscard]] __host__ __device__ auto nt_ionisation_maxupperion(int element, int lowerion) -> int;
[[nodiscard]] __host__ __device__ auto nt_random_upperion(int nonemptymgi, int element, int lowerion,
                                                          bool energyweighted) -> int;
void calculate_deposition_rate_density(int nonemptymgi, int timestep, HeatingCoolingRates& heatingcoolingrates);
[[nodiscard]] __host__ __device__ auto get_deposition_rate_density(int nonemptymgi) -> double;
[[nodiscard]] auto get_nt_frac_heating(int nonemptymgi) -> float;

[[gnu::pure]] [[nodiscard]] __host__ __device__ auto nt_excitation_ratecoeff(int nonemptymgi, int lowerlevel,
                                                                             int upperlevel, int alltransindex)
    -> double;
__host__ __device__ void do_ntalpha_fisprod_deposit(Packet& pkt);
__host__ __device__ void do_ntlepton_deposit(Packet& pkt);
void write_restart_data(FILE* gridsave_file);
void read_restart_data(FILE* gridsave_file);
void nt_MPI_Bcast(ptrdiff_t nonemptymgi, int root_node_id);
void reset_stats();
void print_stats(double modelvolume, double deltat);
}  // namespace nonthermal

#endif  // NONTHERMAL_H
