#ifndef NONTHERMAL_H
#define NONTHERMAL_H

#include <cstdio>

#include "packet.h"
#include "input.h"

namespace nonthermal {
void init(int my_rank, int ndo_nonempty);
void close_file();
void solve_spencerfano(int modelgridindex, int timestep, int iteration);
auto nt_ionization_ratecoeff(int modelgridindex, int element, int ion) -> double;
auto nt_ionization_upperion_probability(int modelgridindex, int element, int lowerion, int upperion,
                                        bool energyweighted) -> double;
auto nt_ionisation_maxupperion(int element, int lowerion) -> int;
auto nt_random_upperion(int modelgridindex, int element, int lowerion, bool energyweighted) -> int;
void calculate_deposition_rate_density(int modelgridindex, int timestep);
auto get_deposition_rate_density(int modelgridindex) -> double;
auto get_nt_frac_heating(int modelgridindex) -> float;
auto nt_excitation_ratecoeff(int modelgridindex, int element, int ion, int lowerlevel, int uptransindex,
                             double epsilon_trans, int lineindex) -> double;
void do_ntlepton(struct packet *pkt_ptr);
void write_restart_data(FILE *gridsave_file);
void read_restart_data(FILE *gridsave_file);
void nt_MPI_Bcast(int modelgridindex, int root);
void nt_reset_stats();
void nt_print_stats(double modelvolume, double deltat);
}  // namespace nonthermal

#endif  // NONTHERMAL_H
