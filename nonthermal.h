#ifndef NONTHERMAL_H
#define NONTHERMAL_H

#include <cstdio>

#include "packet.h"

namespace nonthermal {
void init(int my_rank, int ndo, int ndo_nonempty);
void close_file();
void solve_spencerfano(int modelgridindex, int timestep, int iteration);
double nt_ionization_ratecoeff(int modelgridindex, int element, int ion);
double nt_ionization_upperion_probability(int modelgridindex, int element, int lowerion, int upperion,
                                          bool energyweighted);
int nt_ionisation_maxupperion(int element, int lowerion);
int nt_random_upperion(int modelgridindex, int element, int lowerion, bool energyweighted);
void calculate_deposition_rate_density(int modelgridindex, int timestep);
double get_deposition_rate_density(int modelgridindex);
float get_nt_frac_heating(int modelgridindex);
double nt_excitation_ratecoeff(int modelgridindex, int element, int ion, int lower, int upper, double epsilon_trans,
                               int lineindex);
void do_ntlepton(struct packet *pkt_ptr);
void write_restart_data(FILE *gridsave_file);
void read_restart_data(FILE *gridsave_file);
void nt_MPI_Bcast(int modelgridindex, int root);
void nt_reset_stats();
void nt_print_stats(int timestep, double modelvolume, double deltat);
}  // namespace nonthermal

#endif  // NONTHERMAL_H
