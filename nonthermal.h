#ifndef NONTHERMAL_H
#define NONTHERMAL_H

// #include "sn3d.h"
#include "types.h"
#include <stdio.h>

void nt_init(int my_rank);
void nt_close_file(void);
void nt_solve_spencerfano(int modelgridindex, int timestep, int iteration);
double nt_ionization_ratecoeff(int modelgridindex, int element, int ion);
double nt_ionization_upperion_probability(int modelgridindex, int element, int lowerion, int upperion, bool energyweighted);
int nt_ionisation_maxupperion(int element, int lowerion);
int nt_random_upperion(int modelgridindex, int element, int lowerion, bool energyweighted);
void calculate_deposition_rate_density(int modelgridindex, int timestep);
double get_deposition_rate_density(int modelgridindex);
float get_nt_frac_heating(int modelgridindex);
double nt_excitation_ratecoeff(int modelgridindex, int element, int ion, int lower, int upper, double epsilon_trans, int lineindex);
double do_ntlepton(PKT *pkt_ptr, const double t_current, const double t2, const int timestep);
void nt_write_restart_data(FILE *gridsave_file);
void nt_read_restart_data(FILE *gridsave_file);
void nt_MPI_Bcast(const int modelgridindex, const int root);
void nt_reset_stats(void);
void nt_print_stats(const int timestep, const double modelvolume, const double deltat);

#endif //NONTHERMAL_H
