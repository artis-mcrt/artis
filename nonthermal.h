#ifndef NONTHERMAL_H
#define NONTHERMAL_H

// #include "sn3d.h"
#include <stdio.h>

void nt_init(int my_rank);
void nt_close_file(void);
void nt_solve_spencerfano(int modelgridindex, int timestep);
double nt_ionization_ratecoeff(int modelgridindex, int element, int ion);
void calculate_deposition_rate_density(const int modelgridindex, const int timestep);
double get_deposition_rate_density(const int modelgridindex);
float get_nt_frac_heating(int modelgridindex);
double calculate_nt_excitation_rate(int modelgridindex, int element, int ion, int lowerlevel, int upperlevel);
void nt_write_restart_data(FILE *gridsave_file);
void nt_read_restart_data(FILE *gridsave_file);
void nt_MPI_Bcast(const int my_rank, const int root, const int root_nstart, const int root_ndo);

#endif //NONTHERMAL_H
