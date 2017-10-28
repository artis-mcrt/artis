#ifndef NONTHERMAL_H
#define NONTHERMAL_H

// #include "sn3d.h"
#include "types.h"
#include <stdio.h>

void nt_init(int my_rank);
void nt_close_file(void);
void nt_solve_spencerfano(int modelgridindex, int timestep);
double nt_ionization_ratecoeff(int modelgridindex, int element, int ion);
void calculate_deposition_rate_density(int modelgridindex, int timestep);
double get_deposition_rate_density(int modelgridindex);
float get_nt_frac_heating(int modelgridindex);
float get_nt_frac_excitation(const int modelgridindex);
float get_nt_frac_ionization(const int modelgridindex);
double calculate_nt_excitation_rate(int modelgridindex, int element, int ion, int lowerlevel, int upperlevel);
void do_nt_electron(PKT *pkt_ptr);
void nt_write_restart_data(FILE *gridsave_file);
void nt_read_restart_data(FILE *gridsave_file);
void nt_MPI_Bcast(int my_rank, int root, int root_nstart, int root_ndo);

#endif //NONTHERMAL_H
