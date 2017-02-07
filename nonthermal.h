#ifndef NONTHERMAL_H
#define NONTHERMAL_H

// #include "sn3d.h"

void nonthermal_init(int my_rank);
void nonthermal_close_file(void);
void nt_solve_spencerfano(int modelgridindex, int timestep);
double nt_ionization_ratecoeff(int modelgridindex, int element, int ion);
double get_deposition_rate_density(int modelgridindex);
float get_nt_frac_heating(int modelgridindex);
double calculate_nt_excitation_rate(int modelgridindex, int element, int ion, int lowerlevel, int upperlevel);
void nonthermal_MPI_Bcast(int root, int my_rank, int nstart, int ndo);


#endif //NONTHERMAL_H
