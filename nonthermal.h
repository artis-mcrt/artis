#ifndef NONTHERMAL_H
#define NONTHERMAL_H

// #include "sn3d.h"

void nonthermal_init(void);
void nonthermal_close_file(void);
void nt_solve_spencerfano(int modelgridindex, int timestep);
double nt_ionization_ratecoeff(int modelgridindex, int element, int ion);
double get_deposition_rate_density(int modelgridindex);
double get_nt_frac_heating(int modelgridindex);
double get_nt_excitation_rate(int modelgridindex, int element, int ion, int lowerlevel, int upperlevel);

#endif //NONTHERMAL_H
