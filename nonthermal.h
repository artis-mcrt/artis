#ifndef NONTHERMAL_H
#define NONTHERMAL_H

// #include "sn3d.h"

void nonthermal_init(void);
void nonthermal_close_file(void);
void nt_solve_spencerfano(int modelgridindex, int timestep, double deltaV);
double nt_ionization_ratecoeff(int modelgridindex, int element, int ion);

#endif //NONTHERMAL_H
