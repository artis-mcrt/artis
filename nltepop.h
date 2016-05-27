#ifndef NLTEPOP_H
#define NLTEPOP_H

void nlte_pops_element(int element, int modelgridindex, int timestep);
double nlte_pops(int element, int ion, int modelgridindex, int timestep);
void read_binding_energies(void);
double nt_ionization_rate(int modelgridindex, int element, int ion);
double superlevel_boltzmann(int modelgridindex, int element, int ion, int level);

#endif //NLTEPOP_H
