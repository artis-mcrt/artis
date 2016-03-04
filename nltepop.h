#ifndef NLTEPOP_H
#define NLTEPOP_H

double nlte_pops(int element, int ion, int modelgridindex, int timestep);
double get_tot_nion(int modelgridindex);
double get_oneoverw(int element, int ion, int modelgridindex);
double get_mean_binding_energy(int element, int ion);
int read_binding_energies(void);
double nt_ionization_rate(int modelgridindex, int element, int ion);
double superlevel_boltzmann(int modelgridindex, int element, int ion, int level);

#endif //NLTEPOP_H
