#ifndef NLTEPOP_H
#define NLTEPOP_H

void nlte_pops_element(int element, int modelgridindex, int timestep);
double nlte_pops(int element, int ion, int modelgridindex, int timestep);
double superlevel_boltzmann(int modelgridindex, int element, int ion, int level);
void write_to_nlte_file(int n, int timestep);

#endif //NLTEPOP_H
