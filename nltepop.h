#ifndef NLTEPOP_H
#define NLTEPOP_H

#include <cstdio>

void solve_nlte_pops_element(int element, int modelgridindex, int timestep, int nlte_iter);
double solve_nlte_pops_ion(int element, int ion, int modelgridindex, int timestep);
double superlevel_boltzmann(int modelgridindex, int element, int ion, int level);
void nltepop_write_to_file(int n, int timestep);
void nltepop_open_file(int my_rank);
void nltepop_close_file(void);
void nltepop_write_restart_data(FILE *restart_file);
void nltepop_read_restart_data(FILE *restart_file);

#endif //NLTEPOP_H
