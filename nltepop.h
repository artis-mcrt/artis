#ifndef NLTEPOP_H
#define NLTEPOP_H

#include <stdio.h>

void solve_nlte_pops_element(const int element, const int modelgridindex, const int timestep);
double solve_nlte_pops(const int element, const int ion, const int modelgridindex, const int timestep);
double superlevel_boltzmann(const int modelgridindex, const int element, const int ion, const int level);
void nltepop_write_to_file(const int n, const int timestep);
void nltepop_open_file(const int my_rank);
void nltepop_close_file(void);
void nltepop_write_restart_data(FILE *restart_file);
void nltepop_read_restart_data(FILE *restart_file);

#endif //NLTEPOP_H
