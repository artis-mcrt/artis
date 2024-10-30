#ifndef NLTEPOP_H
#define NLTEPOP_H

#include <cstdio>

void solve_nlte_pops_element(int element, int modelgridindex, int timestep, int nlte_iter);
[[nodiscard]] auto superlevel_boltzmann(int modelgridindex, int element, int ion, int level) -> double;
[[nodiscard]] auto get_nlte_levelpop_over_rho(int modelgridindex, int element, int ion, int level) -> double;
[[nodiscard]] auto get_nlte_superlevelpop_over_rho(int modelgridindex, int element, int ion) -> double;
void set_nlte_levelpop_over_rho(int modelgridindex, int element, int ion, int level, double value);
void set_nlte_superlevelpop_over_rho(int modelgridindex, int element, int ion, double value);
void nltepop_write_to_file(int modelgridindex, int timestep);
void nltepop_open_file(int my_rank);
void nltepop_close_file();
void nltepop_write_restart_data(FILE *restart_file);
void nltepop_read_restart_data(FILE *restart_file);

#endif  // NLTEPOP_H
