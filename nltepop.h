#ifndef NLTEPOP_H
#define NLTEPOP_H

#include <cstdio>

#include "constants.h"

void solve_nlte_pops_element(int element, int nonemptymgi, int timestep, int nlte_iter);
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto superlevel_boltzmann(int nonemptymgi, int element, int ion, int level)
    -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nlte_levelpop_over_rho(int nonemptymgi, int element, int ion,
                                                                        int level) -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nlte_superlevelpop_over_rho_over_slpartfunc(int nonemptymgi,
                                                                                             int element, int ion)
    -> double;
void nltepop_write_to_file(int nonemptymgi, int timestep);
void nltepop_open_file(int my_rank);
void nltepop_write_restart_data(FILE* restart_file);
void nltepop_read_restart_data(FILE* restart_file);

#endif  // NLTEPOP_H
