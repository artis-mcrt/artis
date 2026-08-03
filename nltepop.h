// Declarations for the NLTE level population solver (nltepop.cc).

#ifndef NLTEPOP_H
#define NLTEPOP_H

#include <cstddef>
#include <cstdio>
#include <optional>
#include <span>

#include "constants.h"
#include "mpi_logging.h"

inline MPI_shared_array<double> nltepops_allcells;

void solve_nlte_pops_element(int element, int nonemptymgi, int timestep, int nlte_iter);
// GTH solve for the stationary distribution of the NLTE rate matrix, exposed here so that unittests.cc can test
// it (see the definition in nltepop.cc for the full contract)
[[nodiscard]] auto gth_stationary_distribution(std::span<double> rate_matrix, std::span<double> vec_x)
    -> std::optional<std::ptrdiff_t>;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto superlevel_boltzmann(int nonemptymgi, int element, int ion, int level)
    -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nlte_levelpop_over_rho(int nonemptymgi, int element, int ion,
                                                                        int level) -> double;
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nlte_superlevelpop_over_rho_over_slpartfunc(int nonemptymgi,
                                                                                             int element, int ion)
    -> double;
void nltepop_write_to_file(int nonemptymgi, int timestep);
void nltepop_open_file();
void nltepop_write_restart_data(FILE* restart_file);
void nltepop_read_restart_data(FILE* restart_file);

#endif  // NLTEPOP_H
