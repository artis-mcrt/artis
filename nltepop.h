// Declarations for the NLTE level population solver (nltepop.cc).

#ifndef NLTEPOP_H
#define NLTEPOP_H

#include <cstddef>
#include <cstdio>
#include <optional>
#include <span>
#include <utility>

#include "constants.h"
#include "mpi_logging.h"

inline MPI_shared_array<double> nltepops_allcells;

void solve_nlte_pops_element(int element, int nonemptymgi, int timestep, int nlte_iter);
// the ion range of the last NLTE matrix solution of an element in a cell is the lowermost and the
// uppermost ion of the cell (see grid.h). A lowermost ion of -1 means that the cell holds no NLTE
// solution for the element. The charge transfer reactions read the range, and the NLTE iteration
// loop watches it for a change.
void nltepop_reset_solution_ranges(int nonemptymgi);
[[nodiscard]] auto get_nlte_solution_range(int nonemptymgi, int element) -> std::pair<int, int>;
[[nodiscard]] auto ion_in_nlte_solution(int nonemptymgi, int element, int ion) -> bool;
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
