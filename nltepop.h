// Declarations for the NLTE level population solver (nltepop.cc).

#ifndef NLTEPOP_H
#define NLTEPOP_H

#include <cstddef>
#include <cstdio>
#include <optional>
#include <span>
#include <utility>

#include "artisoptions.h"
#include "constants.h"
#include "mpi_logging.h"

inline MPI_shared_array<double> nltepops_allcells;

// The state of the previous grid update for the time-dependent ionisation and thermal balance (see
// NLTE_TIME_DEPENDENT_FIRST_TIMESTEP). The arrays exist only when that option has a value. Only the rank that
// owns a cell reads and writes the entries of the cell.
inline MPI_shared_array<float> nnion_prev_allcells;  // ion population over the sum of the ion populations
inline MPI_shared_array<float> Te_prev_allcells;
inline MPI_shared_array<double> prev_solution_time_allcells;  // t_mid of the stored previous state, or -1 for none
inline MPI_shared_array<double> solution_time_allcells;  // t_mid of the solution that the grid holds, or -1

static_assert(!NLTE_TIME_DEPENDENT_FIRST_TIMESTEP.has_value() || *NLTE_TIME_DEPENDENT_FIRST_TIMESTEP >= 0);

// The timestep uses the time-dependent ionisation and thermal balance equations
[[nodiscard]] constexpr auto timestep_is_time_dependent(const int nts) -> bool {
  return NLTE_TIME_DEPENDENT_FIRST_TIMESTEP.has_value() && nts >= *NLTE_TIME_DEPENDENT_FIRST_TIMESTEP;
}

void nltepop_allocate_time_dependent_arrays();
void nltepop_store_previous_state(int nonemptymgi);
void nltepop_update_solution_time(int nonemptymgi, int nts);
void nltepop_clear_solution_time(int nonemptymgi);
[[nodiscard]] auto get_time_dependent_dt(int nonemptymgi, int nts) -> std::optional<double>;
[[nodiscard]] auto get_nne_prev(int nonemptymgi) -> double;

void solve_nlte_pops_element(int element, int nonemptymgi, int timestep, int nlte_iter);
// the ion range of the last NLTE matrix solution of an element in a cell is the lowermost and the
// uppermost ion of the cell (see grid.h). The cell holds no NLTE solution for the element before
// the first solve and after a fallback to LTE, and the range is then {-1, -1}. The charge transfer
// reactions read the range, and the NLTE iteration loop watches it for a change.
// nltepop_reset_cell() clears the NLTE solution of every element of the cell. The level
// populations get the -1 marker, and the ranges and the solution time are reset.
void nltepop_reset_cell(int nonemptymgi);
// An ion with a smaller fraction of the element population does not change the map of the outer
// iteration. The ion reduction removes only an ion below the smaller limit. A reduced range
// therefore counts as a change only when it takes a significant ion.
constexpr double NLTE_SIGNIFICANT_ION_FRACTION = 1e-4;
static_assert(NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION < NLTE_SIGNIFICANT_ION_FRACTION);
[[nodiscard]] auto elem_has_nlte_solution(int nonemptymgi, int element) -> bool;
[[nodiscard]] auto get_nlte_solution_range(int nonemptymgi, int element) -> std::pair<int, int>;
// Give each ion of the element the factor of ion_factors, and hold the element population at its
// abundance value with a common factor. The outer iteration uses this to inject the ion populations
// that its accelerator gives (see the definition in nltepop.cc for the full contract).
[[nodiscard]] auto nltepop_scale_ion_pops(int nonemptymgi, int element, std::span<const double> ion_factors) -> bool;
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
