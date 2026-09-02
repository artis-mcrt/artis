// Declarations for non-thermal energy deposition (nonthermal.cc).

#ifndef NONTHERMAL_H
#define NONTHERMAL_H

#include <cstddef>
#include <cstdio>
#include <span>

#include "constants.h"
#include "packet.h"
#include "random.h"
#include "thermalbalance.h"

namespace decay {
struct AnaEmissionPowerPerMass;
}  // namespace decay

namespace nonthermal {
// number of energy points in the Spencer-Fano solution vector
inline constexpr int SFPTS = 4096;

// Solve U * x = b for x, where U is a compacted upper triangular matrix: only the upper triangle is stored, with
// element (i, j) at flattened index (SFPTS * i) - (i * (i + 1) / 2) + j. Externally visible so that unittests.cc can
// check the solution.
void solve_upper_triangular(std::span<const double> uppertri, std::span<const double, SFPTS> bvec,
                            std::span<double, SFPTS> xvec);

void init();
// Assemble and solve the cell's discretised Spencer-Fano degradation balance, including Coulomb heating, excitation,
// ionisation, secondary electrons, and Auger electrons. Convert the electron flux into deposition fractions and
// non-thermal rate coefficients, reusing a recent solution while the ionisation state remains sufficiently similar.
// Shingles et al. (2020), Section 2.5, doi:10.1093/mnras/stz3412.
// Return true when the solver makes a new solution for the cell. A kept solution, and the fallback
// values of a cell without deposition, give false.
auto solve_spencerfano(int nonemptymgi, int timestep, int iteration) -> bool;
[[nodiscard]] DEVICE_FUNC auto nt_ionisation_ratecoeff(int nonemptymgi, int element, int ion) -> double;
[[nodiscard]] DEVICE_FUNC auto nt_ionisation_upperion_probability(int nonemptymgi, int element, int lowerion,
                                                                  int upperion, bool energyweighted) -> double;
[[nodiscard]] DEVICE_FUNC auto nt_ionisation_maxupperion(int element, int lowerion) -> int;
[[nodiscard]] DEVICE_FUNC auto nt_random_upperion(int nonemptymgi, int element, int lowerion, bool energyweighted,
                                                  rngstate_type& rngstate) -> int;
void calculate_deposition_rate_density(int nonemptymgi, HeatingCoolingRates& heatingcoolingrates,
                                       const decay::AnaEmissionPowerPerMass& emission_power_per_mass);
// deposition rate density [erg/s/cm3] of the non-thermal leptons handled by the Spencer-Fano solver, i.e.
// gamma + positron + electron. Alpha and spontaneous fission deposition are tracked separately.
// This returns a stored value, so calculate_deposition_rate_density() must already have run for this cell in
// this timestep; otherwise what comes back is the previous timestep's number.
[[nodiscard]] DEVICE_FUNC auto get_ntlepton_deposition_rate_density(int nonemptymgi) -> double;
[[nodiscard]] auto get_nt_frac_heating(int nonemptymgi) -> float;

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto nt_excitation_ratecoeff(int nonemptymgi, int lowerlevel, int upperlevel,
                                                                     int alltransindex) -> double;
DEVICE_FUNC void do_ntalpha_fisprod_deposit(Packet& pkt);
DEVICE_FUNC void do_ntlepton_deposit(Packet& pkt);
void write_restart_data(FILE* gridsave_file);
void read_restart_data(FILE* gridsave_file);
void nt_MPI_Bcast(ptrdiff_t nstart_nonempty, ptrdiff_t ndo_nonempty, int root_node_id);
void reset_stats();
void print_stats(double modelvolume, double deltat);
}  // namespace nonthermal

#endif  // NONTHERMAL_H
