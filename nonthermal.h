// Declarations for non-thermal energy deposition (nonthermal.cc).

#ifndef NONTHERMAL_H
#define NONTHERMAL_H

#include <cstddef>
#include <cstdio>

#include "constants.h"
#include "decay.h"
#include "packet.h"
#include "random.h"
#include "thermalbalance.h"

namespace nonthermal {
void init();
// Assemble and solve the cell's discretised Spencer-Fano degradation balance, including Coulomb heating, excitation,
// ionisation, secondary electrons, and Auger electrons. Convert the electron flux into deposition fractions and
// non-thermal rate coefficients, reusing a recent solution while the ionisation state remains sufficiently similar.
// Shingles et al. (2020), Section 2.5, doi:10.1093/mnras/stz3412.
void solve_spencerfano(int nonemptymgi, int timestep, int iteration);
[[nodiscard]] DEVICE_FUNC auto nt_ionisation_ratecoeff(int nonemptymgi, int element, int ion) -> double;
[[nodiscard]] DEVICE_FUNC auto nt_ionisation_upperion_probability(int nonemptymgi, int element, int lowerion,
                                                                  int upperion, bool energyweighted) -> double;
[[nodiscard]] DEVICE_FUNC auto nt_ionisation_maxupperion(int element, int lowerion) -> int;
[[nodiscard]] DEVICE_FUNC auto nt_random_upperion(int nonemptymgi, int element, int lowerion, bool energyweighted,
                                                  rngstate_type& rngstate) -> int;
void calculate_deposition_rate_density(int nonemptymgi, HeatingCoolingRates& heatingcoolingrates,
                                       const decay::AnaEmissionRateCoeffs& emissionratecoeffs);
// deposition rate density [erg/s/cm3] of the non-thermal leptons handled by the Spencer-Fano solver, i.e.
// gamma + positron + electron. Alpha and spontaneous fission deposition are tracked separately.
[[nodiscard]] DEVICE_FUNC auto get_ntlepton_deposition_rate_density(int nonemptymgi) -> double;
[[nodiscard]] auto get_nt_frac_heating(int nonemptymgi) -> float;

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto nt_excitation_ratecoeff(int nonemptymgi, int lowerlevel, int upperlevel,
                                                                     int alltransindex) -> double;
DEVICE_FUNC void do_ntalpha_fisprod_deposit(Packet& pkt);
DEVICE_FUNC void do_ntlepton_deposit(Packet& pkt);
void write_restart_data(FILE* gridsave_file);
void read_restart_data(FILE* gridsave_file);
void nt_MPI_Bcast(ptrdiff_t nonemptymgi, int root_node_id);
void reset_stats();
void print_stats(double modelvolume, double deltat);
}  // namespace nonthermal

#endif  // NONTHERMAL_H
