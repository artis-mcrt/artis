// Declarations for the energy source of THERMAL_ENERGY_INPUT_FROM_FILES (energyinput.cc).

#ifndef ENERGYINPUT_H
#define ENERGYINPUT_H

#include <utility>

#include "packet.h"

namespace energyinput {

// Read energyrate.txt: the cumulative fraction of the cell energies that the cells release by each time.
void read_release_table();

// The cumulative fraction of the cell energies that the cells release by the time t [s].
[[nodiscard]] auto get_cumulative_fraction_released(double t) -> double;

// The fraction of the cell energies that the cells release before the time t [s] and that remains at t. A release
// at the time t' keeps t' / t of its energy, as for a decay.
[[nodiscard]] auto get_released_fraction_after_expansion(double t) -> double;

// The release rate of the cell energy over the timestep nts, per volume at the midpoint [erg/s/cm^3].
[[nodiscard]] auto get_release_rate_density(int nonemptymgi, int nts) -> double;

// The times [s] between which the pellets release the cell energy: from tmin, or with INITIAL_PACKETS_ON from 0, to
// tmax. With INITIAL_PACKETS_ON, the release before tmin seeds the cells at tmin.
[[nodiscard]] auto get_pellet_release_time_range() -> std::pair<double, double>;

// Give the pellet pkt the energy e_cmf_per_packet and a release time that samples the cumulative fraction in the
// range of get_pellet_release_time_range().
void setup_thermal_input_pellet(double e_cmf_per_packet, Packet& pkt);

}  // namespace energyinput

#endif  // ENERGYINPUT_H
