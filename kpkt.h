// Declarations for k-packet (thermal energy packet) handling (kpkt.cc).

#ifndef KPKT_H
#define KPKT_H

#include <cstddef>
#include <span>

#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "mpi_logging.h"
#include "packet.h"
#include "thermalbalance.h"

namespace kpkt {

inline MPI_shared_array<double> ion_cooling_contribs_allcells{};

// the share of the thermal energy flux of each cell that becomes radiation (see
// NLTE_TIME_DEPENDENT_FIRST_TIMESTEP). The array exists only when that option has a value.
inline MPI_shared_array<float> radiative_energy_factor_allcells{};
inline int ncoolingterms{0};

void setup_coolinglist();
void calculate_cooling_rates(int nonemptymgi, HeatingCoolingRates* heatingcoolingrates);
void set_radiative_energy_factor(int nonemptymgi, const HeatingCoolingRates& heatingcoolingrates);
DEVICE_FUNC void do_kpkt_blackbody(Packet& pkt);
DEVICE_FUNC void do_kpkt(Packet& pkt, double t2, int nts);

// prepopulate one ion's cooling-rate contributions into the cellcache. Used in GPU mode, where the lazy
// mutex-guarded calculation in do_kpkt() cannot run safely on the device.
DEVICE_FUNC void calculate_cellcache_cooling_rates_ion(int nonemptymgi, int uniqueionindex);

[[nodiscard]] inline auto get_coolinglistoffset(int element, int ion) -> int {
  return globals::elements[element].ions[ion].coolingoffset;
}

[[nodiscard]] inline auto get_ncoolingterms_ion(int element, int ion) -> int {
  return globals::elements[element].ions[ion].ncoolingterms;
}

// get an array of per-ion cumulative cooling contributions for a given cell
[[nodiscard]] inline auto get_cell_ion_cooling_contribs(const std::ptrdiff_t nonemptymgi) -> std::span<double> {
  return kpkt::ion_cooling_contribs_allcells.subspan(nonemptymgi * get_includedions(), get_includedions());
}
}  // namespace kpkt

#endif  // KPKT_H
