#ifndef KPKT_H
#define KPKT_H

#include <cstddef>
#include <span>

#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "packet.h"
#include "thermalbalance.h"

constexpr double COOLING_UNDEFINED = -99;

namespace kpkt {

inline std::span<double> ion_cooling_contribs_allcells{};
inline int ncoolingterms{0};

auto get_ncoolingterms() -> int;
void setup_coolinglist();
void set_kpktdiffusion(float kpktdiffusion_timescale_in, int n_kpktdiffusion_timesteps_in);
void calculate_cooling_rates(int nonemptymgi, HeatingCoolingRates* heatingcoolingrates);
__host__ __device__ void do_kpkt_blackbody(Packet& pkt);
__host__ __device__ void do_kpkt(Packet& pkt, double t2, int nts);

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
