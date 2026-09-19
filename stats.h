// Monte Carlo event counters (stats.cc): the Counter enum and functions to increment, read,
// and print the per-timestep event counts.

#ifndef STATS_H
#define STATS_H

#include <cstddef>

#include "constants.h"

namespace stats {

// global statistics (all cells combined)
enum class Counter {
  MA_STAT_ACTIVATION_COLLEXC,
  MA_STAT_ACTIVATION_COLLION,
  MA_STAT_ACTIVATION_NTCOLLEXC,
  MA_STAT_ACTIVATION_NTCOLLION,
  MA_STAT_ACTIVATION_BB,
  MA_STAT_ACTIVATION_BF,
  MA_STAT_DEACTIVATION_COLLDEEXC,
  MA_STAT_DEACTIVATION_COLLRECOMB,
  MA_STAT_DEACTIVATION_BB,
  MA_STAT_DEACTIVATION_FB,
  MA_STAT_INTERNALUPHIGHER,
  MA_STAT_INTERNALUPHIGHERNT,
  MA_STAT_INTERNALDOWNLOWER,
  K_STAT_TO_MA_COLLEXC,
  K_STAT_TO_MA_COLLION,
  K_STAT_TO_R_FF,
  K_STAT_TO_R_FB,
  K_STAT_TO_R_BB,
  K_STAT_FROM_FF,
  K_STAT_FROM_BF,
  NT_STAT_FROM_GAMMA,
  NT_STAT_TO_IONISATION,
  NT_STAT_TO_EXCITATION,
  NT_STAT_TO_KPKT,
  K_STAT_FROM_EARLIERDECAY,
  INTERACTIONS,
  ELECTRON_SCATTERINGS,
  RESONANCESCATTERINGS,
  CELLCROSSINGS,
  UPSCATTER,
  DOWNSCATTER,
  UPDATECELL,
  PKTESCAPES,
  // Keep COUNT as the last entry. It gives the number of the counters.
  COUNT,
};

DEVICE_FUNC void increment(Counter i);

void pkt_action_counters_reset();

[[nodiscard]] auto get_counter(Counter i) -> ptrdiff_t;

void pkt_action_counters_printout(int nts);

}  // namespace stats

#endif  // STATS_H
