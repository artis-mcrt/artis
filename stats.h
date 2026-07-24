// Monte Carlo event counters (stats.cc): the Counter enum and functions to increment, read,
// and print the per-timestep event counts.

#ifndef STATS_H
#define STATS_H

#include <cstddef>

#include "constants.h"

namespace stats {

// global statistics (all cells combined)
enum class Counter {
  MA_STAT_ACTIVATION_COLLEXC = 0,
  MA_STAT_ACTIVATION_COLLION = 1,
  MA_STAT_ACTIVATION_NTCOLLEXC = 2,
  MA_STAT_ACTIVATION_NTCOLLION = 3,
  MA_STAT_ACTIVATION_BB = 4,
  MA_STAT_ACTIVATION_BF = 5,
  MA_STAT_ACTIVATION_FB = 6,
  MA_STAT_DEACTIVATION_COLLDEEXC = 7,
  MA_STAT_DEACTIVATION_COLLRECOMB = 8,
  MA_STAT_DEACTIVATION_BB = 9,
  MA_STAT_DEACTIVATION_FB = 10,
  MA_STAT_INTERNALUPHIGHER = 11,
  MA_STAT_INTERNALUPHIGHERNT = 12,
  MA_STAT_INTERNALDOWNLOWER = 13,
  K_STAT_TO_MA_COLLEXC = 14,
  K_STAT_TO_MA_COLLION = 15,
  K_STAT_TO_R_FF = 16,
  K_STAT_TO_R_FB = 17,
  K_STAT_TO_R_BB = 18,
  K_STAT_FROM_FF = 19,
  K_STAT_FROM_BF = 20,
  NT_STAT_FROM_GAMMA = 21,
  NT_STAT_TO_IONISATION = 22,
  NT_STAT_TO_EXCITATION = 23,
  NT_STAT_TO_KPKT = 24,
  K_STAT_FROM_EARLIERDECAY = 25,
  INTERACTIONS = 26,
  ELECTRON_SCATTERINGS = 27,
  RESONANCESCATTERINGS = 28,
  CELLCROSSINGS = 29,
  UPSCATTER = 30,
  DOWNSCATTER = 31,
  UPDATECELL = 32,
  PKTESCAPES = 33,
  COUNT = 34,
};

DEVICE_FUNC void increment(Counter i);

void pkt_action_counters_reset();

[[nodiscard]] auto get_counter(Counter i) -> ptrdiff_t;

void pkt_action_counters_printout(int nts);

}  // namespace stats

#endif  // STATS_H
