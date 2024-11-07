#include <cstddef>
#ifndef STATS_H
#define STATS_H

#include "packet.h"

namespace stats {
// number of ion stats counters that should be divided by the ion populations
constexpr int nstatcounters_ratecoeff = 18;

// one counter per ion per cell
enum ionstattypes {
  ION_RADRECOMB_MACROATOM = 0,
  ION_RADRECOMB_KPKT = 1,
  ION_RADRECOMB_ABSORBED = 2,
  ION_BOUNDBOUND_MACROATOM = 3,
  ION_BOUNDBOUND_ABSORBED = 4,
  ION_NTION = 5,
  ION_PHOTOION = 6,
  ION_PHOTOION_FROMBOUNDFREE = 7,
  ION_PHOTOION_FROMBFSAMEELEMENT = 8,
  ION_PHOTOION_FROMBFIONPLUSONE = 9,
  ION_PHOTOION_FROMBFIONPLUSTWO = 10,
  ION_PHOTOION_FROMBFIONPLUSTHREE = 11,
  ION_PHOTOION_FROMBFLOWERSUPERLEVEL = 12,
  ION_PHOTOION_FROMBOUNDBOUND = 13,
  ION_PHOTOION_FROMBOUNDBOUNDIONPLUSONE = 14,
  ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTWO = 15,
  ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTHREE = 16,
  ION_MACROATOM_ENERGYOUT_RADDEEXC = 17,
  ION_MACROATOM_ENERGYOUT_RADRECOMB = 18,
  ION_MACROATOM_ENERGYOUT_COLLDEEXC = 19,
  ION_MACROATOM_ENERGYOUT_COLLRECOMB = 20,
  ION_MACROATOM_ENERGYIN_RADEXC = 21,
  ION_MACROATOM_ENERGYIN_PHOTOION = 22,
  ION_MACROATOM_ENERGYIN_COLLEXC = 23,
  ION_MACROATOM_ENERGYIN_COLLION = 24,
  ION_MACROATOM_ENERGYIN_NTCOLLION = 26,
  ION_MACROATOM_ENERGYIN_TOTAL = 27,
  ION_MACROATOM_ENERGYOUT_TOTAL = 28,
  ION_MACROATOM_ENERGYIN_INTERNAL = 29,
  ION_MACROATOM_ENERGYOUT_INTERNAL = 30,
  ION_STAT_COUNT = 31,
};

// global statistics (all cells combined)
enum eventcounters {
  COUNTER_MA_STAT_ACTIVATION_COLLEXC = 0,
  COUNTER_MA_STAT_ACTIVATION_COLLION = 1,
  COUNTER_MA_STAT_ACTIVATION_NTCOLLEXC = 2,
  COUNTER_MA_STAT_ACTIVATION_NTCOLLION = 3,
  COUNTER_MA_STAT_ACTIVATION_BB = 4,
  COUNTER_MA_STAT_ACTIVATION_BF = 5,
  COUNTER_MA_STAT_ACTIVATION_FB = 6,
  COUNTER_MA_STAT_DEACTIVATION_COLLDEEXC = 7,
  COUNTER_MA_STAT_DEACTIVATION_COLLRECOMB = 8,
  COUNTER_MA_STAT_DEACTIVATION_BB = 9,
  COUNTER_MA_STAT_DEACTIVATION_FB = 10,
  COUNTER_MA_STAT_INTERNALUPHIGHER = 11,
  COUNTER_MA_STAT_INTERNALUPHIGHERNT = 12,
  COUNTER_MA_STAT_INTERNALDOWNLOWER = 13,
  COUNTER_K_STAT_TO_MA_COLLEXC = 14,
  COUNTER_K_STAT_TO_MA_COLLION = 15,
  COUNTER_K_STAT_TO_R_FF = 16,
  COUNTER_K_STAT_TO_R_FB = 17,
  COUNTER_K_STAT_TO_R_BB = 18,
  COUNTER_K_STAT_FROM_FF = 19,
  COUNTER_K_STAT_FROM_BF = 20,
  COUNTER_NT_STAT_FROM_GAMMA = 21,
  COUNTER_NT_STAT_TO_IONIZATION = 22,
  COUNTER_NT_STAT_TO_EXCITATION = 23,
  COUNTER_NT_STAT_TO_KPKT = 24,
  COUNTER_K_STAT_FROM_EARLIERDECAY = 25,
  COUNTER_INTERACTIONS = 26,
  COUNTER_ESCOUNTER = 27,
  COUNTER_RESONANCESCATTERINGS = 28,
  COUNTER_CELLCROSSINGS = 29,
  COUNTER_UPSCATTER = 30,
  COUNTER_DOWNSCATTER = 31,
  COUNTER_UPDATECELL = 32,
  COUNTER_COUNT = 33,
};

void init();

void increment_ion_stats(int nonemptymgi, int element, int ion, enum ionstattypes ionstattype, double increment);

void increment_ion_stats_contabsorption(const Packet &pkt, int nonemptymgi, int element, int ion);

[[nodiscard]] auto get_ion_stats(int nonemptymgi, int element, int ion, enum ionstattypes ionstattype) -> double;

void set_ion_stats(int nonemptymgi, int element, int ion, enum ionstattypes ionstattype, double newvalue);

void reset_ion_stats(int nonemptymgi);

void normalise_ion_estimators(int nonemptymgi, double deltat, double deltaV);

void increment(enum eventcounters);

void pkt_action_counters_reset();

[[nodiscard]] auto get_counter(enum eventcounters i) -> ptrdiff_t;

void pkt_action_counters_printout(int nts);

void reduce_estimators();
}  // namespace stats

#endif  // STATS_H
