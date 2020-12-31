#ifndef STATS_H
#define STATS_H

#include "types.h"

namespace stats {
  // number of ion stats counters that should be divided by the ion populations
  #define nstatcounters_ratecoeff 18

  // one counter per ion per cell
  enum ionstattypes {
    ION_RADRECOMB_MACROATOM = 0,
    ION_RADRECOMB_KPKT = 1,
    ION_RADRECOMB_ABSORBED = 2,
    ION_RADRECOMB_ESCAPED = 3,
    ION_BOUNDBOUND_MACROATOM = 4,
    ION_BOUNDBOUND_ABSORBED = 5,
    ION_NTION = 6,
    ION_PHOTOION = 7,
    ION_PHOTOION_FROMBOUNDFREE = 8,
    ION_PHOTOION_FROMBFSAMEELEMENT = 9,
    ION_PHOTOION_FROMBFIONPLUSONE = 10,
    ION_PHOTOION_FROMBFIONPLUSTWO = 11,
    ION_PHOTOION_FROMBFIONPLUSTHREE = 12,
    ION_PHOTOION_FROMBFLOWERSUPERLEVEL = 13,
    ION_PHOTOION_FROMBOUNDBOUND = 14,
    ION_PHOTOION_FROMBOUNDBOUNDIONPLUSONE = 15,
    ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTWO = 16,
    ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTHREE = 17,
    ION_MACROATOM_ENERGYOUT_RADDEEXC = 18,
    ION_MACROATOM_ENERGYOUT_RADRECOMB = 19,
    ION_MACROATOM_ENERGYOUT_COLLDEEXC = 20,
    ION_MACROATOM_ENERGYOUT_COLLRECOMB = 21,
    ION_MACROATOM_ENERGYIN_RADEXC = 22,
    ION_MACROATOM_ENERGYIN_PHOTOION = 23,
    ION_MACROATOM_ENERGYIN_COLLEXC = 24,
    ION_MACROATOM_ENERGYIN_COLLION = 25,
    ION_MACROATOM_ENERGYIN_NTCOLLION = 27,
    ION_MACROATOM_ENERGYIN_TOTAL = 28,
    ION_MACROATOM_ENERGYOUT_TOTAL = 29,
    ION_MACROATOM_ENERGYIN_INTERNAL = 30,
    ION_MACROATOM_ENERGYOUT_INTERNAL = 31,
    ION_STAT_COUNT = 32,
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
    COUNTER_ESCOUNTER = 26,
    COUNTER_RESONANCESCATTERINGS = 27,
    COUNTER_CELLCROSSINGS = 28,
    COUNTER_UPSCATTER = 29,
    COUNTER_DOWNSCATTER = 30,
    COUNTER_UPDATECELL = 31,
    COUNTER_COOLINGRATECALCCOUNTER = 32,
    COUNTER_NESC = 33,
    COUNTER_COUNT = 34,
  };

  void init(void);

  void cleanup(void);

  void increment_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstattypes ionstattype, const double increment);

  void increment_ion_stats_contabsorption(const PKT *const pkt_ptr, const int modelgridindex, const int element, const int ion);

  double get_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstattypes ionstattype);

  void set_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstattypes ionstattype, const double newvalue);

  void reset_ion_stats(int modelgridindex);

  void normalise_ion_estimators(const int mgi, const double deltat, const double deltaV);

  void increment(enum eventcounters);

  void pkt_action_counters_reset(void);

  int get_counter(enum eventcounters i);

  void pkt_action_counters_printout(const PKT *const pkt, const int nts);
}

#endif //STATS_H
