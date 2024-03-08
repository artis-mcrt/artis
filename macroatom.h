#pragma once
#ifndef MACROATOM_H
#define MACROATOM_H

#include <cmath>

enum ma_action {
  /// Radiative deexcitation rate from this level.
  MA_ACTION_RADDEEXC = 0,
  /// Collisional deexcitation rate from this level.
  MA_ACTION_COLDEEXC = 1,
  /// Radiative recombination from this level.
  MA_ACTION_RADRECOMB = 2,
  /// Collisional recombination rate from this level.
  MA_ACTION_COLRECOMB = 3,
  /// Rate for internal downward transitions to same ionisation stage.
  MA_ACTION_INTERNALDOWNSAME = 4,
  /// Rate for internal upward transitions to same ionisation stage.
  MA_ACTION_INTERNALDOWNLOWER = 5,
  /// Rate for internal downward transitions to lower ionisation stage.
  MA_ACTION_INTERNALUPSAME = 6,
  /// Rate for internal upward transitions to higher ionisation stage.
  MA_ACTION_INTERNALUPHIGHER = 7,
  /// Rate for internal upward transitions to higher ionisation stage due to non-thermal collisions.
  MA_ACTION_INTERNALUPHIGHERNT = 8,
  MA_ACTION_COUNT = 9,
};

#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "packet.h"

void macroatom_open_file(int my_rank);
void macroatom_close_file();

void do_macroatom(Packet &pkt_ptr, const MacroAtomState &pktmastate);

[[nodiscard]] auto rad_deexcitation_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower,
                                              double epsilon_trans, float A_ul, double upperstatweight,
                                              double t_current) -> double;
[[nodiscard]] auto rad_excitation_ratecoeff(int modelgridindex, int element, int ion, int lower, int uptransindex,
                                            double epsilon_trans, int lineindex, double t_current) -> double;
[[nodiscard]] auto rad_recombination_ratecoeff(float T_e, float nne, int element, int upperion, int upperionlevel,
                                               int lowerionlevel, int modelgridindex) -> double;
[[nodiscard]] auto stim_recombination_ratecoeff(float nne, int element, int upperion, int upper, int lower,
                                                int modelgridindex) -> double;

[[nodiscard]] auto col_recombination_ratecoeff(int modelgridindex, int element, int upperion, int upper, int lower,
                                               double epsilon_trans) -> double;
[[nodiscard]] auto col_ionization_ratecoeff(float T_e, float nne, int element, int ion, int lower, int phixstargetindex,
                                            double epsilon_trans) -> double;

[[nodiscard]] auto col_deexcitation_ratecoeff(float T_e, float nne, double epsilon_trans, int element, int ion,
                                              int upper, const LevelTransition &downtransition) -> double;

[[nodiscard]] auto col_excitation_ratecoeff(float T_e, float nne, int element, int ion, int lower, int uptransindex,
                                            double epsilon_trans, double lowerstatweight) -> double;

#endif  // MACROATOM_H
