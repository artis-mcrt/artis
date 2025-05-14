#ifndef MACROATOM_H
#define MACROATOM_H

#include "globals.h"
#include "packet.h"

void macroatom_open_file(int my_rank);
void macroatom_close_file();

void do_macroatom(Packet &pkt, const MacroAtomState &pktmastate);

[[nodiscard]] auto rad_deexcitation_ratecoeff(int nonemptymgi, int lower_uniquelevelindex, double epsilon_trans,
                                              float A_ul, double upperstatweight, double lowerstatweight,
                                              double nnlevelupper, double t_current) -> double;

[[nodiscard]] auto rad_excitation_ratecoeff(int nonemptymgi, int upper_uniquelevelindex, double upper_statweight,
                                            const LevelTransition &uptrans, double epsilon_trans, double nnlevel_lower,
                                            double statweight_lower, int lineindex, double t_current) -> double;

[[nodiscard]] auto rad_recombination_ratecoeff(float T_e, float nne, int element, int upperion, int upperionlevel,
                                               int lowerionlevel, int nonemptymgi) -> double;

[[nodiscard]] auto stim_recombination_ratecoeff(float nne, int element, int upperion, int upper, int lower,
                                                int nonemptymgi) -> double;

[[nodiscard]] auto col_recombination_ratecoeff(float T_e, float nne, int element, int upperion, int upper, int lower,
                                               double epsilon_trans) -> double;

[[nodiscard]] auto col_ionization_ratecoeff(float T_e, float nne, int element, int ion, int lower, int phixstargetindex,
                                            double epsilon_trans) -> double;

[[nodiscard]] auto col_deexcitation_ratecoeff(float T_e, float nne, double epsilon_trans, double upperstatweight,
                                              double lowerstatweight, float coll_str, float osc_strength,
                                              bool forbidden) -> double;

[[nodiscard]] auto col_excitation_ratecoeff(float T_e, float nne, double upperstatweight,
                                            const LevelTransition &uptrans, double epsilon_trans,
                                            double lowerstatweight) -> double;

#endif  // MACROATOM_H
