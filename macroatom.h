#ifndef MACROATOM_H
#define MACROATOM_H

#include "packet.h"

void macroatom_open_file(int my_rank);
void macroatom_close_file();

void do_macroatom(Packet& pkt, const MacroAtomState& pktmastate);

[[gnu::pure]] [[nodiscard]] auto rad_deexcitation_ratecoeff(int nonemptymgi, int lower_uniquelevelindex,
                                                            double epsilon_trans, float A_ul, double upperstatweight,
                                                            double lowerstatweight, double nnlevelupper,
                                                            double t_current) -> double;

[[gnu::pure]] [[nodiscard]] auto rad_excitation_ratecoeff(int nonemptymgi, int upper_uniquelevelindex,
                                                          double upper_statweight, double einstein_A,
                                                          double epsilon_trans, double nnlevel_lower,
                                                          double statweight_lower, int alltransindex, double t_current)
    -> double;

[[gnu::pure]] [[nodiscard]] auto rad_recombination_ratecoeff(float T_e, float nne, int element, int upperion,
                                                             int upperionlevel, int lowerionlevel, int nonemptymgi)
    -> double;

[[gnu::pure]] [[nodiscard]] auto stim_recombination_ratecoeff(float nne, int element, int upperion, int upper,
                                                              int lower, int nonemptymgi) -> double;

[[gnu::pure]] [[nodiscard]] auto col_recombination_ratecoeff(float T_e, float nne, int element, int upperion, int upper,
                                                             int lower, double epsilon_trans) -> double;

[[gnu::pure]] [[nodiscard]] auto col_ionisation_ratecoeff(float T_e, float nne, int element, int ion, int lower,
                                                          int phixstargetindex, double epsilon_trans) -> double;

[[gnu::pure]] [[nodiscard]] auto col_deexcitation_ratecoeff(float T_e, float nne, double epsilon_trans,
                                                            double upperstatweight, double lowerstatweight,
                                                            int alltransindex) -> double;

[[gnu::pure]] [[nodiscard]] auto col_excitation_ratecoeff(float T_e, float nne, double upperstatweight,
                                                          int alltransindex, double epsilon_trans,
                                                          double lowerstatweight) -> double;

#endif  // MACROATOM_H
