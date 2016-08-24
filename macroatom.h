#ifndef MACROATOM_H
#define MACROATOM_H

#include "grid_init.h"
#include "sn3d.h"
#include "types.h"

double do_ma(PKT *restrict pkt_ptr, double t1, double t2, int timestep);

double rad_deexcitation(int modelgridindex, int lower, double epsilon_trans, int lineindex, double t_current);
double rad_deexcitation_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower, double epsilon_trans, int lineindex, double t_current);
double rad_excitation(int modelgridindex, int upper, double epsilon_trans, int lineindex, double t_current);
double rad_excitation_ratecoeff(int modelgridindex, int element, int ion, int lower, int upper, double epsilon_trans, int lineindex, double t_current);
double rad_recombination_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower);
double rad_recombination(int modelgridindex, int lower);

double col_deexcitation(int modelgridindex, int lower, double epsilon_trans, int lineindex);
double col_deexcitation_ratecoeff(int modelgridindex, int upper, int lower, double epsilon_trans, int lineindex);
double col_excitation_ratecoeff(int modelgridindex, int lineindex, double epsilon_trans);
double col_excitation(int modelgridindex, int upper, int lineindex, double epsilon_trans);
double col_recombination_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower, double epsilon_trans);
double col_recombination(int modelgridindex, int lower, double epsilon_trans);
double col_ionization_ratecoeff(int modelgridindex, int element, int ion, int lower, int phixstargetindex, double epsilon_trans);
double col_ionization(int modelgridindex, int phixstargetindex, double epsilon_trans);

inline double get_individ_rad_deexc(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[i];
}

inline double get_individ_internal_down_same(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[i];
}

inline double get_individ_internal_up_same(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[i];
}

#endif //MACROATOM_H
