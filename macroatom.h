#ifndef MACROATOM_H
#define MACROATOM_H

#include "grid_init.h"
#include "sn3d.h"
#include "types.h"

double do_ma(PKT *pkt_ptr, double t1, double t2, int timestep);

double rad_deexcitation(int modelgridindex, int lower, double epsilon_trans,
                        int lineindex, double t_current);
double rad_excitation(int modelgridindex, int upper, double epsilon_trans,
                      int lineindex, double t_current);
double rad_recombination(int modelgridindex, int lower);
double photoionization(int modelgridindex, int phixstargetindex,
                       double epsilon_trans);

double col_deexcitation(int modelgridindex, int lower, double epsilon_trans,
                        int lineindex);
double col_excitation(int modelgridindex, int upper, int lineindex,
                      double epsilon_trans);
double col_recombination(int modelgridindex, int lower, double epsilon_trans);
double col_ionization(int modelgridindex, int phixstargetindex,
                      double epsilon_trans);

static inline
double get_individ_rad_deexc(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].
         individ_rad_deexc[i];
}

static inline
double get_individ_internal_down_same(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].
         individ_internal_down_same[i];
}

static inline
double get_individ_internal_up_same(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].
         individ_internal_up_same[i];
}

#endif //MACROATOM_H
