#ifndef MACROATOM_H
#define MACROATOM_H

#include "grid_init.h"
#include "types.h"
#include "sn3d.h"

double do_ma(PKT *pkt_ptr, double t1, double t2, int timestep);

double rad_deexcitation(int modelgridindex, int lower, double epsilon_trans, int lineindex, double t_current);
double rad_excitation(int modelgridindex, int upper, double epsilon_trans, int lineindex, double t_current);//, double T_R, double W);
double rad_recombination(int modelgridindex, int lower, double epsilon_trans);
double photoionization(int modelgridindex, int phixstargetindex, double epsilon_trans);

double col_deexcitation(int modelgridindex, int lower, double epsilon_trans, int lineindex);
double col_excitation(int modelgridindex, int upper, int lineindex, double epsilon_trans);
double col_recombination(int modelgridindex, int lower, double epsilon_trans);
double col_ionization(int modelgridindex, int phixstargetindex, double epsilon_trans);


static inline
double radfield(double nu, int modelgridindex)
/// calculates ambient radiation field, which is parameterised as a diluted black body
{
  float T_R = get_TR(modelgridindex);
  float W   = get_W(modelgridindex);

  return W * TWOHOVERCLIGHTSQUARED * pow(nu,3) * 1.0/(exp(HOVERKB*nu/T_R) - 1);
}


static inline
double radfield2(double nu, double T, double W)
/// calculates ambient radiation field, which is parameterised as a diluted black body
{
  return W * TWOHOVERCLIGHTSQUARED * pow(nu,3) * 1.0/(exp(HOVERKB*nu/T) - 1);
}


static inline
double get_individ_rad_deexc(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[i];
}


static inline
double get_individ_internal_down_same(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[i];
}


static inline
double get_individ_internal_up_same(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[i];
}


#endif //MACROATOM_H
