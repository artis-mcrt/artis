#ifndef ATOMIC_H
#define ATOMIC_H

#include "sn3d.h"

int get_element(int element);
int get_elementindex(int Z);
int get_nions(int element);
int get_ionstage(int element, int ion);
int get_nlevels(int element, int ion);
int get_nlevels_nlte(int element, int ion);
int get_ionisinglevels(int element, int ion);
int get_bfcontinua(int element, int ion);
double epsilon(int element, int ion, int level);
double stat_weight(int element, int ion, int level);
short is_nlte(int element, int ion, int level);
int get_continuumindex(int element, int ion, int level);
int get_nphixstargets(int element, int ion, int level);
int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
float get_phixsprobability(int element, int ion, int level, int phixstargetindex);

//inline functions:

inline
int transitioncheck(int upper, int lower)
/// reads A_ul from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
  int index = (upper - lower) - 1;
  int flag = transitions[upper].to[index];

  return flag;
}


inline
double einstein_spontaneous_emission(int lineindex)
//double einstein_spontaneous_emission(int element, int ion, int upper, int lower)
/// reads A_ul from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
/*
  int index = (upper-lower) - 1;
  double A_ul = elements[element].ions[ion].levels[upper].transitions[index].einstein_A;
*/
  double A_ul = linelist[lineindex].einstein_A;

  return A_ul;
}


inline
double osc_strength(int lineindex)
//double osc_strength(int element, int ion, int upper, int lower)
/// reads f_lu from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
  return linelist[lineindex].osc_strength;
}

inline
double coll_str(int lineindex)
{
  return linelist[lineindex].coll_str;
}

inline
double statw_up(int lineindex)
{
  return elements[linelist[lineindex].elementindex].ions[linelist[lineindex].ionindex].levels[linelist[lineindex].upperlevelindex].stat_weight;
}

inline
double statw_down(int lineindex)
{
  return elements[linelist[lineindex].elementindex].ions[linelist[lineindex].ionindex].levels[linelist[lineindex].lowerlevelindex].stat_weight;
}

/*inline
double osc_strength_old(int lineindex)
//double osc_strength(int element, int ion, int upper, int lower)
/// reads f_lu from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{

  int index = (upper-lower) - 1;
  double f_ul = elements[element].ions[ion].levels[upper].transitions[index].oscillator_strength;

  return f_ul;
}*/

double photoionization_crosssection(double nu_edge, double nu);

#endif //ATOMIC_H
