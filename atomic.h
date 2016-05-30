#ifndef ATOMIC_H
#define ATOMIC_H

#include "sn3d.h"

double last_phixs_nuovernuedge; // last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge


static inline
int get_element(int element)
/// Returns the atomic number associated with a given elementindex.
{
  return elements[element].anumber;
}


static inline
int get_elementindex(int Z)
/// Returns the elementindex associated with a given atomic number.
/// If there is no element with the given atomic number in the atomic data
/// a negative value is returned to flag this event.
{
  for (int i = 0; i < nelements; i++)
  {
    //printf("i %d, Z %d, elements[i].anumber %d\n",i,Z,elements[i].anumber);
    if (Z == elements[i].anumber)
      return i;
  }

  //printout("[debug] get_elementindex: element Z=%d was not found in atomic data ... skip readin of cross sections for this element\n",Z);
  //printout("[fatal] get_elementindex: element Z=%d was not found in atomic data ... abort\n");
  //exit(0);
  return -100;
}


static inline
int get_nions(int element)
/// Returns the number of ions associated with a specific element given by
/// its elementindex.
{
  return elements[element].nions;
}


static inline
int get_ionstage(int element, int ion)
/// Returns the ionisation stage of an ion specified by its elementindex and
/// ionindex.
{
  return elements[element].ions[ion].ionstage;
}


static inline
int get_nlevels(int element, int ion)
/// Returns the number of levels associated with with a specific ion given
/// its elementindex and ionindex.
{
  return elements[element].ions[ion].nlevels;
}


static inline
int get_nlevels_nlte(int element, int ion)
/// Returns the number of levels associated with with a specific ion given
/// its elementindex and ionindex.
{
  return elements[element].ions[ion].nlevels_nlte;
}


static inline
int get_ionisinglevels(int element, int ion)
/// Returns the number of levels associated with an ion that
/// have energies below the ionisation threshold.
{
  return elements[element].ions[ion].ionisinglevels;
}


static inline
double epsilon(int element, int ion, int level)
/// Returns the energy of (element,ion,level).
{
  return elements[element].ions[ion].levels[level].epsilon;
}


static inline
double stat_weight(int element, int ion, int level)
/// Returns the statistical weight of (element,ion,level).
{
  #ifdef DEBUG_ON
  if (level > elements[element].ions[ion].nlevels)
  {
    printout("[fatal] stat_weight: level %d greater than nlevels=%d ... abort\n",level,elements[element].ions[ion].nlevels);
    exit(0);
  }
  #endif
  return elements[element].ions[ion].levels[level].stat_weight;
}


static inline
int get_bfcontinua(int element, int ion)
/// Returns the number of bf-continua associated with ion ion of element element.
{
  int nionisinglevels = get_ionisinglevels(element,ion);

  if (nionisinglevels < max_bf_continua)
    return nionisinglevels;
  else
    return max_bf_continua;
}


static inline
bool is_nlte(int element, int ion, int level)
// Returns true if (element,ion,level) is to be treated in nlte.
// (note this function returns true for the ground state)
{
  if (level < 100) //TODO: change back to 200
    elements[element].ions[ion].levels[level].is_nlte = true;
  else
    elements[element].ions[ion].levels[level].is_nlte = false;

  return elements[element].ions[ion].levels[level].is_nlte;
}


static inline
int get_continuumindex(int element, int ion, int level)
/// Returns the index of the continuum associated to the given level.
{
  return elements[element].ions[ion].levels[level].cont_index;
}


static inline
int get_nphixstargets(int element, int ion, int level)
/// Returns the number of target states for photoionization of (element,ion,level).
{
  int nions = get_nions(element);
  int nionisinglevels = get_ionisinglevels(element,ion);
  if ((ion < nions-1) && (level < nionisinglevels))
    return elements[element].ions[ion].levels[level].nphixstargets;
  else
  {
    return 0;
  }
}


static inline
int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex)
/// Returns the level index of a target state for photoionization of (element,ion,level).
{
  #ifdef DEBUG_ON
    if ((phixstargetindex < 0) || (phixstargetindex > get_nphixstargets(element,ion,level)-1))
    {
      printout("[fatal]   get_phixsupperlevel called with invalid phixstargetindex");
      printout("arguments: element %d, ion %d, level %d phixstargetindex %d, nphixstargets %d\n",element,ion,level,phixstargetindex,get_nphixstargets(element,ion,level));
      abort();
    }
  #endif

  return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].levelindex;
}


static inline
double get_phixsprobability(int element, int ion, int level, int phixstargetindex)
/// Returns the probability of a target state for photoionization of (element,ion,level).
{
  #ifdef DEBUG_ON
    if ((phixstargetindex < 0) || (phixstargetindex >= get_nphixstargets(element,ion,level)))
    {
      printout("[fatal]   get_phixsprobability called with invalid phixstargetindex");
      printout("arguments: element %d, ion %d, level %d phixstargetindex %g, nphixstargets %g\n",element,ion,level,phixstargetindex,get_nphixstargets(element,ion,level));
      abort();
    }
  #endif

  return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].probability;
}


static inline
int transitioncheck(int upper, int lower)
/// reads A_ul from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
  int index = (upper - lower) - 1;
  int flag = transitions[upper].to[index];

  return flag;
}


static inline
double einstein_spontaneous_emission(int lineindex)
//double einstein_spontaneous_emission(int element, int ion, int upper, int lower)
/// reads A_ul from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
/*
  int index = (upper-lower) - 1;
  double A_ul = elements[element].ions[ion].levels[upper].transitions[index].einstein_A;
*/
  return linelist[lineindex].einstein_A;
}


static inline
double osc_strength(int lineindex)
//double osc_strength(int element, int ion, int upper, int lower)
/// reads f_lu from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
  return linelist[lineindex].osc_strength;
}

static inline
double coll_str(int lineindex)
{
  return linelist[lineindex].coll_str;
}

static inline
double statw_up(int lineindex)
{
  return elements[linelist[lineindex].elementindex].ions[linelist[lineindex].ionindex].levels[linelist[lineindex].upperlevelindex].stat_weight;
}

static inline
double statw_down(int lineindex)
{
  return elements[linelist[lineindex].elementindex].ions[linelist[lineindex].ionindex].levels[linelist[lineindex].lowerlevelindex].stat_weight;
}

/*static inline
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
