#ifndef ATOMIC_H
#define ATOMIC_H

#include "sn3d.h"

double last_phixs_nuovernuedge; // last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge

double get_tau_sobolev(int modelgridindex, int lineindex, double t_current);
int get_tot_nions(void);
double photoionization_crosssection_fromtable(float *photoion_xs, double nu_edge, double nu);

inline int get_element(int element)
/// Returns the atomic number associated with a given elementindex.
{
  return elements[element].anumber;
}


inline int get_elementindex(int Z)
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
  //abort();;
  return -100;
}


inline int get_nions(int element)
/// Returns the number of ions associated with a specific element given by
/// its elementindex.
{
  return elements[element].nions;
}


inline int get_ionstage(int element, int ion)
/// Returns the ionisation stage of an ion specified by its elementindex and
/// ionindex.
{
  return elements[element].ions[ion].ionstage;
}


inline int get_nlevels(int element, int ion)
/// Returns the number of levels associated with with a specific ion given
/// its elementindex and ionindex.
{
  return elements[element].ions[ion].nlevels;
}


inline int get_nlevels_nlte(int element, int ion)
// Returns the number of levels associated with with a specific ion given
// its elementindex and ionindex.
{
  return elements[element].ions[ion].nlevels_nlte;
}


inline int get_nlevels_groundterm(int element, int ion)
{
  return elements[element].ions[ion].nlevels_groundterm;
}


inline int get_ionisinglevels(int element, int ion)
/// Returns the number of levels associated with an ion that
/// have energies below the ionisation threshold.
{
  return elements[element].ions[ion].ionisinglevels;
}


inline int get_uniqueionindex(const int element, const int ion)
// Get an index for an ion that is unique for every ion of every element
{
  int index = 0;
  for (int e = 0; e < element; e++)
  {
    index += get_nions(element);
  }
  index += ion;

  // assert(index < get_tot_nions());
  return index;
}


inline void get_ionfromuniqueionindex(const int allionsindex, int *element, int *ion)
{
  int allionsindex_thiselementfirstion = 0;
  for (int e = 0; e < nelements; e++)
  {
    if ((allionsindex - allionsindex_thiselementfirstion) >= get_nions(e))
    {
      allionsindex_thiselementfirstion += get_nions(e);
    }
    else
    {
      *element = e;
      *ion = allionsindex - allionsindex_thiselementfirstion;
      assert(get_uniqueionindex(*element, *ion) == allionsindex);
      return;
    }
  }
  *element = -1;
  *ion = -1;
}


inline double epsilon(int element, int ion, int level)
/// Returns the energy of (element,ion,level).
{
  return elements[element].ions[ion].levels[level].epsilon;
}


inline double stat_weight(int element, int ion, int level)
/// Returns the statistical weight of (element,ion,level).
{
  #ifdef DEBUG_ON
  if (level > elements[element].ions[ion].nlevels)
  {
    printout("[fatal] stat_weight: level %d greater than nlevels=%d ... abort\n",level,elements[element].ions[ion].nlevels);
    abort();
  }
  #endif
  return elements[element].ions[ion].levels[level].stat_weight;
}


inline int get_maxrecombininglevel(int element, int ion)
/// Returns the number of bf-continua associated with ion ion of element element.
{
  return elements[element].ions[ion].maxrecombininglevel;
}


inline bool is_nlte(int element, int ion, int level)
// Returns true if (element,ion,level) is to be treated in nlte.
// (note this function returns true for the ground state,
//  although it is stored separately from the excited NLTE states)
{
//  if (get_element(element) == 26 && get_ionstage(element, ion) == 2)
//    return (level <= 197);
//  else
//    return (level <= 80);
   return (level <= 300);
}


inline int get_continuumindex(int element, int ion, int level)
/// Returns the index of the continuum associated to the given level.
{
  return elements[element].ions[ion].levels[level].cont_index;
}


inline int get_ndowntrans(int element, int ion, int level)
// the number of downward bound-bound transitions from the specified level
{
  return elements[element].ions[ion].levels[level].ndowntrans;
}


inline int get_nuptrans(int element, int ion, int level)
// the number of upward bound-bound transitions from the specified level
{
  return elements[element].ions[ion].levels[level].nuptrans;
}


inline void set_ndowntrans(const int element, const int ion, const int level, const int ndowntrans)
// the number of downward bound-bound transitions from the specified level
{
  elements[element].ions[ion].levels[level].ndowntrans = ndowntrans;
}


inline void set_nuptrans(const int element, const int ion, const int level, const int nuptrans)
// the number of upward bound-bound transitions from the specified level
{
  elements[element].ions[ion].levels[level].nuptrans = nuptrans;
}


inline int get_nphixstargets(const int element, const int ion, const int level)
/// Returns the number of target states for photoionization of (element,ion,level).
{
  const int nions = get_nions(element);
  const int nionisinglevels = get_ionisinglevels(element,ion);
  if ((ion < nions-1) && (level < nionisinglevels))
    return elements[element].ions[ion].levels[level].nphixstargets;
  else
  {
    return 0;
  }
}


inline int get_phixsupperlevel(const int element, const int ion, const int level, const int phixstargetindex)
/// Returns the level index of a target state for photoionization of (element,ion,level).
{
  #ifdef DEBUG_ON
    if ((phixstargetindex < 0) || (phixstargetindex > get_nphixstargets(element,ion,level)-1))
    {
      printout("[fatal]   get_phixsupperlevel called with invalid phixstargetindex\n");
      printout("arguments: element %d, ion %d, level %d phixstargetindex %d, nphixstargets %d\n",element,ion,level,phixstargetindex,get_nphixstargets(element,ion,level));
      abort();
    }
  #endif

  return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].levelindex;
}


inline double get_phixs_threshold(int element, int ion, int level, int phixstargetindex)
/// Returns the energy of (element,ion,level).
{
  // const double phixs_threshold_stored = elements[element].ions[ion].levels[level].phixs_threshold;
  // if (phixs_threshold_stored > 0.)
  //   return phixs_threshold_stored;
  // else
  {
    const int upperlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
    const double E_threshold = epsilon(element, ion + 1, upperlevel) - epsilon(element, ion, level);
    return E_threshold;
  }
}


inline double get_phixsprobability(int element, int ion, int level, int phixstargetindex)
/// Returns the probability of a target state for photoionization of (element,ion,level).
{
  #ifdef DEBUG_ON
    if ((phixstargetindex < 0) || (phixstargetindex >= get_nphixstargets(element,ion,level)))
    {
      printout("[fatal]   get_phixsprobability called with invalid phixstargetindex");
      //printout("arguments: element %d, ion %d, level %d phixstargetindex %g, nphixstargets %g\n",element,ion,level,phixstargetindex,get_nphixstargets(element,ion,level));
      abort();
    }
  #endif

  return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].probability;
}


inline double einstein_spontaneous_emission(int lineindex)
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


inline double osc_strength(int lineindex)
//double osc_strength(int element, int ion, int upper, int lower)
/// reads f_lu from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
  return linelist[lineindex].osc_strength;
}


inline double get_coll_str(int lineindex)
{
  return linelist[lineindex].coll_str;
}


inline double statw_upper(int lineindex)
{
  const int element = linelist[lineindex].elementindex;
  const int ion = linelist[lineindex].ionindex;
  const int upper = linelist[lineindex].upperlevelindex;
  return elements[element].ions[ion].levels[upper].stat_weight;
}


inline double statw_lower(int lineindex)
{
  const int element = linelist[lineindex].elementindex;
  const int ion = linelist[lineindex].ionindex;
  const int lower = linelist[lineindex].lowerlevelindex;
  return elements[element].ions[ion].levels[lower].stat_weight;
}


inline double photoionization_crosssection(int element, int ion, int level, double nu_edge, double nu)
{
  return photoionization_crosssection_fromtable(elements[element].ions[ion].levels[level].photoion_xs, nu_edge, nu);
}


inline double photoionization_crosssection_macroatom(double nu_edge, double nu)
///        - BE AWARE: the elements of the global structure variable mastate
///                    must fit to the bound state of the desired bf-continuum!!!
{
  const int element = mastate[tid].element;
  const int ion = mastate[tid].ion;
  const int level = mastate[tid].level;

  return photoionization_crosssection_fromtable(elements[element].ions[ion].levels[level].photoion_xs, nu_edge, nu);
}

/*static double osc_strength_old(int lineindex)
//double osc_strength(int element, int ion, int upper, int lower)
/// reads f_lu from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{

  int index = (upper-lower) - 1;
  double f_ul = elements[element].ions[ion].levels[upper].transitions[index].oscillator_strength;

  return f_ul;
}*/


#endif //ATOMIC_H
