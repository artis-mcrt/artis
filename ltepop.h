#ifndef LTEPOP_H
#define LTEPOP_H

#include "atomic.h"

double nne_solution_f(double x, void *restrict paras);
double ionfract(int element, int ion, int modelgridindex, double nne);
double phi(int element, int ion, int modelgridindex);
double calculate_partfunct(int element, int ion, int modelgridindex);
double get_groundlevelpop(int modelgridindex, int element, int ion);
double calculate_levelpop_lte(int modelgridindex, int element, int ion, int level);
double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
double get_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
void initialise_photoionestimators(void);

inline double calculate_sahafact(int element, int ion, int level, int upperionlevel, double T, double E_threshold)
/// calculates saha factor in LTE: Phi_level,ion,element = nn_level,ion,element/(nne*nn_upper,ion+1,element)
{
  double sf = SAHACONST * stat_weight(element,ion,level) / (stat_weight(element,ion+1,upperionlevel) * T * sqrt(T)) * exp(E_threshold/(KB*T));
  //printout("element %d, ion %d, level %d, T, %g, E %g has sf %g (g_l %g g_u %g)\n", element, ion, level, T, E_threshold, sf,stat_weight(element,ion,level),stat_weight(element,ion+1,0) );
  if (sf < 0)
  {
    printout("[fatal] sahafact: negative saha factor");
    abort();
  }
  return sf;
}


inline double ionstagepop(int modelgridindex, int element, int ion)
/// Calculates the given ionstages total population in nebular approximation for modelgridindex
/// The precalculated ground level population and partition function are used.
{
  return get_groundlevelpop(modelgridindex,element,ion) * modelgrid[modelgridindex].composition[element].partfunct[ion]
          / stat_weight(element,ion,0);
}


inline double get_levelpop(int modelgridindex, int element, int ion, int level)
/// Returns the given levels occupation number, which are stored in the active
/// entry of the cellhistory.
{
//printout("get_levelpop histindex %d\n",histindex);
  if (use_cellhist)
  {
    double pop = cellhistory[tid].chelements[element].chions[ion].chlevels[level].population;
    if (pop > -1)
      return pop;
    else
    {
      printout("Abort: get_levelpop called, but no population in cellhistory for element %d ion %d level %d",element,ion,level);
      abort();
    }
    int cellmgi = cell[cellhistory[tid].cellnumber].modelgridindex;

    if (cellmgi != modelgridindex)
    {
      printout("Abort: get_levelpop called, but cellhistory mgi %d != argument modelgridindex %d",
               cellmgi,modelgridindex);
      abort();
    }
  }
  else
  {
    return calculate_exclevelpop(modelgridindex,element,ion,level);
  }
}

#endif //LTEPOP_H
