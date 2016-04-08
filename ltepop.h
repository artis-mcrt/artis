#ifndef LTEPOP_H
#define LTEPOP_H

#include "atomic.h"

double nne_solution_f(double x, void *paras);
double ionfract(int element, int ion, int modelgridindex, double nne);
double phi(int element, int ion, int modelgridindex);
double calculate_partfunct_old(int element, int ion, int modelgridindex);
double calculate_partfunct(int element, int ion, int modelgridindex);
double get_groundlevelpop(int modelgridindex, int element, int ion);
double calculate_exclevelpop_old(int modelgridindex, int element, int ion, int level);
double calculate_levelpop_lte(int modelgridindex, int element, int ion, int level);
double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);

static inline
double ionstagepop(int modelgridindex, int element, int ion)
/// Calculates the given ionstages total population in nebular approximation for modelgridindex
/// The precalculated ground level population and partition function are used.
{
  return get_groundlevelpop(modelgridindex,element,ion) * modelgrid[modelgridindex].composition[element].partfunct[ion]
          / stat_weight(element,ion,0);
}

//void calculate_levelpops(int modelgridindex);

static inline
double get_levelpop(int modelgridindex, int element, int ion, int level)
/// Returns the given levels occupation number, which are stored in the active
/// entry of the cellhistory.
{
//printout("get_levelpop histindex %d\n",histindex);
  if (use_cellhist >= 0)
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

double calculate_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
double get_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
void initialise_photoionestimators(void);

#endif //LTEPOP_H
