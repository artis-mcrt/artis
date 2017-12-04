#ifndef LTEPOP_H
#define LTEPOP_H

#include "atomic.h"
#include "sn3d.h"

double nne_solution_f(double x, void *restrict paras);
double ionfract(int element, int ion, int modelgridindex, double nne);
double phi(int element, int ion, int modelgridindex);
double calculate_partfunct(int element, int ion, int modelgridindex);
double get_groundlevelpop(int modelgridindex, int element, int ion);
double get_levelpop(int modelgridindex, int element, int ion, int level);
double calculate_levelpop_lte(int modelgridindex, int element, int ion, int level);
double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
double get_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
double get_groundmultiplet_pop(
  int modelgridindex, float T_e, int element, int ion, bool assume_lte);

void initialise_photoionestimators(void);

inline double calculate_sahafact(int element, int ion, int level, int upperionlevel, double T, double E_threshold)
/// calculates saha factor in LTE: Phi_level,ion,element = nn_level,ion,element/(nne*nn_upper,ion+1,element)
{
  const double sf = SAHACONST * stat_weight(element, ion, level) / stat_weight(element, ion + 1, upperionlevel) * pow(T, -1.5) * exp(E_threshold / (KB * T));
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


#endif //LTEPOP_H
