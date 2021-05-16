#ifndef LTEPOP_H
#define LTEPOP_H

#include "atomic.h"
#include "sn3d.h"

double nne_solution_f(double x, void *paras);
void get_ionfractions(int element, int modelgridindex, double nne, double ionfractions[], int uppermost_ion);
double phi(int element, int ion, int modelgridindex);
double calculate_partfunct(int element, int ion, int modelgridindex);
__host__ __device__ double get_groundlevelpop(int modelgridindex, int element, int ion);
__host__ __device__ double calculate_levelpop_lte(int modelgridindex, int element, int ion, int level);
__host__ __device__ double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
__host__ __device__ double get_groundmultiplet_pop(int modelgridindex, int element, int ion);


__host__ __device__
inline double calculate_sahafact(int element, int ion, int level, int upperionlevel, double T, double E_threshold)
/// calculates saha factor in LTE: Phi_level,ion,element = nn_level,ion,element/(nne*nn_upper,ion+1,element)
{
  const double g_lower = stat_weight(element, ion, level);
  const double g_upper = stat_weight(element, ion + 1, upperionlevel);
  const double sf = SAHACONST * g_lower / g_upper * pow(T, -1.5) * exp(E_threshold / KB / T);
  //printout("element %d, ion %d, level %d, T, %g, E %g has sf %g (g_l %g g_u %g)\n", element, ion, level, T, E_threshold, sf,stat_weight(element,ion,level),stat_weight(element,ion+1,0) );
  if (sf < 0)
  {
    printout("[fatal] calculate_sahafact: Negative Saha factor. sfac %g element %d ion %d level %d upperionlevel %d g_lower %g g_upper %g T %g E_threshold %g exppart %g\n",
             sf, element, ion, level, upperionlevel, g_lower, g_upper, T, E_threshold, exp(E_threshold / KB / T));
    abort();
  }
  return sf;
}


__host__ __device__
inline double ionstagepop(int modelgridindex, int element, int ion)
/// Calculates the given ionstages total population in nebular approximation for modelgridindex
/// The precalculated ground level population and partition function are used.
{
  return get_groundlevelpop(modelgridindex,element,ion) * grid::modelgrid[modelgridindex].composition[element].partfunct[ion]
          / stat_weight(element,ion,0);
}


#endif //LTEPOP_H
