#ifndef LTEPOP_H
#define LTEPOP_H

#include "atomic.h"
#include "sn3d.h"

double nne_solution_f(double x, void *restrict paras);
void get_ionfractions(int element, int modelgridindex, double nne, double ionfractions[], int uppermost_ion);
double phi(int element, int ion, int modelgridindex);
double calculate_partfunct(int element, int ion, int modelgridindex);
double get_groundlevelpop(int modelgridindex, int element, int ion);
double get_levelpop(int modelgridindex, int element, int ion, int level);
double calculate_levelpop_lte(int modelgridindex, int element, int ion, int level);
double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
double get_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
double get_groundmultiplet_pop(
  int modelgridindex, float T_e, int element, int ion, bool assume_lte);

typedef enum {
// ALEXEI:
// enum to indicate fields within cellhistory_struct which could be reset
  sahafact_mask = 1,
  //spontaneousrecombrate_mask = 2,
  //bfcooling_mask = 4,
  //bfheatingcoeff_mask = 8,
  //corrphotoioncoeff_mask = 16,
  //population_mask = 32,
  //processrates_mask = 64,
  //cooling_contribution_mask = 128,
  count,
  cell_reset_mask = 1
} cellhist_reset_field;

__attribute__((always_inline)) inline bool check_cellhist_param_reset(int reset_mask, cellhist_reset_field field_mask) {
/// ALEXEI:
/// Checks whether relevant field within cellhist has been reset (e.g. for a new timestep) 
/// by comparing with bitmask.
  
  bool reset = false;
  
  // compare to reset mask
  reset = (reset_mask & field_mask);

  return reset;
}

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


inline double ionstagepop(int modelgridindex, int element, int ion)
/// Calculates the given ionstages total population in nebular approximation for modelgridindex
/// The precalculated ground level population and partition function are used.
{
  return get_groundlevelpop(modelgridindex,element,ion) * modelgrid[modelgridindex].composition[element].partfunct[ion]
          / stat_weight(element,ion,0);
}


#endif //LTEPOP_H
