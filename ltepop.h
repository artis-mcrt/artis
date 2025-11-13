#ifndef LTEPOP_H
#define LTEPOP_H

#include <cmath>
#include <cstddef>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "sn3d.h"

[[gnu::pure]] [[nodiscard]] auto get_groundlevelpop(int nonemptymgi, int element, int ion) -> double;

[[gnu::pure]] [[nodiscard]] auto calculate_levelpop(int nonemptymgi, int element, int ion, int level) -> double;

[[gnu::pure]] [[nodiscard]] auto calculate_levelpop_boltzmann(int nonemptymgi, int element, int ion, int level)
    -> double;

[[nodiscard]] auto find_uppermost_ion(int nonemptymgi, int element, double nne_hi, bool force_saha) -> int;
void calculate_ion_balance_nne(int nonemptymgi);
void calculate_cellpartfuncts(int nonemptymgi, int element);
[[nodiscard]] auto calculate_ionfractions(int element, int nonemptymgi, double nne, bool use_phi_saha)
    -> std::vector<double>;
void set_groundlevelpops(int nonemptymgi, int element, float nne, bool force_saha);

// calculates saha factor in LTE: Phi_level,ion,element = nn_level,ion,element/(nne*nn_upper,ion+1,element)
[[gnu::const]] [[nodiscard]] constexpr auto calculate_sahafact(const double g_lower, const double g_upper,
                                                               const double T, const double E_threshold) -> double {
  const double sahafact = SAHACONST * g_lower / g_upper * std::pow(T, -1.5) * std::exp(E_threshold / KB / T);
  assert_testmodeonly(std::isfinite(sahafact));
  return sahafact;
}
[[gnu::pure]] [[nodiscard]] inline auto get_cellcache_levelpop(const int nonemptymgi, const int uniquelevelindex)
    -> double {
  assert_testmodeonly(use_cellcache);
  assert_testmodeonly(globals::cellcache[cellcacheslotid].nonemptymgi == nonemptymgi);
  const auto nn = globals::cellcache[cellcacheslotid].alllevels_pops[uniquelevelindex];

  assert_testmodeonly(nn >= 0.);

  return nn;
}

// Calculate the population of a level from either LTE or NLTE information
[[gnu::pure]] [[nodiscard]] inline auto get_cellcache_levelpop(const int nonemptymgi, const int element, const int ion,
                                                               const int level) -> double {
  assert_testmodeonly(use_cellcache);
  assert_testmodeonly(globals::cellcache[cellcacheslotid].nonemptymgi == nonemptymgi);
  const auto nn = globals::cellcache[cellcacheslotid].alllevels_pops[get_uniquelevelindex(element, ion, level)];

  assert_testmodeonly(nn >= 0.);

  return nn;
}

// Use the ground level population and partition function to get an ion population
[[gnu::pure]] [[nodiscard]] inline auto get_nnion(const int nonemptymgi, const int element, const int ion) -> double {
  const auto nnion = get_groundlevelpop(nonemptymgi, element, ion) *
                     grid::ion_partfuncts_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) +
                                                   get_uniqueionindex(element, ion)] /
                     stat_weight(element, ion, 0);
  assert_testmodeonly(nnion >= 0.);
  assert_testmodeonly(std::isfinite(nnion));
  return nnion;
}

// Return the given ions groundlevel population for modelgridindex which was precalculated
// during update_grid and stored to the grid.
[[gnu::pure]] [[nodiscard]] inline auto get_groundlevelpop(const int nonemptymgi, const int element, const int ion)
    -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  const double nn = grid::ion_groundlevelpops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) +
                                                       get_uniqueionindex(element, ion)];
  if (nn < MINPOP) {
    if (grid::get_elem_abundance(nonemptymgi, element) > 0) {
      return MINPOP;
    }
    return 0.;
  }
  return nn;
}

#endif  // LTEPOP_H
