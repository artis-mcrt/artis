#ifndef LTEPOP_H
#define LTEPOP_H

#include <cmath>
#include <cstddef>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "mpi_logging.h"

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto calculate_levelpop(int nonemptymgi, int element, int ion, int level)
    -> double;

[[gnu::pure]] [[nodiscard]] auto calculate_levelpop_boltzmann(int nonemptymgi, int element, int ion, int level)
    -> double;

[[nodiscard]] auto find_uppermost_ion(int nonemptymgi, int element, double nne_hi, bool force_saha) -> int;
void calculate_ion_balance_nne(int nonemptymgi);
void calculate_cellpartfuncts(int nonemptymgi, int element);
[[nodiscard]] auto calculate_ionfractions(int element, int nonemptymgi, double nne, bool use_phi_saha)
    -> std::vector<double>;
void set_groundlevelpops(int nonemptymgi, int element, float nne, bool force_saha);

[[gnu::pure]] [[nodiscard]] inline DEVICE_FUNC auto get_cellcache_levelpop(const int nonemptymgi,
                                                                           const int uniquelevelindex) -> double {
  assert_testmodeonly(get_cellcache(nonemptymgi).nonemptymgi == nonemptymgi);

  const auto nn = get_cellcache(nonemptymgi).alllevels_pops[uniquelevelindex];
  assert_testmodeonly(nn >= 0.);
  return nn;
}

// Calculate the population of a level from either LTE or NLTE information
[[gnu::pure]] [[nodiscard]] inline DEVICE_FUNC auto get_cellcache_levelpop(const int nonemptymgi, const int element,
                                                                           const int ion, const int level) -> double {
  return get_cellcache_levelpop(nonemptymgi, get_uniquelevelindex(element, ion, level));
}

// Return the given ions groundlevel population for modelgridindex which was precalculated
// during update_grid and stored to the grid.
[[gnu::pure]] [[nodiscard]] inline DEVICE_FUNC auto get_groundlevelpop(const int nonemptymgi, const int element,
                                                                       const int ion) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));

  const double nn = grid::ion_groundlevelpops_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_includedions()) +
                                                       get_uniqueionindex(element, ion)];
  if (nn < MINPOP) {
    if (grid::get_elem_massfrac(nonemptymgi, element) > 0) {
      return MINPOP;
    }
    return 0.;
  }
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

#endif  // LTEPOP_H
