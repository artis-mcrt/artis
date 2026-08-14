// Declarations for LTE/approximate-NLTE level populations and ionisation balance (ltepop.cc).

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

// Level population [cm^-3]: NLTE where the solver provided one, else Boltzmann off the ground level.
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto calculate_levelpop(int nonemptymgi, int element, int ion, int level)
    -> double;

// Level population [cm^-3] forced to LTE, ignoring any NLTE solution. Unlike calculate_levelpop(), no MINPOP
// floor is applied to excited levels, so this can return an arbitrarily small population.
[[gnu::pure]] [[nodiscard]] auto calculate_levelpop_boltzmann(int nonemptymgi, int element, int ion, int level)
    -> double;

// Highest ion stage worth following in this cell, as an ion index; higher stages are treated as unpopulated.
// nne_hi is an assumed maximum electron density: the ion list is truncated where the running product of
// nne_hi * phi overflows. force_saha selects Saha even for elements with NLTE levels (which otherwise return
// all ions untruncated); pass it when the populations that follow will also be LTE.
[[nodiscard]] auto find_uppermost_ion(int nonemptymgi, int element, double nne_hi, bool force_saha) -> int;

// Solve for the free electron density self-consistently (the ion balance depends on nne, which is itself the
// sum of the ionic charges) and store nne and the ion populations in the grid.
void calculate_ion_balance_nne(int nonemptymgi);

void calculate_cellpartfuncts(int nonemptymgi, int element);

// Fractional ion abundances normalised to sum to one, indexed by ion index up to the cell's uppermost_ion.
// use_phi_saha selects the Saha balance rather than the photoionisation/recombination balance.
[[nodiscard]] auto calculate_ionfractions(int element, int nonemptymgi, double nne, bool use_phi_saha)
    -> std::vector<double>;

// Set ground level populations from Saha LTE or photoionisation equilibrium. Writes every ion of the element
// unconditionally, so callers must themselves skip elements whose populations the NLTE solver owns (or call it
// deliberately to overwrite them). force_saha as for find_uppermost_ion(); works up to the cell's stored
// uppermost_ion, so that must already be consistent with the balance chosen here.
void set_groundlevelpops(int nonemptymgi, int element, float nne, bool force_saha);

// Level population [cm^-3] read from the cell cache rather than recomputed: the fast path for packet
// propagation. Requires the cache to already hold this cell.
[[gnu::pure]] [[nodiscard]] inline DEVICE_FUNC auto get_cellcache_levelpop(const int nonemptymgi,
                                                                           const int uniquelevelindex) -> double {
  assert_testmodeonly(get_cellcache(nonemptymgi).nonemptymgi == nonemptymgi);

  const auto nn = get_cellcache(nonemptymgi).alllevels_pops[uniquelevelindex];
  assert_testmodeonly(nn >= 0.);
  return nn;
}

// as above, by (element, ion, level)
[[gnu::pure]] [[nodiscard]] inline DEVICE_FUNC auto get_cellcache_levelpop(const int nonemptymgi, const int element,
                                                                           const int ion, const int level) -> double {
  return get_cellcache_levelpop(nonemptymgi, get_uniquelevelindex(element, ion, level));
}

// Return the given ions groundlevel population for modelgridindex which was precalculated
// during update_grid and stored to the grid.
[[gnu::pure]] [[nodiscard]] inline DEVICE_FUNC auto get_groundlevelpop(const int nonemptymgi, const int element,
                                                                       const int ion) -> double {
  testmodeassert_valid_ion(element, ion);

  const double nn = grid::ion_groundlevelpops_allcells[get_cellionindex(nonemptymgi, element, ion)];
  if (nn < MINPOP) {
    if (grid::get_elem_massfrac(nonemptymgi, element) > 0) {
      return MINPOP;
    }
    return 0.;
  }
  return nn;
}

// Store an ion's ground level population for a cell (usually calculated during update_grid).
inline void set_groundlevelpop(const ptrdiff_t nonemptymgi, const int element, const int ion,
                               const float groundlevelpop) {
  grid::ion_groundlevelpops_allcells[get_cellionindex(nonemptymgi, element, ion)] = groundlevelpop;
}

// Return an ion's partition function for a cell, which was precalculated during update_grid.
[[gnu::pure]] [[nodiscard]] inline auto get_ion_partfunct(const ptrdiff_t nonemptymgi, const int element, const int ion)
    -> float {
  return grid::ion_partfuncts_allcells[get_cellionindex(nonemptymgi, element, ion)];
}

// Store an ion's partition function for a cell.
inline void set_ion_partfunct(const int nonemptymgi, const int element, const int ion, const float partfunct) {
  grid::ion_partfuncts_allcells[get_cellionindex(nonemptymgi, element, ion)] = partfunct;
}

// Use the ground level population and partition function to get an ion population
[[gnu::pure]] [[nodiscard]] inline auto get_nnion(const int nonemptymgi, const int element, const int ion) -> double {
  const auto nnion = get_groundlevelpop(nonemptymgi, element, ion) * get_ion_partfunct(nonemptymgi, element, ion) /
                     stat_weight(element, ion, 0);
  assert_testmodeonly(nnion >= 0.);
  assert_testmodeonly(std::isfinite(nnion));
  return nnion;
}

#endif  // LTEPOP_H
