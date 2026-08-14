// Inline accessors for the atomic dataset: elements, ions, levels, transitions, and
// photoionisation cross-section tables.

#ifndef ATOMIC_H
#define ATOMIC_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <tuple>

#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "mpi_logging.h"

// highest number of ions for any element
inline int maxnions = 0;

// number of ions of all elements combined
inline int includedions = 0;

// number of levels of all elements combined
inline int includedlevels = 0;

// last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge
inline double last_phixs_nuovernuedge;

// first value in this array is not used but exists so the indexes match those of the phixsdata_filenames array
inline std::array<bool, 3> phixs_file_version_exists;

inline const std::array phixsdata_filenames{"IGNORE", "phixsdata.txt", "phixsdata_v2.txt"};

// return the number of levels of all elements combined
[[gnu::pure]] inline auto get_includedlevels() -> int { return includedlevels; }

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto get_nelements() -> int {
  return static_cast<int>(globals::elements.size());
}

// Bounds check (active only in TESTMODE builds) for an elementindex
DEVICE_FUNC inline void testmodeassert_valid_element([[maybe_unused]] const int element) {
  assert_testmodeonly(element >= 0);
  assert_testmodeonly(element < get_nelements());
}

// total density of nuclei
[[gnu::pure]] DEVICE_FUNC inline auto get_nnion_tot(int nonemptymgi) -> double {
  double nntot = 0.;
  for (int element = 0; element < get_nelements(); element++) {
    nntot += grid::get_elem_numberdens(nonemptymgi, element);
  }

  return nntot;
}

// Return the number of ions associated with a specific element given by its elementindex.
[[gnu::pure]] inline auto get_nions(const int element) -> int {
  testmodeassert_valid_element(element);
  return static_cast<int>(globals::elements[element].ions.size());
}

// Bounds check (active only in TESTMODE builds) for an elementindex and ionindex
DEVICE_FUNC inline void testmodeassert_valid_ion([[maybe_unused]] const int element, [[maybe_unused]] const int ion) {
  testmodeassert_valid_element(element);
  assert_testmodeonly(ion >= 0);
  assert_testmodeonly(ion < get_nions(element));
}

// Return the number of levels associated with a specific ion given its elementindex and ionindex.
[[gnu::pure]] DEVICE_FUNC inline auto get_nlevels(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].ions[ion].nlevels;
}

// Bounds check (active only in TESTMODE builds) for elementindex, ionindex, and levelindex
DEVICE_FUNC inline void testmodeassert_valid_level([[maybe_unused]] const int element, [[maybe_unused]] const int ion,
                                                   [[maybe_unused]] const int level) {
  testmodeassert_valid_ion(element, ion);
  assert_testmodeonly(level >= 0);
  assert_testmodeonly(level < get_nlevels(element, ion));
}

// get the uniquelevelindex of the lowest level of an ion (excited levels are offset from this)
[[nodiscard]] DEVICE_FUNC inline auto get_ionuniquelevelindexstart(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].ions[ion].uniquelevelindexstart;
}

// get the index into the nltepops_allcells array of the lowest NLTE level of an ion
// higher levels are consecutive offsets from this
[[nodiscard]] DEVICE_FUNC inline auto get_allnltelevelsindexstart(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].ions[ion].allnltelevelsindexstart;
}

// Get an index for level of an ionstage of an element that is unique across every ion of every element
[[nodiscard]] inline auto get_uniquelevelindex(const int element, const int ion, const int level) -> int {
  testmodeassert_valid_level(element, ion, level);
  const auto uniquelevelindex = get_ionuniquelevelindexstart(element, ion) + level;
  assert_testmodeonly(uniquelevelindex < get_includedlevels());

  return uniquelevelindex;
}

// Return the ionisation stage of an ion specified by its elementindex and ionindex.
[[nodiscard]] DEVICE_FUNC inline auto get_ionstage(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].lowest_ionstage + ion;
}

// Return the ionisation potential of an ion in erg.
[[nodiscard]] DEVICE_FUNC inline auto get_ionpot(const int element, const int ion) -> double {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].ions[ion].ionpot;
}

// Return the index in the groundcont array for the ground state photoion continuum of an ion.
[[nodiscard]] DEVICE_FUNC inline auto get_groundcontindex(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].ions[ion].groundcontindex;
}

// Return the number of levels associated with an ion that have energies below the ionisation threshold.
[[nodiscard]] DEVICE_FUNC inline auto get_nlevels_ionising(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].ions[ion].nlevels_ionising;
}

// Get the number of target states for photoionisation of a level, by unique level index or by
// (element, ion, level).
DEVICE_FUNC inline auto get_nphixstargets(const int uniquelevelindex) -> int {
  return globals::alllevels.nphixstargets[uniquelevelindex];
}

DEVICE_FUNC inline auto get_nphixstargets(const int element, const int ion, const int level) -> int {
  testmodeassert_valid_ion(element, ion);
  const auto nphixstargets = get_nphixstargets(get_uniquelevelindex(element, ion, level));
  assert_testmodeonly(nphixstargets == 0 ||
                      ((ion < (get_nions(element) - 1)) && (level < get_nlevels_ionising(element, ion))));
  return nphixstargets;
}

// Return the index into the allphixstargets arrays for a target state for photoionisation of (element,ion,level).
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto get_allphixstargetindex(const int uniquelevelindex,
                                                                            const int phixstargetindex) -> int {
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(uniquelevelindex));

  return globals::alllevels.phixstargetstart[uniquelevelindex] + phixstargetindex;
}
// Return the upper-ion level index of a photoionisation target state, by unique level index or by
// (element, ion, level).
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto get_phixsupperlevel(const int uniquelevelindex,
                                                                        const int phixstargetindex) -> int {
  return globals::allphixstargets_levelindex[get_allphixstargetindex(uniquelevelindex, phixstargetindex)];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto get_phixsupperlevel(const int element, const int ion,
                                                                        const int level, const int phixstargetindex)
    -> int {
  return get_phixsupperlevel(get_uniquelevelindex(element, ion, level), phixstargetindex);
}

// Return the branching probability of a photoionisation target state (the targets of a level sum to one),
// by unique level index or by (element, ion, level).
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto get_phixsprobability(const int uniquelevelindex,
                                                                         const int phixstargetindex) -> double {
  return globals::allphixstargets_probability[get_allphixstargetindex(uniquelevelindex, phixstargetindex)];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto get_phixsprobability(const int element, const int ion,
                                                                         const int level, const int phixstargetindex)
    -> double {
  testmodeassert_valid_level(element, ion, level);

  return get_phixsprobability(get_uniquelevelindex(element, ion, level), phixstargetindex);
}

// Return the level index of the highest level with a non-zero recombination rate for ion ion of element element.
[[gnu::pure]] [[nodiscard]] inline auto get_maxrecombininglevel(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].ions[ion].maxrecombininglevel;
}

[[gnu::pure]] [[nodiscard]] inline DEVICE_FUNC auto get_phixs_table(const int uniquelevelindex)
    -> std::span<const float> {
  const ptrdiff_t phixsstart = globals::alllevels.phixsstart[uniquelevelindex];
  assert_testmodeonly(phixsstart >= 0);
  return globals::allphixs.subspan(phixsstart * globals::NPHIXSPOINTS, globals::NPHIXSPOINTS);
}

[[gnu::pure]] [[nodiscard]] inline DEVICE_FUNC auto get_phixs_table(const int element, const int ion, const int level)
    -> std::span<const float> {
  return get_phixs_table(get_uniquelevelindex(element, ion, level));
}

// Calculate the photoionisation cross-section at frequency nu out of the atomic data.
[[gnu::pure]] [[nodiscard]] inline auto photoionisation_crosssection_fromtable(std::span<const float> photoion_xs,
                                                                               const double nu_edge, const double nu)
    -> float {
  float sigma_bf = 0.;

  if constexpr (PHIXS_CLASSIC_NO_INTERPOLATION) {
    // classic mode: no interpolation
    if (nu < nu_edge) {
      // below the threshold the cross section is zero (and the table index below would be negative)
      sigma_bf = 0.;
    } else if (nu == nu_edge) {
      sigma_bf = photoion_xs[0];
    } else if (nu < nu_edge * (1 + (globals::NPHIXSNUINCREMENT * globals::NPHIXSPOINTS))) {
      // the range guard above and the index below are computed with different floating-point
      // expressions, so for nu just under the bound the division can round up to exactly
      // NPHIXSPOINTS. Clamp so that the read stays inside this level's table
      const int i = std::min(static_cast<int>((nu - nu_edge) / (globals::NPHIXSNUINCREMENT * nu_edge)),
                             globals::NPHIXSPOINTS - 1);
      sigma_bf = photoion_xs[i];
    } else {
      // above the top of the table, extrapolate with the Kramers (1923) nu^-3 scaling. It is anchored to the
      // highest tabulated point rather than to the threshold value so that the cross-section stays continuous
      // across the end of the table.
      sigma_bf = static_cast<float>(photoion_xs[globals::NPHIXSPOINTS - 1] *
                                    pow(nu_edge * (1 + (globals::NPHIXSNUINCREMENT * globals::NPHIXSPOINTS)) / nu, 3));
    }
    return sigma_bf;
  }

  const double ireal = ((nu / nu_edge) - 1.0) / globals::NPHIXSNUINCREMENT;
  // floor() so that nu < nu_edge always gives i < 0 (zero cross-section below the threshold);
  // truncation toward zero would map ireal in (-1, 0) to the first table point instead
  const int i = static_cast<int>(std::floor(ireal));

  if (i < 0) {
    sigma_bf = 0.;
  } else if (i < globals::NPHIXSPOINTS - 1) {
    const double sigma_bf_a = photoion_xs[i];
    const double sigma_bf_b = photoion_xs[i + 1];
    const double factor_b = ireal - i;
    sigma_bf = static_cast<float>(((1. - factor_b) * sigma_bf_a) + (factor_b * sigma_bf_b));
  } else {
    // above the top of the table, extrapolate with the Kramers (1923) nu^-3 scaling. It is anchored to the
    // highest tabulated point rather than to the threshold value so that the cross-section stays continuous
    // across the end of the table.
    const double nu_max_phixs = nu_edge * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table
    sigma_bf = static_cast<float>(photoion_xs[globals::NPHIXSPOINTS - 1] * pow3(nu_max_phixs / nu));
  }

  return sigma_bf;
}

// inverse of get_uniquelevelindex(). get the element/ion/level from a unique level index
[[gnu::pure]] [[nodiscard]] inline auto get_levelfromuniquelevelindex(const int uniquelevelindex)
    -> std::tuple<int, int, int> {
  assert_testmodeonly(uniquelevelindex < get_includedlevels());
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      if (get_nlevels(element, ion) == 0) {
        continue;
      }
      const int level = uniquelevelindex - get_ionuniquelevelindexstart(element, ion);
      if (level >= 0 && level < get_nlevels(element, ion)) {
        assert_testmodeonly(get_uniquelevelindex(element, ion, level) == uniquelevelindex);
        return {element, ion, level};
      }
    }
  }
  assert_always(false);  // uniquelevelindex too high to be valid
  return {-1, -1, -1};
}
// Return the statistical weight of a level, by unique level index or by (element, ion, level).
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto stat_weight(const int uniquelevelindex) -> double {
  return globals::alllevels.statweight[uniquelevelindex];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto stat_weight(const int element, const int ion, const int level)
    -> double {
  testmodeassert_valid_level(element, ion, level);
  return stat_weight(get_uniquelevelindex(element, ion, level));
}

// Return the energy of (element,ion,level).
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto epsilon(const int uniquelevelindex) -> double {
  return globals::alllevels.epsilon[uniquelevelindex];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto epsilon(const int element, const int ion, const int level)
    -> double {
  return epsilon(get_uniquelevelindex(element, ion, level));
}

// Get the atomic number associated with a given elementindex.
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC inline auto get_atomicnumber(const int element) -> int {
  testmodeassert_valid_element(element);
  return globals::elements[element].anumber;
}

// Return true if (element,ion,level) is to be treated in NLTE.
// (note this function returns true for the ground state,
//  although it is stored separately from the excited NLTE states)
[[gnu::pure]] [[nodiscard]] inline auto is_nlte(const int element, const int ion, const int level) -> bool {
  testmodeassert_valid_level(element, ion, level);
  return level <= ION_NLEVELS_EXCITED_NLTE(get_atomicnumber(element), get_ionstage(element, ion));
}

// Return the elementindex associated with a given atomic number, or -1 if there is no element with the
// given atomic number in the atomic data.
[[gnu::pure]] [[nodiscard]] inline auto get_elementindex(const int Z) -> int {
  const auto elem =
      std::ranges::find_if(globals::elements, [Z](const Element& element) { return element.anumber == Z; });
  if (elem != globals::elements.end()) {
    return static_cast<int>(elem - globals::elements.begin());
  }

  return -1;
}

inline void update_includedionslevels_maxnions() {
  includedions = 0;
  includedlevels = 0;
  maxnions = 0;
  for (int element = 0; element < get_nelements(); element++) {
    includedions += get_nions(element);
    maxnions = std::max(maxnions, get_nions(element));
    for (int ion = 0; ion < get_nions(element); ion++) {
      includedlevels += get_nlevels(element, ion);
    }
  }
}

// return the number of ions of all elements combined
[[gnu::pure]] [[nodiscard]] inline auto get_includedions() -> int {
  assert_testmodeonly(includedions > 0);
  return includedions;
}

// get a number greater than or equal to nions(element) for all elements
[[gnu::pure]] [[nodiscard]] inline auto get_max_nions() -> int { return maxnions; }

[[gnu::pure]] [[nodiscard]] inline auto elem_has_nlte_levels(const int element) -> bool {
  return globals::elements[element].has_nlte_levels;
}

[[gnu::pure]] [[nodiscard]] inline auto elem_has_nlte_levels_search(const int element) -> bool {
  for (int ion = 0; ion < get_nions(element); ion++) {
    for (int level = 1; level < get_nlevels(element, ion); level++) {
      if (is_nlte(element, ion, level)) {
        return true;
      }
    }
  }
  return false;
}

// Return the number of NLTE levels associated with with a specific ion given
// its elementindex and ionindex. Does not include the ground state or superlevel
[[gnu::pure]] [[nodiscard]] inline auto get_nlevels_excited_nlte(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].ions[ion].nlevels_excited_nlte;
}

// Get the number of autoionising levels for an ion
[[gnu::pure]] [[nodiscard]] inline auto get_nlevels_autoion(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].ions[ion].nlevels_autoion;
}

// get the number of levels that are not autoionising
[[gnu::pure]] [[nodiscard]] inline auto get_nlevels_nonautoion(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return get_nlevels(element, ion) - get_nlevels_autoion(element, ion);
}

// the number of downward autoionisation transitions from the specified level
[[gnu::pure]] [[nodiscard]] inline auto get_nautoiondowntrans(const int uniquelevelindex) -> int {
  return globals::alllevels.nautoiondowntrans[uniquelevelindex];
}

[[gnu::pure]] [[nodiscard]] inline auto get_nautoiondowntrans(const int element, const int ion, const int level)
    -> int {
  testmodeassert_valid_level(element, ion, level);
  return get_nautoiondowntrans(get_uniquelevelindex(element, ion, level));
}

// level has autoionising de-excitation transition
[[gnu::pure]] [[nodiscard]] inline auto level_isautoionising(const int element, const int ion, const int level)
    -> bool {
  assert_testmodeonly(get_nlevels_autoion(element, ion) < get_nlevels(element, ion));
  const bool is_autoionising = level >= get_nlevels_nonautoion(element, ion);
  assert_testmodeonly(is_autoionising == (get_nautoiondowntrans(element, ion, level) > 0));
  return is_autoionising;
}

// ion has NLTE levels, but this one is not NLTE and not autoionising => is in the superlevel
[[gnu::pure]] [[nodiscard]] inline auto level_isinsuperlevel(const int element, const int ion, const int level)
    -> bool {
  return (!is_nlte(element, ion, level) && level != 0 && (get_nlevels_excited_nlte(element, ion) > 0) &&
          !level_isautoionising(element, ion, level));
}

[[gnu::pure]] [[nodiscard]] inline auto get_nlevels_groundterm(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);
  return globals::elements[element].ions[ion].nlevels_groundterm;
}

// Get an index for an ionstage of an element that is unique for every ion of every element
[[gnu::pure]] [[nodiscard]] inline auto get_uniqueionindex(const int element, const int ion) -> int {
  testmodeassert_valid_ion(element, ion);

  const auto uniqueionindex = globals::elements[element].uniqueionindexstart + ion;
  assert_testmodeonly(uniqueionindex < get_includedions());

  return uniqueionindex;
}

// Index into the flattened per-cell-per-ion arrays (e.g. grid::ion_partfuncts_allcells) for an ion
// of a non-empty model grid cell.
[[gnu::pure]] [[nodiscard]] inline auto get_cellionindex(const ptrdiff_t nonemptymgi, const int uniqueionindex)
    -> ptrdiff_t {
  return (nonemptymgi * get_includedions()) + uniqueionindex;
}

[[gnu::pure]] [[nodiscard]] inline auto get_cellionindex(const ptrdiff_t nonemptymgi, const int element, const int ion)
    -> ptrdiff_t {
  return get_cellionindex(nonemptymgi, get_uniqueionindex(element, ion));
}

[[gnu::pure]] [[nodiscard]] inline auto get_ionfromuniqueionindex(const int uniqueionindex) -> std::tuple<int, int> {
  assert_testmodeonly(uniqueionindex < get_includedions());

  for (int element = 0; element < get_nelements(); element++) {
    if (get_nions(element) == 0) {
      continue;
    }
    const int ion = uniqueionindex - globals::elements[element].uniqueionindexstart;
    if (ion < get_nions(element)) {
      assert_testmodeonly(get_uniqueionindex(element, ion) == uniqueionindex);
      return {element, ion};
    }
  }
  assert_always(false);  // uniqueionindex too high to be valid
  return {-1, -1};
}

[[gnu::pure]] [[nodiscard]] inline auto ion_has_superlevel(const int element, const int ion) -> bool {
  testmodeassert_valid_ion(element, ion);
  return (get_nlevels(element, ion) > (get_nlevels_excited_nlte(element, ion) + get_nlevels_autoion(element, ion) + 1));
}

// the number of downward bound-bound transitions from a level, by unique level index or by (element, ion, level)
[[gnu::pure]] [[nodiscard]] inline auto get_ndowntrans(const int uniquelevelindex) -> int {
  return globals::alllevels.ndowntrans[uniquelevelindex];
}

[[gnu::pure]] [[nodiscard]] inline auto get_ndowntrans(const int element, const int ion, const int level) -> int {
  testmodeassert_valid_level(element, ion, level);
  return get_ndowntrans(get_uniquelevelindex(element, ion, level));
}

// index into globals::alltrans for first down transition from this level
[[gnu::pure]] [[nodiscard]] inline auto get_alltrans_startdown(const int uniquelevelindex) -> int {
  return globals::alllevels.alltrans_startdown[uniquelevelindex];
}

// index into globals::alltrans for first up transition from this level
[[gnu::pure]] [[nodiscard]] inline auto get_alltrans_startup(const int uniquelevelindex) -> int {
  return get_alltrans_startdown(uniquelevelindex) + get_ndowntrans(uniquelevelindex);
}

[[gnu::pure]] [[nodiscard]] inline auto get_alltrans_startup(const int element, const int ion, const int level) -> int {
  const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
  return get_alltrans_startup(uniquelevelindex);
}

// the number of upward bound-bound transitions from a level, by unique level index or by (element, ion, level)
[[gnu::pure]] [[nodiscard]] inline auto get_nuptrans(const int uniquelevelindex) -> int {
  return globals::alllevels.nuptrans[uniquelevelindex];
}

[[gnu::pure]] [[nodiscard]] inline auto get_nuptrans(const int element, const int ion, const int level) -> int {
  testmodeassert_valid_level(element, ion, level);
  return get_nuptrans(get_uniquelevelindex(element, ion, level));
}

// the number of upward autoionisation transitions from the specified level
[[gnu::pure]] [[nodiscard]] inline auto get_nautoionuptrans(const int uniquelevelindex) -> int {
  return globals::alllevels.nautoionuptrans[uniquelevelindex];
}

[[gnu::pure]] [[nodiscard]] inline auto get_nautoionuptrans(const int element, const int ion, const int level) -> int {
  testmodeassert_valid_level(element, ion, level);
  return globals::alllevels.nautoionuptrans[get_uniquelevelindex(element, ion, level)];
}

// Index of upperionlevel in the level's photoionisation target list, or -1 if it is not a target
[[gnu::pure]] [[nodiscard]] inline auto find_phixstargetindex(const int uniquelevelindex, const int upperionlevel)
    -> int {
  const auto nphixstargets = get_nphixstargets(uniquelevelindex);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
    if (upperionlevel == get_phixsupperlevel(uniquelevelindex, phixstargetindex)) {
      return phixstargetindex;
    }
  }
  return -1;
}

// Return the emissiontype index of the continuum associated to the given level, by unique level index or by
// (element, ion, level). Will be negative and ordered by element/ion/level/phixstargetindex. (NOTE! this is not
// an index into globals::allcont, which is ordered by ascending nu_edge)
[[gnu::pure]] [[nodiscard]] inline auto get_emtype_continuum(const int uniquelevelindex, const int phixstargetindex)
    -> int {
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(uniquelevelindex));
  return -1 - globals::alllevels.bflist_start[uniquelevelindex] - phixstargetindex;
}

[[gnu::pure]] [[nodiscard]] inline auto get_emtype_continuum(const int element, const int ion, const int level,
                                                             const int phixstargetindex) -> int {
  return get_emtype_continuum(get_uniquelevelindex(element, ion, level), phixstargetindex);
}

// Inverse of get_emtype_continuum(): decode a (negative) bound-free continuum emission type into its index into
// globals::bflist. Only valid for actual continuum emission types as returned by get_emtype_continuum(); the
// caller must first exclude the EMTYPE_NOTSET and EMTYPE_FREEFREE sentinel values (large negative numbers that
// would decode to out-of-range bflist indices).
[[nodiscard]] constexpr auto get_bflistindex_from_emtype_continuum(const int emissiontype) -> int {
  return -1 - emissiontype;
}

// Return the photionisation threshold energy [erg]
[[gnu::pure]] [[nodiscard]] inline auto get_phixs_threshold(const int element, const int ion, const int level,
                                                            const int phixstargetindex) -> double {
  testmodeassert_valid_level(element, ion, level);
  const int uniquelevelindex = get_uniquelevelindex(element, ion, level);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(uniquelevelindex));
  const int upperlevel = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
  const double E_threshold = epsilon(element, ion + 1, upperlevel) - epsilon(uniquelevelindex);
  return E_threshold;
}

#endif  // ATOMIC_H
