#ifndef ATOMIC_H
#define ATOMIC_H

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <string>
#include <tuple>

#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "sn3d.h"

// highest number of ions for any element
inline int maxnions = 0;

// number of ions of any element
inline int includedions = 0;

// total number of levels of any element
inline int includedlevels = 0;

// last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge
inline double last_phixs_nuovernuedge;

// first value in this array is not used but exists so the indexes match those of the phixsdata_filenames array
inline std::array<bool, 3> phixs_file_version_exists;

const std::array<const std::string, 3> phixsdata_filenames = {"IGNORE", "phixsdata.txt", "phixsdata_v2.txt"};

// return the number of ions of all elements combined
inline auto get_includedlevels() -> int { return includedlevels; }

[[nodiscard]] __host__ __device__ inline auto get_nelements() -> int {
  return static_cast<int>(globals::elements.size());
}

// total density of nuclei
__host__ __device__ inline auto get_nnion_tot(int nonemptymgi) -> double {
  double nntot = 0.;
  for (int element = 0; element < get_nelements(); element++) {
    nntot += grid::get_elem_numberdens(nonemptymgi, element);
  }

  return nntot;
}

// Return the number of ions associated with a specific element given by its elementindex.
inline auto get_nions(const int element) -> int {
  assert_testmodeonly(element < get_nelements());
  return globals::elements[element].nions;
}

// Return the number of levels associated with a specific ion given its elementindex and ionindex.
__host__ __device__ inline auto get_nlevels(const int element, const int ion) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].nlevels;
}

// Get an index for level of an ionstage of an element that is unique across every ion of every element
[[nodiscard]] inline auto get_uniquelevelindex(const int element, const int ion, const int level) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  const auto uniquelevelindex = globals::elements[element].ions[ion].uniquelevelindexstart + level;
  assert_testmodeonly(uniquelevelindex < get_includedlevels());

  return uniquelevelindex;
}

// Return the ionisation stage of an ion specified by its elementindex and ionindex.
[[nodiscard]] __host__ __device__ inline auto get_ionstage(const int element, const int ion) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].ionstage;
}

// Return the number of levels associated with an ion that have energies below the ionisation threshold.
[[nodiscard]] __host__ __device__ inline auto get_nlevels_ionising(const int element, const int ion) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].nlevels_ionising;
}

// Returns the number of target states for photoionization of (element,ion,level).
__host__ __device__ inline auto get_nphixstargets(const int uniquelevelindex) -> int {
  return globals::alllevels.nphixstargets[uniquelevelindex];
}

// Returns the number of target states for photoionization of (element,ion,level).
__host__ __device__ inline auto get_nphixstargets(const int element, const int ion, const int level) -> int {
  const auto nphixstargets = get_nphixstargets(get_uniquelevelindex(element, ion, level));
  assert_testmodeonly(nphixstargets == 0 ||
                      ((ion < (get_nions(element) - 1)) && (level < get_nlevels_ionising(element, ion))));
  return nphixstargets;
}

// Return the index into the allphixstargets arrays for a target state for photoionization of (element,ion,level).
[[nodiscard]] __host__ __device__ inline auto get_allphixstargetindex(const int uniquelevelindex,
                                                                      const int phixstargetindex) -> int {
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(uniquelevelindex));

  return globals::alllevels.phixstargetstart[uniquelevelindex] + phixstargetindex;
}
// Return the level index of a target state for photoionization of (element,ion,level).
[[nodiscard]] __host__ __device__ inline auto get_phixsupperlevel(const int uniquelevelindex,
                                                                  const int phixstargetindex) -> int {
  return globals::allphixstargets_levelindex[get_allphixstargetindex(uniquelevelindex, phixstargetindex)];
}

// Return the level index of a target state for photoionization of (element,ion,level).
[[nodiscard]] __host__ __device__ inline auto get_phixsupperlevel(const int element, const int ion, const int level,
                                                                  const int phixstargetindex) -> int {
  return get_phixsupperlevel(get_uniquelevelindex(element, ion, level), phixstargetindex);
}

// Return the probability of a target state for photoionization of (element,ion,level).
[[nodiscard]] __host__ __device__ inline auto get_phixsprobability(const int uniquelevelindex,
                                                                   const int phixstargetindex) -> double {
  return globals::allphixstargets_probability[get_allphixstargetindex(uniquelevelindex, phixstargetindex)];
}

// Return the probability of a target state for photoionization of (element,ion,level).
[[nodiscard]] __host__ __device__ inline auto get_phixsprobability(const int element, const int ion, const int level,
                                                                   const int phixstargetindex) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));

  return get_phixsprobability(get_uniquelevelindex(element, ion, level), phixstargetindex);
}

// Return the number of bf-continua associated with ion ion of element element.
[[nodiscard]] inline auto get_maxrecombininglevel(const int element, const int ion) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].maxrecombininglevel;
}

[[nodiscard]] inline __host__ __device__ auto get_phixs_table(const int uniquelevelindexl) -> float * {
  const auto phixsstart = globals::alllevels.phixsstart[uniquelevelindexl];
  assert_testmodeonly(phixsstart >= 0);
  return globals::allphixs.data() + (phixsstart * globals::NPHIXSPOINTS);
}

[[nodiscard]] inline __host__ __device__ auto get_phixs_table(const int element, const int ion, const int level)
    -> float * {
  return get_phixs_table(get_uniquelevelindex(element, ion, level));
}

// Calculate the photoionisation cross-section at frequency nu out of the atomic data.
[[nodiscard]] inline auto photoionization_crosssection_fromtable(const float *const photoion_xs, const double nu_edge,
                                                                 const double nu) -> float {
  // if (nu < nu_edge || nu > nu_edge * 1.05)
  //   return 0;
  // else
  //   return 1.;
  // return 1. * pow(nu_edge / nu, 3);

  float sigma_bf = 0.;

  if constexpr (PHIXS_CLASSIC_NO_INTERPOLATION) {
    // classic mode: no interpolation
    if (nu == nu_edge) {
      sigma_bf = photoion_xs[0];
    } else if (nu <= nu_edge * (1 + globals::NPHIXSNUINCREMENT * globals::NPHIXSPOINTS)) {
      const int i = static_cast<int>((nu - nu_edge) / (globals::NPHIXSNUINCREMENT * nu_edge));
      sigma_bf = photoion_xs[i];
    } else {
      // use a parameterization of sigma_bf by the Kramers formula
      // which anchor point should we take ??? the cross-section at the edge or at the highest grid point ???
      // so far the highest grid point, otherwise the cross-section is not continuous
      sigma_bf = static_cast<float>(photoion_xs[globals::NPHIXSPOINTS - 1] *
                                    pow(nu_edge * (1 + globals::NPHIXSNUINCREMENT * globals::NPHIXSPOINTS) / nu, 3));
    }
    return sigma_bf;
  }

  const double ireal = (nu / nu_edge - 1.0) / globals::NPHIXSNUINCREMENT;
  const int i = static_cast<int>(ireal);

  if (i < 0) {
    sigma_bf = 0.;
  } else if (i < globals::NPHIXSPOINTS - 1) {
    const double sigma_bf_a = photoion_xs[i];
    const double sigma_bf_b = photoion_xs[i + 1];
    const double factor_b = ireal - i;
    sigma_bf = static_cast<float>(((1. - factor_b) * sigma_bf_a) + (factor_b * sigma_bf_b));
  } else {
    // use a parameterization of sigma_bf by the Kramers formula
    // which anchor point should we take ??? the cross-section at the edge or at the highest grid point ???
    // so far the highest grid point, otherwise the cross-section is not continuous
    const double nu_max_phixs = nu_edge * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table
    sigma_bf = static_cast<float>(photoion_xs[globals::NPHIXSPOINTS - 1] * pow(nu_max_phixs / nu, 3));
  }

  return sigma_bf;
}

// inverse of get_uniquelevelindex(). get the element/ion/level from a unique level index
[[nodiscard]] inline auto get_levelfromuniquelevelindex(const int uniquelevelindex) -> std::tuple<int, int, int> {
  assert_testmodeonly(uniquelevelindex < get_includedlevels());
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      if (get_nlevels(element, ion) == 0) {
        continue;
      }
      const int level = uniquelevelindex - globals::elements[element].ions[ion].uniquelevelindexstart;
      if (level >= 0 && level < get_nlevels(element, ion)) {
        assert_testmodeonly(get_uniquelevelindex(element, ion, level) == uniquelevelindex);
        return {element, ion, level};
      }
    }
  }
  assert_always(false);  // uniquelevelindex too high to be valid
  return {-1, -1, -1};
}
// Return the statistical weight of (element,ion,level).
[[nodiscard]] __host__ __device__ inline auto stat_weight(const int uniquelevelindex) -> double {
  return globals::alllevels.statweight[uniquelevelindex];
}

// Return the statistical weight of (element,ion,level).
[[nodiscard]] __host__ __device__ inline auto stat_weight(const int element, const int ion, const int level) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return stat_weight(get_uniquelevelindex(element, ion, level));
}

// Return the energy of (element,ion,level).
[[nodiscard]] __host__ __device__ inline auto epsilon(const int uniquelevelindex) -> double {
  return globals::alllevels.epsilon[uniquelevelindex];
}

[[nodiscard]] __host__ __device__ inline auto epsilon(const int element, const int ion, const int level) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return epsilon(get_uniquelevelindex(element, ion, level));
}

[[nodiscard]] inline auto get_tau_sobolev(const int nonemptymgi, const int lineindex, const double t_current)
    -> double {
  const auto &line = globals::linelist[lineindex];
  const int element = line.elementindex;
  const int ion = line.ionindex;
  const int lower = line.lowerlevelindex;
  const int upper = line.upperlevelindex;
  const auto ionuniquelevelindexstart = globals::elements[element].ions[ion].uniquelevelindexstart;

  const double n_l = get_levelpop(nonemptymgi, ionuniquelevelindexstart + lower);

  const double A_ul = line.einstein_A;
  const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(line.nu, 3) * A_ul;
  const double B_lu =
      stat_weight(ionuniquelevelindexstart + upper) / stat_weight(ionuniquelevelindexstart + lower) * B_ul;

  return std::max(B_lu * n_l * HCLIGHTOVERFOURPI * t_current, 0.);
}

[[nodiscard]] inline auto get_tau_sobolev_subupdown(const int nonemptymgi, const TransitionLine &line,
                                                    const double t_current) -> double {
  const int element = line.elementindex;
  const int ion = line.ionindex;
  const int lower = line.lowerlevelindex;
  const int upper = line.upperlevelindex;

  const auto ionuniquelevelindexstart = globals::elements[element].ions[ion].uniquelevelindexstart;
  const double n_l = get_levelpop(nonemptymgi, ionuniquelevelindexstart + lower);

  const double A_ul = line.einstein_A;
  const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(line.nu, 3) * A_ul;
  const double B_lu =
      stat_weight(ionuniquelevelindexstart + upper) / stat_weight(ionuniquelevelindexstart + lower) * B_ul;

  const double n_u = get_levelpop(nonemptymgi, ionuniquelevelindexstart + upper);
  return std::max((B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current, 0.);
}

// Returns the atomic number associated with a given elementindex.
inline auto get_atomicnumber(const int element) -> int {
  assert_testmodeonly(element >= 0);
  assert_testmodeonly(element < get_nelements());
  return globals::elements[element].anumber;
}

// Returns true if (element,ion,level) is to be treated in nlte.
// (note this function returns true for the ground state,
//  although it is stored separately from the excited NLTE states)
[[nodiscard]] inline auto is_nlte(const int element, const int ion, const int level) -> bool {
  return LEVEL_IS_NLTE(get_atomicnumber(element), get_ionstage(element, ion),
                       level);  // defined in artisoptions.h
}

// Return the elementindex associated with a given atomic number.
// If there is no element with the given atomic number in the atomic data,
// then return a negative value
inline auto get_elementindex(const int Z) -> int {
  const auto elem =
      std::ranges::find_if(globals::elements, [Z](const Element &element) { return element.anumber == Z; });
  if (elem != globals::elements.end()) {
    return static_cast<int>(elem - globals::elements.begin());
  }

  return -100;
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
inline auto get_includedions() -> int {
  assert_testmodeonly(includedions > 0);
  return includedions;
}

// get a number greater than or equal to nions(element) for all elements
[[nodiscard]] inline auto get_max_nions() -> int { return maxnions; }

[[nodiscard]] inline auto elem_has_nlte_levels(const int element) -> bool {
  return globals::elements[element].has_nlte_levels;
}

[[nodiscard]] inline auto elem_has_nlte_levels_search(const int element) -> bool {
  for (int ion = 0; ion < get_nions(element); ion++) {
    for (int level = 1; level < get_nlevels(element, ion); level++) {
      if (is_nlte(element, ion, level)) {
        return true;
      }
    }
  }
  return false;
}

// Returns the number of NLTE levels associated with with a specific ion given
// its elementindex and ionindex. Includes the superlevel if there is one but does not include the ground state
[[nodiscard]] inline auto get_nlevels_nlte(const int element, const int ion) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].nlevels_nlte;
}

// ion has NLTE levels, but this one is not NLTE => is in the superlevel
[[nodiscard]] inline auto level_isinsuperlevel(const int element, const int ion, const int level) -> bool {
  return (!is_nlte(element, ion, level) && level != 0 && (get_nlevels_nlte(element, ion) > 0));
}

[[nodiscard]] inline auto get_nlevels_groundterm(const int element, const int ion) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].nlevels_groundterm;
}

// Get an index for an ionstage of an element that is unique for every ion of every element
[[nodiscard]] inline auto get_uniqueionindex(const int element, const int ion) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));

  const auto uniqueionindex = globals::elements[element].uniqueionindexstart + ion;
  assert_testmodeonly(uniqueionindex < get_includedions());

  return uniqueionindex;
}

[[nodiscard]] inline auto get_ionfromuniqueionindex(const int allionsindex) -> std::tuple<int, int> {
  assert_testmodeonly(allionsindex < get_includedions());

  for (int element = 0; element < get_nelements(); element++) {
    if (get_nions(element) == 0) {
      continue;
    }
    const int ion = allionsindex - globals::elements[element].uniqueionindexstart;
    if (ion < get_nions(element)) {
      assert_testmodeonly(get_uniqueionindex(element, ion) == allionsindex);
      return {element, ion};
    }
  }
  assert_always(false);  // allionsindex too high to be valid
  return {-1, -1};
}

[[nodiscard]] inline auto ion_has_superlevel(const int element, const int ion) -> bool {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return (get_nlevels(element, ion) > get_nlevels_nlte(element, ion) + 1);
}

// the number of downward bound-bound transitions from the specified level
[[nodiscard]] inline auto get_ndowntrans(const int uniquelevelindex) -> int {
  return globals::alllevels.ndowntrans[uniquelevelindex];
}

// the number of downward bound-bound transitions from the specified level
[[nodiscard]] inline auto get_ndowntrans(const int element, const int ion, const int level) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return get_ndowntrans(get_uniquelevelindex(element, ion, level));
}

// index into globals::alltrans for first down transition from this level
[[nodiscard]] inline auto get_alltrans_startdown(const int uniquelevelindex) -> int {
  return globals::alllevels.alltrans_startdown[uniquelevelindex];
}

// index into globals::alltrans for first up transition from this level
[[nodiscard]] inline auto get_alltrans_startup(const int uniquelevelindex) -> int {
  return get_alltrans_startdown(uniquelevelindex) + get_ndowntrans(uniquelevelindex);
}

[[nodiscard]] inline auto get_alltrans_startup(const int element, const int ion, const int level) -> int {
  const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
  return get_alltrans_startup(uniquelevelindex);
}

// the number of upward bound-bound transitions from the specified level
[[nodiscard]] inline auto get_nuptrans(const int uniquelevelindex) -> int {
  return globals::alllevels.nuptrans[uniquelevelindex];
}

// the number of upward bound-bound transitions from the specified level
[[nodiscard]] inline auto get_nuptrans(const int element, const int ion, const int level) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return get_nuptrans(get_uniquelevelindex(element, ion, level));
}

[[nodiscard]] inline auto get_phixtargetindex(const int uniquelevelindex, const int upperionlevel) -> int {
  const auto nphixstargets = get_nphixstargets(uniquelevelindex);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
    if (upperionlevel == get_phixsupperlevel(uniquelevelindex, phixstargetindex)) {
      return phixstargetindex;
    }
  }
  assert_testmodeonly(false);
  if constexpr (!TESTMODE) {
    __builtin_unreachable();
  }
  return -1;
}

// Return the emissiontype index of the continuum associated to the given level. Will be negative and ordered by
// element/ion/level/phixstargetindex. (NOTE! this is not an index into globals::allcont, which is ordered by ascending
// nu_edge)
[[nodiscard]] inline auto get_emtype_continuum(const int element, const int ion, const int level,
                                               const int upperionlevel) -> int {
  const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
  const int phixstargetindex = get_phixtargetindex(uniquelevelindex, upperionlevel);
  return -1 - globals::alllevels.bflist_start[uniquelevelindex] - phixstargetindex;
}

// Return the photionisation threshold energy [erg]
[[nodiscard]] inline auto get_phixs_threshold(const int uniquelevelindex, const int phixstargetindex) -> double {
  assert_testmodeonly(phixstargetindex < get_nphixstargets(uniquelevelindex));
  const int upperlevel = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
  const auto [element, ion, _] = get_levelfromuniquelevelindex(uniquelevelindex);
  const double E_threshold = epsilon(element, ion + 1, upperlevel) - epsilon(uniquelevelindex);
  return E_threshold;
}

// Return the photionisation threshold energy [erg]
[[nodiscard]] inline auto get_phixs_threshold(const int element, const int ion, const int level,
                                              const int phixstargetindex) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  const int uniquelevelindex = get_uniquelevelindex(element, ion, level);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(uniquelevelindex));
  const int upperlevel = get_phixsupperlevel(uniquelevelindex, phixstargetindex);
  const double E_threshold = epsilon(element, ion + 1, upperlevel) - epsilon(uniquelevelindex);
  return E_threshold;
}

[[nodiscard]] inline auto get_bflutindex(const int temperatureindex, const int uniquelevelindex,
                                         const int phixstargetindex) -> int {
  const int contindex = globals::alllevels.bflist_start[uniquelevelindex] + phixstargetindex;
  const int bflutindex = (temperatureindex * globals::nbfcontinua) + contindex;
  assert_testmodeonly(bflutindex >= 0);
  assert_testmodeonly(bflutindex <= TABLESIZE * globals::nbfcontinua);
  return bflutindex;
}

[[nodiscard]] inline auto get_bflutindex(const int temperatureindex, const int element, const int ion, const int level,
                                         const int phixstargetindex) -> int {
  return get_bflutindex(temperatureindex, get_uniquelevelindex(element, ion, level), phixstargetindex);
}

#endif  // ATOMIC_H
