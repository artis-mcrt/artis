#pragma once
#ifndef ATOMIC_H
#define ATOMIC_H

#include <algorithm>
#include <array>
#include <tuple>

#include "grid.h"
#include "ltepop.h"
#include "sn3d.h"

// highest number of ions for any element
inline int maxnions = 0;

// number of ions of any element
inline int includedions = 0;

// total number of levels of any element
inline int includedlevels = 0;

inline double
    last_phixs_nuovernuedge;  // last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge
inline std::array<bool, 3> phixs_file_version_exists;  // first value in this array is not used but exists so the
                                                       // indexes match those of the phixsdata_filenames array
constexpr std::array<const char *, 3> phixsdata_filenames = {"version0ignore", "phixsdata.txt", "phixsdata_v2.txt"};

[[nodiscard]] inline auto get_nelements() -> int { return static_cast<int>(globals::elements.size()); }

inline auto get_nnion_tot(int modelgridindex) -> double
// total density of nuclei
{
  double nntot = 0.;
  for (int element = 0; element < get_nelements(); element++) {
    nntot += grid::get_elem_numberdens(modelgridindex, element);
  }

  return nntot;
}

inline auto get_nions(const int element) -> int
/// Returns the number of ions associated with a specific element given by
/// its elementindex.
{
  assert_testmodeonly(element < get_nelements());
  return globals::elements[element].nions;
}

inline auto get_nlevels(const int element, const int ion) -> int
/// Returns the number of levels associated with with a specific ion given
/// its elementindex and ionindex.
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].nlevels;
}

[[nodiscard]] inline auto epsilon(const int element, const int ion, const int level) -> double
/// Returns the energy of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return globals::elements[element].ions[ion].levels[level].epsilon;
}

[[nodiscard]] inline auto get_ionstage(const int element, const int ion) -> int
/// Returns the ionisation stage of an ion specified by its elementindex and
/// ionindex.
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].ionstage;
}

[[nodiscard]] inline auto get_ionisinglevels(const int element, const int ion) -> int
/// Returns the number of levels associated with an ion that
/// have energies below the ionisation threshold.
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].ionisinglevels;
}

inline auto get_nphixstargets(const int element, const int ion, const int level) -> int
/// Returns the number of target states for photoionization of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(level < get_nlevels(element, ion));
  const int nions = get_nions(element);
  assert_testmodeonly(ion < nions);
  if ((ion < nions - 1) && (level < get_ionisinglevels(element, ion))) {
    return globals::elements[element].ions[ion].levels[level].nphixstargets;
  }
  return 0;
}

[[nodiscard]] inline auto get_phixsupperlevel(const int element, const int ion, const int level,
                                              const int phixstargetindex) -> int
/// Returns the level index of a target state for photoionization of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(element, ion, level));

  return globals::elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].levelindex;
}

[[nodiscard]] inline auto get_phixsprobability(const int element, const int ion, const int level,
                                               const int phixstargetindex) -> double
/// Returns the probability of a target state for photoionization of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(element, ion, level));

  return globals::elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].probability;
}

[[nodiscard]] inline auto einstein_spontaneous_emission(const int lineindex) -> double
// double einstein_spontaneous_emission(int element, int ion, int upper, int lower)
/// reads A_ul from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... |
/// (A_upper,1; f_upper,1)
{
  return globals::linelist[lineindex].einstein_A;
}

[[nodiscard]] inline auto stat_weight(const int element, const int ion, const int level) -> double
/// Returns the statistical weight of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return globals::elements[element].ions[ion].levels[level].stat_weight;
}

[[nodiscard]] inline auto get_maxrecombininglevel(const int element, const int ion) -> int
/// Returns the number of bf-continua associated with ion ion of element element.
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].maxrecombininglevel;
}

[[nodiscard]] inline auto photoionization_crosssection_fromtable(const float *const photoion_xs, const double nu_edge,
                                                                 const double nu) -> double
/// Calculates the photoionisation cross-section at frequency nu out of the atomic data.
/// Input: - edge frequency nu_edge of the desired bf-continuum
///        - nu
{
  // if (nu < nu_edge || nu > nu_edge * 1.05)
  //   return 0;
  // else
  //   return 1.;
  // return 1. * pow(nu_edge / nu, 3);

  float sigma_bf = 0.;

  if (phixs_file_version_exists[1] && !phixs_file_version_exists[2]) {
    // classic mode: no interpolation
    if (nu == nu_edge) {
      sigma_bf = photoion_xs[0];
    } else if (nu <= nu_edge * (1 + globals::NPHIXSNUINCREMENT * globals::NPHIXSPOINTS)) {
      const int i = static_cast<int>(floor((nu - nu_edge) / (globals::NPHIXSNUINCREMENT * nu_edge)));
      sigma_bf = photoion_xs[i];
    } else {
      /// use a parameterization of sigma_bf by the Kramers formula
      /// which anchor point should we take ??? the cross-section at the edge or at the highest grid point ???
      /// so far the highest grid point, otherwise the cross-section is not continuous
      sigma_bf = photoion_xs[globals::NPHIXSPOINTS - 1] *
                 pow(nu_edge * (1 + globals::NPHIXSNUINCREMENT * globals::NPHIXSPOINTS) / nu, 3);
    }
    return sigma_bf;
  }

  const double ireal = (nu / nu_edge - 1.0) / globals::NPHIXSNUINCREMENT;
  const int i = floor(ireal);

  if (i < 0) {
    sigma_bf = 0.;
    // printout("[warning] photoionization_crosssection was called with nu=%g < nu_edge=%g\n",nu,nu_edge);
    // printout("[warning]   element %d, ion %d, level %d, epsilon %g, ionpot
    // %g\n",element,ion,level,epsilon(element,ion,level),elements[element].ions[ion].ionpot); printout("[warning]
    // element %d, ion+1 %d, level %d epsilon %g, ionpot
    // %g\n",element,ion+1,0,epsilon(element,ion+1,0),elements[element].ions[ion].ionpot); printout("[warning]
    // photoionization_crosssection %g\n",sigma_bf); abort();
  } else if (i < globals::NPHIXSPOINTS - 1) {
    // sigma_bf = globals::elements[element].ions[ion].levels[level].photoion_xs[i];

    const double sigma_bf_a = photoion_xs[i];
    const double sigma_bf_b = photoion_xs[i + 1];
    const double factor_b = ireal - i;
    sigma_bf = ((1. - factor_b) * sigma_bf_a) + (factor_b * sigma_bf_b);
  } else {
    /// use a parameterization of sigma_bf by the Kramers formula
    /// which anchor point should we take ??? the cross-section at the edge or at the highest grid point ???
    /// so far the highest grid point, otherwise the cross-section is not continuous
    const double nu_max_phixs = nu_edge * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table
    sigma_bf = photoion_xs[globals::NPHIXSPOINTS - 1] * pow(nu_max_phixs / nu, 3);
  }

  // if (sigma_bf < 0)
  // {
  //   printout("[warning] photoionization_crosssection returns negative cross-section %g\n",sigma_bf);
  //   printout("[warning]   nu=%g,  nu_edge=%g\n",nu,nu_edge);
  //   printout("[warning]   xs@edge=%g,
  //   xs@maxfreq\n",elements[element].ions[ion].levels[level].photoion_xs[0],elements[element].ions[ion].levels[level].photoion_xs[NPHIXSPOINTS-1]);
  //   printout("[warning]   element %d, ion %d, level %d, epsilon %g, ionpot
  //   %g\n",element,ion,level,epsilon(element,ion,level),elements[element].ions[ion].ionpot);
  // }

  return sigma_bf;
}

[[nodiscard]] inline auto photoionization_crosssection(const int element, const int ion, const int level,
                                                       const double nu_edge, const double nu) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return photoionization_crosssection_fromtable(globals::elements[element].ions[ion].levels[level].photoion_xs, nu_edge,
                                                nu);
}

[[nodiscard]] inline auto get_tau_sobolev(const int modelgridindex, const int lineindex, const double t_current,
                                          bool sub_updown) -> double {
  const int element = globals::linelist[lineindex].elementindex;
  const int ion = globals::linelist[lineindex].ionindex;
  const int lower = globals::linelist[lineindex].lowerlevelindex;
  const int upper = globals::linelist[lineindex].upperlevelindex;

  const double n_l = get_levelpop(modelgridindex, element, ion, lower);

  const double nu_trans = (epsilon(element, ion, upper) - epsilon(element, ion, lower)) / H;
  const double A_ul = einstein_spontaneous_emission(lineindex);
  const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
  const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

  if (sub_updown) {
    const double n_u = get_levelpop(modelgridindex, element, ion, upper);
    return std::max((B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current, 0.);
  }
  return std::max(B_lu * n_l * HCLIGHTOVERFOURPI * t_current, 0.);
}

inline void set_nelements(const int nelements_in) { globals::elements.resize(nelements_in); }

inline auto get_atomicnumber(const int element) -> int
/// Returns the atomic number associated with a given elementindex.
{
  assert_testmodeonly(element >= 0);
  assert_testmodeonly(element < get_nelements());
  return globals::elements[element].anumber;
}

[[nodiscard]] inline auto is_nlte(const int element, const int ion, const int level) -> bool
// Returns true if (element,ion,level) is to be treated in nlte.
// (note this function returns true for the ground state,
//  although it is stored separately from the excited NLTE states)
{
  return LEVEL_IS_NLTE(get_atomicnumber(element), get_ionstage(element, ion),
                       level);  // defined in artisoptions.h
}

inline auto get_elementindex(const int Z) -> int
/// Returns the elementindex associated with a given atomic number.
/// If there is no element with the given atomic number in the atomic data
/// a negative value is returned to flag this event.
{
  const auto elem =
      std::ranges::find_if(globals::elements, [Z](const Element &element) { return element.anumber == Z; });
  if (elem != globals::elements.end()) {
    return std::distance(globals::elements.begin(), elem);
  }

  // printout("[debug] get_elementindex: element Z=%d was not found in atomic data ... skip readin of cross sections
  // for this element\n",Z); printout("[fatal] get_elementindex: element Z=%d was not found in atomic data ...
  // abort\n"); abort();;
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

inline auto get_includedions() -> int
// returns the number of ions of all elements combined
{
  return includedions;
}

inline auto get_includedlevels() -> int
// returns the number of ions of all elements combined
{
  return includedlevels;
}

[[nodiscard]] inline auto get_max_nions() -> int {
  // number greater than or equal to nions(element) for all elements
  return maxnions;
}

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

[[nodiscard]] inline auto get_nlevels_nlte(const int element, const int ion) -> int
// Returns the number of NLTE levels associated with with a specific ion given
// its elementindex and ionindex. Includes the superlevel if there is one but does not include the ground state
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].nlevels_nlte;
}

[[nodiscard]] inline auto level_isinsuperlevel(const int element, const int ion, const int level) -> bool
// ion has NLTE levels, but this one is not NLTE => is in the superlevel
{
  return (!is_nlte(element, ion, level) && level != 0 && (get_nlevels_nlte(element, ion) > 0));
}

[[nodiscard]] inline auto get_nlevels_groundterm(const int element, const int ion) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].nlevels_groundterm;
}

[[nodiscard]] inline auto get_uniqueionindex(const int element, const int ion) -> int
// Get an index for an ionstage of an element that is unique for every ion of every element
{
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
}

[[nodiscard]] inline auto get_uniquelevelindex(const int element, const int ion, const int level) -> int
// Get an index for level of an ionstage of an element that is unique across every ion of every element
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  const auto uniquelevelindex = globals::elements[element].ions[ion].uniquelevelindexstart + level;
  assert_testmodeonly(uniquelevelindex < get_includedlevels());

  return uniquelevelindex;
}

[[nodiscard]] inline auto get_levelfromuniquelevelindex(const int alllevelsindex) -> std::tuple<int, int, int>
// inverse of get_uniquelevelindex(). get the element/ion/level from a unique level index
{
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      if (get_nlevels(element, ion) == 0) {
        continue;
      }
      const int level = alllevelsindex - globals::elements[element].ions[ion].uniquelevelindexstart;
      if (level < get_nlevels(element, ion)) {
        assert_testmodeonly(get_uniquelevelindex(element, ion, level) == alllevelsindex);
        return {element, ion, level};
      }
    }
  }
  assert_always(false);  // alllevelsindex too high to be valid
}

[[nodiscard]] inline auto ion_has_superlevel(const int element, const int ion) -> bool {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return (get_nlevels(element, ion) > get_nlevels_nlte(element, ion) + 1);
}

[[nodiscard]] inline auto get_ndowntrans(const int element, const int ion, const int level) -> int
// the number of downward bound-bound transitions from the specified level
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return globals::elements[element].ions[ion].levels[level].ndowntrans;
}

[[nodiscard]] inline auto get_nuptrans(const int element, const int ion, const int level) -> int
// the number of upward bound-bound transitions from the specified level
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return globals::elements[element].ions[ion].levels[level].nuptrans;
}

inline void set_ndowntrans(const int element, const int ion, const int level, const int ndowntrans)
// the number of downward bound-bound transitions from the specified level
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  globals::elements[element].ions[ion].levels[level].ndowntrans = ndowntrans;
}

inline void set_nuptrans(const int element, const int ion, const int level, const int nuptrans)
// the number of upward bound-bound transitions from the specified level
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  globals::elements[element].ions[ion].levels[level].nuptrans = nuptrans;
}

[[nodiscard]] inline auto get_phixtargetindex(const int element, const int ion, const int level,
                                              const int upperionlevel) -> int {
  const auto nphixstargets = get_nphixstargets(element, ion, level);
  for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
    if (upperionlevel == get_phixsupperlevel(element, ion, level, phixstargetindex)) {
      return phixstargetindex;
    }
  }
  printout("Could not find phixstargetindex\n");
  std::abort();
  return -1;
}

[[nodiscard]] inline auto get_emtype_continuum(const int element, const int ion, const int level,
                                               const int upperionlevel) -> int
// Returns the emissiontype index of the continuum associated to the given level. Will be negative and ordered by
// element/ion/level/phixstargetindex
{
  const int phixstargetindex = get_phixtargetindex(element, ion, level, upperionlevel);
  return globals::elements[element].ions[ion].levels[level].cont_index - phixstargetindex;
}

[[nodiscard]] inline auto get_phixs_threshold(const int element, const int ion, const int level,
                                              const int phixstargetindex) -> double
/// Returns the energy of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  assert_testmodeonly(phixstargetindex < get_nphixstargets(element, ion, level));
  const int upperlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
  const double E_threshold = epsilon(element, ion + 1, upperlevel) - epsilon(element, ion, level);
  return E_threshold;
}

#endif  // ATOMIC_H
