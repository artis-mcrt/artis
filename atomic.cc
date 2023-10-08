#include "atomic.h"

#include <cmath>

#include "artisoptions.h"
#include "grid.h"
#include "ltepop.h"
#include "sn3d.h"

double last_phixs_nuovernuedge =
    -1;                // last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge
int nelements = 0;     // total number of elements included in the simulation
int maxnions = 0;      // highest number of ions for any element
int includedions = 0;  // number of ions of any element
std::array<bool, 3> phixs_file_version_exists;

auto get_continuumindex_phixstargetindex(const int element, const int ion, const int level, const int phixstargetindex)
    -> int
/// Returns the index of the continuum associated to the given level.
{
  return globals::elements[element].ions[ion].levels[level].cont_index - phixstargetindex;
}

auto get_phixtargetindex(const int element, const int ion, const int level, const int upperionlevel) -> int {
  for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++) {
    if (upperionlevel == get_phixsupperlevel(element, ion, level, phixstargetindex)) {
      return phixstargetindex;
    }
  }
  printout("Could not find phixstargetindex\n");
  abort();
  return -1;
}

auto get_continuumindex(const int element, const int ion, const int level, const int upperionlevel) -> int
/// Returns the index of the continuum associated to the given level.
{
  const int phixstargetindex = get_phixtargetindex(element, ion, level, upperionlevel);
  return get_continuumindex_phixstargetindex(element, ion, level, phixstargetindex);
}

auto get_tau_sobolev(const int modelgridindex, const int lineindex, const double t_current) -> double {
  const int element = globals::linelist[lineindex].elementindex;
  const int ion = globals::linelist[lineindex].ionindex;
  const int lower = globals::linelist[lineindex].lowerlevelindex;
  const int upper = globals::linelist[lineindex].upperlevelindex;

  const double n_l = get_levelpop(modelgridindex, element, ion, lower);
  const double n_u = get_levelpop(modelgridindex, element, ion, upper);

  const double nu_trans = (epsilon(element, ion, upper) - epsilon(element, ion, lower)) / H;
  const double A_ul = einstein_spontaneous_emission(lineindex);
  const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
  const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

  const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;
  return tau_sobolev;
}

auto get_tot_nions(int modelgridindex) -> double
// total density of nuclei
{
  double nntot = 0.;
  for (int element = 0; element < get_nelements(); element++) {
    nntot += grid::get_elem_numberdens(modelgridindex, element);
  }

  return nntot;
}

auto is_nlte(const int element, const int ion, const int level) -> bool
// Returns true if (element,ion,level) is to be treated in nlte.
// (note this function returns true for the ground state,
//  although it is stored separately from the excited NLTE states)
{
  if (!NLTE_POPS_ON) {
    return false;
  }
  return LEVEL_IS_NLTE(get_atomicnumber(element), get_ionstage(element, ion),
                       level);  // defined in artisoptions.h
}

auto level_isinsuperlevel(const int element, const int ion, const int level) -> bool
// ion has NLTE levels, but this one is not NLTE => is in the superlevel
{
  return (NLTE_POPS_ON && !is_nlte(element, ion, level) && level != 0 && (get_nlevels_nlte(element, ion) > 0));
}

auto photoionization_crosssection_fromtable(const float *const photoion_xs, const double nu_edge, const double nu)
    -> double
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
    sigma_bf = 0.0;
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

void set_nelements(const int nelements_in) { nelements = nelements_in; }

auto get_nelements() -> int { return nelements; }

auto get_atomicnumber(const int element) -> int
/// Returns the atomic number associated with a given elementindex.
{
  assert_testmodeonly(element >= 0);
  assert_testmodeonly(element < get_nelements());
  return globals::elements[element].anumber;
}

auto get_elementindex(const int Z) -> int
/// Returns the elementindex associated with a given atomic number.
/// If there is no element with the given atomic number in the atomic data
/// a negative value is returned to flag this event.
{
  for (int i = 0; i < get_nelements(); i++) {
    // printf("i %d, Z %d, elements[i].anumber %d\n",i,Z,elements[i].anumber);
    if (Z == globals::elements[i].anumber) {
      return i;
    }
  }

  // printout("[debug] get_elementindex: element Z=%d was not found in atomic data ... skip readin of cross sections for
  // this element\n",Z); printout("[fatal] get_elementindex: element Z=%d was not found in atomic data ... abort\n");
  // abort();;
  return -100;
}

void increase_includedions(const int nions) { includedions += nions; }

auto get_includedions() -> int
// returns the number of ions of all elements combined
{
  return includedions;
}

void update_max_nions(const int nions)
// Will ensure that maxnions is always greater than or equal to the number of nions
// this is called at startup once per element with the number of ions
{
  if (nions > maxnions || maxnions < 0) {
    maxnions = nions;
  }
}

auto get_max_nions() -> int {
  // number greater than or equal to nions(element) for all elements
  return maxnions;
}

auto get_nions(const int element) -> int
/// Returns the number of ions associated with a specific element given by
/// its elementindex.
{
  assert_testmodeonly(element < get_nelements());
  return globals::elements[element].nions;
}

auto get_ionstage(const int element, const int ion) -> int
/// Returns the ionisation stage of an ion specified by its elementindex and
/// ionindex.
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].ionstage;
}

auto get_nlevels(const int element, const int ion) -> int
/// Returns the number of levels associated with with a specific ion given
/// its elementindex and ionindex.
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].nlevels;
}

auto get_nlevels_nlte(const int element, const int ion) -> int
// Returns the number of NLTE levels associated with with a specific ion given
// its elementindex and ionindex. Includes the superlevel if there is one but does not include the ground state
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].nlevels_nlte;
}

auto get_nlevels_groundterm(const int element, const int ion) -> int {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].nlevels_groundterm;
}

auto get_ionisinglevels(const int element, const int ion) -> int
/// Returns the number of levels associated with an ion that
/// have energies below the ionisation threshold.
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].ionisinglevels;
}

auto get_uniqueionindex(const int element, const int ion) -> int
// Get an index for an ionstage of an element that is unique for every ion of every element
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  int index = 0;
  for (int e = 0; e < element; e++) {
    index += get_nions(e);
  }
  index += ion;

  assert_testmodeonly(index == globals::elements[element].ions[ion].uniqueionindex);
  assert_testmodeonly(index < get_includedions());
  return index;
}

void get_ionfromuniqueionindex(const int allionsindex, int *element, int *ion) {
  assert_testmodeonly(allionsindex < get_includedions());
  int allionsindex_thiselementfirstion = 0;
  for (int e = 0; e < get_nelements(); e++) {
    if ((allionsindex - allionsindex_thiselementfirstion) >= get_nions(e)) {
      allionsindex_thiselementfirstion += get_nions(e);  // skip this element
    } else {
      *element = e;
      *ion = allionsindex - allionsindex_thiselementfirstion;
      assert_testmodeonly(get_uniqueionindex(*element, *ion) == allionsindex);
      return;
    }
  }
  assert_always(false);  // allionsindex too high to be valid
}

auto get_uniquelevelindex(const int element, const int ion, const int level) -> int
// Get an index for level of an ionstage of an element that is unique across every ion of every element
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));

  int index = 0;
  for (int e = 0; e < element; e++) {
    const int nions = get_nions(e);
    for (int i = 0; i < nions; i++) {
      index += get_nlevels(e, i);
    }
  }
  // selected element, levels from lower ions
  for (int i = 0; i < ion; i++) {
    index += get_nlevels(element, i);
  }
  // lower levels in selected element/ion
  index += level;

  assert_testmodeonly(index == globals::elements[element].ions[ion].levels[level].uniquelevelindex);
  return index;
}

void get_levelfromuniquelevelindex(const int alllevelsindex, int *element, int *ion, int *level)
// inverse of get_uniquelevelindex(). get the element/ion/level from a unique level index
{
  int allionsindex_thisionfirstlevel = 0;
  for (int e = 0; e < get_nelements(); e++) {
    const int nions = get_nions(e);
    for (int i = 0; i < nions; i++) {
      if ((alllevelsindex - allionsindex_thisionfirstlevel) >= get_nlevels(e, i)) {
        allionsindex_thisionfirstlevel += get_nlevels(e, i);  // skip this ion
      } else {
        *element = e;
        *ion = i;
        *level = alllevelsindex - allionsindex_thisionfirstlevel;
        assert_testmodeonly(get_uniquelevelindex(*element, *ion, *level) == alllevelsindex);
        return;
      }
    }
  }
  assert_always(false);  // alllevelsindex too high to be valid
}

auto epsilon(const int element, const int ion, const int level) -> double
/// Returns the energy of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return globals::elements[element].ions[ion].levels[level].epsilon;
}

auto stat_weight(const int element, const int ion, const int level) -> double
/// Returns the statistical weight of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return globals::elements[element].ions[ion].levels[level].stat_weight;
}

auto get_maxrecombininglevel(const int element, const int ion) -> int
/// Returns the number of bf-continua associated with ion ion of element element.
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return globals::elements[element].ions[ion].maxrecombininglevel;
}

auto ion_has_superlevel(const int element, const int ion) -> bool {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  return (get_nlevels(element, ion) > get_nlevels_nlte(element, ion) + 1);
}

auto get_ndowntrans(const int element, const int ion, const int level) -> int
// the number of downward bound-bound transitions from the specified level
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return globals::elements[element].ions[ion].levels[level].ndowntrans;
}

auto get_nuptrans(const int element, const int ion, const int level) -> int
// the number of upward bound-bound transitions from the specified level
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return globals::elements[element].ions[ion].levels[level].nuptrans;
}

void set_ndowntrans(const int element, const int ion, const int level, const int ndowntrans)
// the number of downward bound-bound transitions from the specified level
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  globals::elements[element].ions[ion].levels[level].ndowntrans = ndowntrans;
}

void set_nuptrans(const int element, const int ion, const int level, const int nuptrans)
// the number of upward bound-bound transitions from the specified level
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  globals::elements[element].ions[ion].levels[level].nuptrans = nuptrans;
}

auto get_nphixstargets(const int element, const int ion, const int level) -> int
/// Returns the number of target states for photoionization of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  const int nions = get_nions(element);
  const int nionisinglevels = get_ionisinglevels(element, ion);
  if ((ion < nions - 1) && (level < nionisinglevels)) {
    return globals::elements[element].ions[ion].levels[level].nphixstargets;
  }
  return 0;
}

auto get_phixsupperlevel(const int element, const int ion, const int level, const int phixstargetindex) -> int
/// Returns the level index of a target state for photoionization of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(element, ion, level));

  return globals::elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].levelindex;
}

auto get_phixs_threshold(const int element, const int ion, const int level, const int phixstargetindex) -> double
/// Returns the energy of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  assert_testmodeonly(phixstargetindex < get_nphixstargets(element, ion, level));
  // const double phixs_threshold_stored = globals::elements[element].ions[ion].levels[level].phixs_threshold;
  // if (phixs_threshold_stored > 0.)
  //   return phixs_threshold_stored;
  // else
  {
    const int upperlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
    const double E_threshold = epsilon(element, ion + 1, upperlevel) - epsilon(element, ion, level);
    return E_threshold;
  }
}

auto get_phixsprobability(const int element, const int ion, const int level, const int phixstargetindex) -> double
/// Returns the probability of a target state for photoionization of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(element, ion, level));

  return globals::elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].probability;
}

auto einstein_spontaneous_emission(const int lineindex) -> double
// double einstein_spontaneous_emission(int element, int ion, int upper, int lower)
/// reads A_ul from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... |
/// (A_upper,1; f_upper,1)
{
  return globals::linelist[lineindex].einstein_A;
}

auto photoionization_crosssection(const int element, const int ion, const int level, const double nu_edge,
                                  const double nu) -> double {
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return photoionization_crosssection_fromtable(globals::elements[element].ions[ion].levels[level].photoion_xs, nu_edge,
                                                nu);
}
