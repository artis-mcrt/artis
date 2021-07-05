#include "artisoptions.h"
#include "sn3d.h"
#include "atomic.h"
#include "grid.h"
#include "ltepop.h"
#include "update_grid.h"

__managed__ double last_phixs_nuovernuedge; // last photoion cross section point as a factor of nu_edge = last_phixs_nuovernuedge
__managed__ int nelements = 0;  // total number of elements included in the simulation
__managed__ int maxnions = 0;  // highest number of ions for any element

__host__ __device__
static int get_continuumindex_phixstargetindex(int element, int ion, int level, int phixstargetindex)
/// Returns the index of the continuum associated to the given level.
{
  return globals::elements[element].ions[ion].levels[level].cont_index - phixstargetindex;
}


__host__ __device__
int get_phixtargetindex(const int element, const int ion, const int level, const int upperionlevel)
{
  for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
  {
    if (upperionlevel == get_phixsupperlevel(element, ion, level, phixstargetindex))
      return phixstargetindex;
  }
  printout("Could not find phixstargetindex\n");
  abort();
  return -1;
}


__host__ __device__
int get_continuumindex(int element, int ion, int level, int upperionlevel)
/// Returns the index of the continuum associated to the given level.
{
  const int phixstargetindex = get_phixtargetindex(element, ion, level, upperionlevel);
  return get_continuumindex_phixstargetindex(element, ion, level, phixstargetindex);
}


__host__ __device__
double get_tau_sobolev(int modelgridindex, int lineindex, double t_current)
{
  const int element = globals::linelist[lineindex].elementindex;
  const int ion = globals::linelist[lineindex].ionindex;
  const int lower = globals::linelist[lineindex].lowerlevelindex;
  const int upper = globals::linelist[lineindex].upperlevelindex;

  const double n_l = calculate_exclevelpop(modelgridindex,element,ion,lower);
  const double n_u = calculate_exclevelpop(modelgridindex,element,ion,upper);

  const double nu_trans = (epsilon(element, ion, upper) - epsilon(element, ion, lower)) / H;
  const double A_ul = einstein_spontaneous_emission(lineindex);
  const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
  const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

  const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;
  return tau_sobolev;
}


__host__ __device__
double get_nntot(int modelgridindex)
// total ion (nuclei) density
{
  const double rho = grid::get_rho(modelgridindex);
  int nntot = 0.;
  for (int element = 0; element < get_nelements(); element++)
  {
    nntot += grid::get_elem_abundance(modelgridindex, element) / globals::elements[element].mass * rho;
  }

  return nntot;
}


__host__ __device__
bool is_nlte(int element, int ion, int level)
// Returns true if (element,ion,level) is to be treated in nlte.
// (note this function returns true for the ground state,
//  although it is stored separately from the excited NLTE states)
{
  if (!NLTE_POPS_ON)
  {
    return false;
  }
  else
  {
    LEVEL_IS_NLTE(element, ion, level);  // macro function defined in artisoptions.h
  }
}


__host__ __device__
bool level_isinsuperlevel(int element, int ion, int level)
// ion has NLTE levels, but this one is not NLTE => is in the superlevel
{
  return (NLTE_POPS_ON && !is_nlte(element,ion,level) && level != 0 && (get_nlevels_nlte(element, ion) > 0));
}


__host__ __device__
double photoionization_crosssection_fromtable(float *photoion_xs, double nu_edge, double nu)
/// Calculates the photoionisation cross-section at frequency nu out of the atomic data.
/// Input: - edge frequency nu_edge of the desired bf-continuum
///        - nu
{
  // if (nu < nu_edge || nu > nu_edge * 1.05)
  //   return 0;
  // else
  //   return 1.;
  // return 1. * pow(nu_edge / nu, 3);

  float sigma_bf;
  const double ireal = (nu / nu_edge - 1.0) / globals::NPHIXSNUINCREMENT;
  const int i = floor(ireal);

  if (i < 0)
  {
    sigma_bf = 0.0;
    //printout("[warning] photoionization_crosssection was called with nu=%g < nu_edge=%g\n",nu,nu_edge);
    //printout("[warning]   element %d, ion %d, level %d, epsilon %g, ionpot %g\n",element,ion,level,epsilon(element,ion,level),elements[element].ions[ion].ionpot);
    //printout("[warning]   element %d, ion+1 %d, level %d epsilon %g, ionpot %g\n",element,ion+1,0,epsilon(element,ion+1,0),elements[element].ions[ion].ionpot);
    //printout("[warning]   photoionization_crosssection %g\n",sigma_bf);
    //abort();
  }
  else if (i < globals::NPHIXSPOINTS - 1)
  {
    // sigma_bf = globals::elements[element].ions[ion].levels[level].photoion_xs[i];

    const double sigma_bf_a = photoion_xs[i];
    const double sigma_bf_b = photoion_xs[i + 1];
    const double factor_b = ireal - i;
    sigma_bf = ((1. - factor_b) * sigma_bf_a) + (factor_b * sigma_bf_b);
  }
  else
  {
    /// use a parameterization of sigma_bf by the Kramers formula
    /// which anchor point should we take ??? the cross-section at the edge or at the highest grid point ???
    /// so far the highest grid point, otherwise the cross-section is not continuous
    const double nu_max_phixs = nu_edge * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table
    sigma_bf = photoion_xs[globals::NPHIXSPOINTS-1] * pow(nu_max_phixs / nu, 3);
  }

#ifdef DEBUG_ON
  // if (sigma_bf < 0)
  // {
  //   printout("[warning] photoionization_crosssection returns negative cross-section %g\n",sigma_bf);
  //   printout("[warning]   nu=%g,  nu_edge=%g\n",nu,nu_edge);
  //   printout("[warning]   xs@edge=%g, xs@maxfreq\n",elements[element].ions[ion].levels[level].photoion_xs[0],elements[element].ions[ion].levels[level].photoion_xs[NPHIXSPOINTS-1]);
  //   printout("[warning]   element %d, ion %d, level %d, epsilon %g, ionpot %g\n",element,ion,level,epsilon(element,ion,level),elements[element].ions[ion].ionpot);
  // }
#endif

  return sigma_bf;
}


__host__ __device__
void set_nelements(const int nelements_in)
{
  nelements = nelements_in;
}


__host__ __device__
int get_nelements(void)
{
  return nelements;
}


__host__ __device__
int get_element(int element)
/// Returns the atomic number associated with a given elementindex.
{
  assert_testmodeonly(element >= 0);
  assert_testmodeonly(element < nelements);
  return globals::elements[element].anumber;
}


__host__ __device__
int get_elementindex(int Z)
/// Returns the elementindex associated with a given atomic number.
/// If there is no element with the given atomic number in the atomic data
/// a negative value is returned to flag this event.
{
  for (int i = 0; i < get_nelements(); i++)
  {
    //printf("i %d, Z %d, elements[i].anumber %d\n",i,Z,elements[i].anumber);
    if (Z == globals::elements[i].anumber)
      return i;
  }

  //printout("[debug] get_elementindex: element Z=%d was not found in atomic data ... skip readin of cross sections for this element\n",Z);
  //printout("[fatal] get_elementindex: element Z=%d was not found in atomic data ... abort\n");
  //abort();;
  return -100;
}


__host__ __device__
void update_max_nions(const int nions)
// Will ensure that maxnions is always greater than or equal to the number of nions
// this is called at startup once per element with the number of ions
{
  if (nions > maxnions || maxnions < 0)
  {
    maxnions = nions;
  }
}


__host__ __device__
int get_max_nions(void)
{
  // number greater than or equal to nions(element) for all elements
  return maxnions;
}


__host__ __device__
int get_nions(int element)
/// Returns the number of ions associated with a specific element given by
/// its elementindex.
{
  return globals::elements[element].nions;
}


__host__ __device__
int get_ionstage(int element, int ion)
/// Returns the ionisation stage of an ion specified by its elementindex and
/// ionindex.
{
  return globals::elements[element].ions[ion].ionstage;
}


__host__ __device__
int get_nlevels(int element, int ion)
/// Returns the number of levels associated with with a specific ion given
/// its elementindex and ionindex.
{
  return globals::elements[element].ions[ion].nlevels;
}


__host__ __device__
int get_nlevels_nlte(int element, int ion)
// Returns the number of NLTE levels associated with with a specific ion given
// its elementindex and ionindex. Includes the superlevel if there is one but does not include the ground state
{
  return globals::elements[element].ions[ion].nlevels_nlte;
}


__host__ __device__
int get_nlevels_groundterm(int element, int ion)
{
  return globals::elements[element].ions[ion].nlevels_groundterm;
}


__host__ __device__
int get_ionisinglevels(int element, int ion)
/// Returns the number of levels associated with an ion that
/// have energies below the ionisation threshold.
{
  return globals::elements[element].ions[ion].ionisinglevels;
}


__host__ __device__
int get_uniqueionindex(const int element, const int ion)
// Get an index for an ionstage of an element that is unique for every ion of every element
{
  int index = 0;
  for (int e = 0; e < element; e++)
  {
    index += get_nions(e);
  }
  index += ion;

  assert_testmodeonly(index == globals::elements[element].ions[ion].uniqueionindex);
  assert_testmodeonly(index < globals::includedions);
  return index;
}


__host__ __device__
void get_ionfromuniqueionindex(const int allionsindex, int *element, int *ion)
{
  int allionsindex_thiselementfirstion = 0;
  for (int e = 0; e < get_nelements(); e++)
  {
    if ((allionsindex - allionsindex_thiselementfirstion) >= get_nions(e))
    {
      allionsindex_thiselementfirstion += get_nions(e); // skip this element
    }
    else
    {
      *element = e;
      *ion = allionsindex - allionsindex_thiselementfirstion;
      assert_testmodeonly(get_uniqueionindex(*element, *ion) == allionsindex);
      return;
    }
  }
  *element = -1;
  *ion = -1;
}


__host__ __device__
double epsilon(int element, int ion, int level)
/// Returns the energy of (element,ion,level).
{
  assert_testmodeonly(element < get_nelements());
  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(level < get_nlevels(element, ion));
  return globals::elements[element].ions[ion].levels[level].epsilon;
}


__host__ __device__
double stat_weight(int element, int ion, int level)
/// Returns the statistical weight of (element,ion,level).
{
  assert_testmodeonly(level <= globals::elements[element].ions[ion].nlevels);
  return globals::elements[element].ions[ion].levels[level].stat_weight;
}


__host__ __device__
int get_maxrecombininglevel(int element, int ion)
/// Returns the number of bf-continua associated with ion ion of element element.
{
  return globals::elements[element].ions[ion].maxrecombininglevel;
}


__host__ __device__
bool ion_has_superlevel(const int element, const int ion)
{
  return (get_nlevels(element, ion) > get_nlevels_nlte(element, ion) + 1);
}


__host__ __device__
int get_ndowntrans(int element, int ion, int level)
// the number of downward bound-bound transitions from the specified level
{
  return globals::elements[element].ions[ion].levels[level].ndowntrans;
}


__host__ __device__
int get_nuptrans(int element, int ion, int level)
// the number of upward bound-bound transitions from the specified level
{
  return globals::elements[element].ions[ion].levels[level].nuptrans;
}


__host__ __device__
void set_ndowntrans(const int element, const int ion, const int level, const int ndowntrans)
// the number of downward bound-bound transitions from the specified level
{
  globals::elements[element].ions[ion].levels[level].ndowntrans = ndowntrans;
}


__host__ __device__
void set_nuptrans(const int element, const int ion, const int level, const int nuptrans)
// the number of upward bound-bound transitions from the specified level
{
  globals::elements[element].ions[ion].levels[level].nuptrans = nuptrans;
}


__host__ __device__
int get_nphixstargets(const int element, const int ion, const int level)
/// Returns the number of target states for photoionization of (element,ion,level).
{
  const int nions = get_nions(element);
  const int nionisinglevels = get_ionisinglevels(element,ion);
  if ((ion < nions-1) && (level < nionisinglevels))
    return globals::elements[element].ions[ion].levels[level].nphixstargets;
  else
  {
    return 0;
  }
}


__host__ __device__
int get_phixsupperlevel(const int element, const int ion, const int level, const int phixstargetindex)
/// Returns the level index of a target state for photoionization of (element,ion,level).
{
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(element,ion,level));

  return globals::elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].levelindex;
}


__host__ __device__
double get_phixs_threshold(int element, int ion, int level, int phixstargetindex)
/// Returns the energy of (element,ion,level).
{
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


__host__ __device__
double get_phixsprobability(int element, int ion, int level, int phixstargetindex)
/// Returns the probability of a target state for photoionization of (element,ion,level).
{
  assert_testmodeonly(phixstargetindex >= 0);
  assert_testmodeonly(phixstargetindex < get_nphixstargets(element,ion,level));

  return globals::elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].probability;
}


__host__ __device__
double einstein_spontaneous_emission(int lineindex)
//double einstein_spontaneous_emission(int element, int ion, int upper, int lower)
/// reads A_ul from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
/*
  int index = (upper-lower) - 1;
  double A_ul = globals::elements[element].ions[ion].levels[upper].transitions[index].einstein_A;
*/
  return globals::linelist[lineindex].einstein_A;
}


__host__ __device__
double osc_strength(int lineindex)
//double osc_strength(int element, int ion, int upper, int lower)
/// reads f_lu from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
  return globals::linelist[lineindex].osc_strength;
}


__host__ __device__
double get_coll_str(int lineindex)
{
  return globals::linelist[lineindex].coll_str;
}


__host__ __device__
double statw_upper(int lineindex)
{
  const int element = globals::linelist[lineindex].elementindex;
  const int ion = globals::linelist[lineindex].ionindex;
  const int upper = globals::linelist[lineindex].upperlevelindex;
  return globals::elements[element].ions[ion].levels[upper].stat_weight;
}


__host__ __device__
double statw_lower(int lineindex)
{
  const int element = globals::linelist[lineindex].elementindex;
  const int ion = globals::linelist[lineindex].ionindex;
  const int lower = globals::linelist[lineindex].lowerlevelindex;
  return globals::elements[element].ions[ion].levels[lower].stat_weight;
}


__host__ __device__
double photoionization_crosssection(int element, int ion, int level, double nu_edge, double nu)
{
  return photoionization_crosssection_fromtable(globals::elements[element].ions[ion].levels[level].photoion_xs, nu_edge, nu);
}


/*static double osc_strength_old(int lineindex)
//double osc_strength(int element, int ion, int upper, int lower)
/// reads f_lu from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{

  int index = (upper-lower) - 1;
  double f_ul = globals::elements[element].ions[ion].levels[upper].transitions[index].oscillator_strength;

  return f_ul;
}*/
