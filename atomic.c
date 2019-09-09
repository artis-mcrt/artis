#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "ltepop.h"
#include "update_grid.h"

extern inline int get_element(int element);
extern inline int get_elementindex(int Z);
extern inline int get_nions(int element);
extern inline int get_ionstage(int element, int ion);
extern inline int get_nlevels(int element, int ion);
extern inline int get_nlevels_nlte(int element, int ion);
extern inline int get_nlevels_groundterm(int element, int ion);
extern inline int get_ionisinglevels(int element, int ion);
extern inline int get_uniqueionindex(int element, int ion);
extern inline void get_ionfromuniqueionindex(int allionsindex, int *element, int *ion);
extern inline double epsilon(int element, int ion, int level);
extern inline double stat_weight(int element, int ion, int level);
extern inline int get_maxrecombininglevel(int element, int ion);
extern inline bool is_nlte(int element, int ion, int level);
extern inline bool ion_has_superlevel(int element, int ion);
extern inline int get_ndowntrans(int element, int ion, int level);
extern inline int get_nuptrans(int element, int ion, int level);
extern inline int get_nphixstargets(int element, int ion, int level);
extern inline int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
extern inline double get_phixsprobability(int element, int ion, int level, int phixstargetindex);
extern inline void set_ndowntrans( int element,  int ion,  int level,  int ndowntrans);
extern inline void set_nuptrans( int element, int ion,  int level,  int nuptrans);
extern inline double einstein_spontaneous_emission(int lineindex);
extern inline double osc_strength(int lineindex);
extern inline double get_coll_str(int lineindex);
extern inline double statw_upper(int lineindex);
extern inline double statw_lower(int lineindex);
extern inline double photoionization_crosssection(int element, int ion, int level, double nu_edge, double nu);
extern inline double get_phixs_threshold(int element, int ion, int level, int phixstargetindex);


static long get_continuumindex_phixstargetindex(int element, int ion, int level, int phixstargetindex)
/// Returns the index of the continuum associated to the given level.
{
  return elements[element].ions[ion].levels[level].cont_index - phixstargetindex;
}


int get_phixtargetindex(const int element, const int ion, const int level, const int upperionlevel)
{
  for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++)
  {
    if (upperionlevel == get_phixsupperlevel(element, ion, level, phixstargetindex))
      return phixstargetindex;
  }
  printout("Could not find phixstargetindex\n");
  abort();
}


long get_continuumindex(int element, int ion, int level, int upperionlevel)
/// Returns the index of the continuum associated to the given level.
{
  const int phixstargetindex = get_phixtargetindex(element, ion, level, upperionlevel);
  return get_continuumindex_phixstargetindex(element, ion, level, phixstargetindex);
}


double get_tau_sobolev(int modelgridindex, int lineindex, double t_current)
{
  const int element = linelist[lineindex].elementindex;
  const int ion = linelist[lineindex].ionindex;
  const int lower = linelist[lineindex].lowerlevelindex;
  const int upper = linelist[lineindex].upperlevelindex;

  const double n_l = calculate_exclevelpop(modelgridindex,element,ion,lower);
  const double n_u = calculate_exclevelpop(modelgridindex,element,ion,upper);

  const double nu_trans = (epsilon(element, ion, upper) - epsilon(element, ion, lower)) / H;
  const double A_ul = einstein_spontaneous_emission(lineindex);
  const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
  const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

  const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;
  return tau_sobolev;
}


int get_tot_nions(void)
{
  int nions = 0.;
  for (int element = 0; element < nelements; element++)
  {
    nions += get_nions(element);
  }

  return nions;
}


double get_nntot(int modelgridindex)
// total ion (nuclei) density
{
  const double rho = get_rho(modelgridindex);
  int nntot = 0.;
  for (int element = 0; element < nelements; element++)
  {
    nntot += get_abundance(modelgridindex, element) / elements[element].mass * rho;
  }

  return nntot;
}

bool is_nlte(int element, int ion, int level)
// Returns true if (element,ion,level) is to be treated in nlte.
// (note this function returns true for the ground state,
//  although it is stored separately from the excited NLTE states)
{
  if (!NLTE_POPS_ON)
  {
    return false;
  }
  else if (get_element(element) == 26 && get_ionstage(element, ion) == 2)
    return (level <= 197);
  else
    return (level <= 80);
}


bool level_isinsuperlevel(int element, int ion, int level)
// ion has NLTE levels, but this one is not NLTE => is in the superlevel
{
  return (NLTE_POPS_ON && !is_nlte(element,ion,level) && level != 0 && (get_nlevels_nlte(element, ion) > 0));
}


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
  const double ireal = (nu / nu_edge - 1.0) / NPHIXSNUINCREMENT;
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
  else if (i < NPHIXSPOINTS - 1)
  {
    // sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[i];

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
    sigma_bf = photoion_xs[NPHIXSPOINTS-1] * pow(nu_max_phixs / nu, 3);
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
