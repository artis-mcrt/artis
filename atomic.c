#include "sn3d.h"
#include "atomic.h"
#include "ltepop.h"

extern inline int get_element(int element);
extern inline int get_elementindex(int Z);
extern inline int get_nions(int element);
extern inline int get_ionstage(int element, int ion);
extern inline int get_nlevels(int element, int ion);
extern inline int get_nlevels_nlte(int element, int ion);
extern inline int get_ionisinglevels(int element, int ion);
extern inline double epsilon(int element, int ion, int level);
extern inline double stat_weight(int element, int ion, int level);
extern inline int get_bfcontinua(int element, int ion);
extern inline bool is_nlte(int element, int ion, int level);
extern inline int get_continuumindex(int element, int ion, int level);
extern inline int get_nphixstargets(int element, int ion, int level);
extern inline int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
extern inline double get_phixsprobability(int element, int ion, int level, int phixstargetindex);
extern inline int transitioncheck(int upper, int lower);
extern inline double einstein_spontaneous_emission(int lineindex);
extern inline double osc_strength(int lineindex);
extern inline double get_coll_str(int lineindex);
extern inline double statw_upper(int lineindex);
extern inline double statw_lower(int lineindex);
extern inline double photoionization_crosssection(double nu_edge, double nu);
extern inline double xs_photoionization(int element, int ion, int level, double nu_edge, double nu);

double get_tau_sobolev(int modelgridindex, int lineindex, double t_current)
{
  const int element = linelist[lineindex].elementindex;
  const int ion = linelist[lineindex].ionindex;
  const int lower = linelist[lineindex].lowerlevelindex;
  const int upper = linelist[lineindex].upperlevelindex;

  const double statweight_target = statw_upper(lineindex);
  const double statweight_lower = statw_lower(lineindex);

  const double n_l = get_levelpop(modelgridindex,element,ion,lower);
  const double n_u = get_levelpop(modelgridindex,element,ion,upper);

  const double nu_trans = (epsilon(element, ion, upper) - epsilon(element, ion, lower)) / H;
  const double A_ul = einstein_spontaneous_emission(lineindex);
  const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
  const double B_lu = statweight_target / statweight_lower * B_ul;

  const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * t_current;
  return tau_sobolev;
}


/*double interpolate_photoionization_crosssection(double nu_edge, double nu)
/// Calculates the photoionisation cross-section at frequency nu out of the atomic data.
/// Input: - edge frequency nu_edge of the desired bf-continuum
///        - nu
///        - BE AWARE: the elements of the global structure variable mastate
///                    must fit to the bound state of the desired bf-continuum!!!
{
  double sigma_bf;

  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  lowerindex = floor((nu/nu_edge - 1.0)/NPHIXSNUINCREMENT);
  if (lowerindex < 0):
  {
    sigma_bf = 0.0;
    #ifdef DEBUG_ON
      printout("[warning] interpolate_photoionization_crosssection was called with nu=%g < nu_edge=%g\n",nu,nu_edge);
      printout("[warning]   element %d, ion %d, level %d, epsilon %g, ionpot %g\n",element,ion,level,epsilon(element,ion,level),elements[element].ions[ion].ionpot);
      printout("[warning]   element %d, ion+1 %d, level %d epsilon %g, ionpot %g\n",element,ion+1,0,epsilon(element,ion+1,0),elements[element].ions[ion].ionpot);
      printout("[warning]   photoionization_crosssection %g\n",sigma_bf);
      //abort();
    #endif
  }
  else if (lowerindex < NPHIXSPOINTS)
  {
    if (lowerindex+1 < NPHIXSPOINTS)
    {
      double f_upper = elements[element].ions[ion].levels[level].photoion_xs[lowerindex+1];
      double f_lower = elements[element].ions[ion].levels[level].photoion_xs[lowerindex];
      double nu_upper = nu_edge * (1.0 + NPHIXSNUINCREMENT*(lowerindex+1));
      double nu_lower = nu_edge * (1.0 + NPHIXSNUINCREMENT*(lowerindex));
      sigma_bf = f_lower + (f_upper-f_lower) / (nu_upper-nu_lower) * (nu-nu_lower);
    }
    else
    {
      /// use a parameterization of sigma_bf by the Kramers formula
      /// which anchor point should we take ??? the cross-section at the edge or at the highest grid point ???
      /// so far the highest grid point, otherwise the cross-section is not continuous
      double nu_max_phixs = nu_edge * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table
      sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[NPHIXSPOINTS-1] * pow(nu_max_phixs/nu, 3);
    }
  }

  return sigma_bf;
}
*/