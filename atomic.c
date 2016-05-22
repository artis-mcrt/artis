#include "sn3d.h"
#include "atomic.h"


inline
double photoionization_crosssection(double nu_edge, double nu)
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

  if (nu == nu_edge)
  {
    sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[0];
  }
  else
  {
    int i = floor((nu/nu_edge - 1.0)/NPHIXSNUINCREMENT);

    #ifdef DEBUG_ON
    if (i < 0)
    {
      sigma_bf = 0.0;
      printout("[warning] photoionization_crosssection was called with nu=%g < nu_edge=%g\n",nu,nu_edge);
      printout("[warning]   element %d, ion %d, level %d, epsilon %g, ionpot %g\n",element,ion,level,epsilon(element,ion,level),elements[element].ions[ion].ionpot);
      printout("[warning]   element %d, ion+1 %d, level %d epsilon %g, ionpot %g\n",element,ion+1,0,epsilon(element,ion+1,0),elements[element].ions[ion].ionpot);
      printout("[warning]   photoionization_crosssection %g\n",sigma_bf);
      //abort();
    }
    else if (i < NPHIXSPOINTS)
    #else
    if (i < NPHIXSPOINTS)
    #endif
    {
      sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[i];
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

  #ifdef DEBUG_ON
    if (sigma_bf < 0)
    {
      printout("[warning] photoionization_crosssection returns negative cross-section %g\n",sigma_bf);
      printout("[warning]   nu=%g,  nu_edge=%g, xs@edge=%g, xs@maxfreq\n",nu,nu_edge,elements[element].ions[ion].levels[level].photoion_xs[0],elements[element].ions[ion].levels[level].photoion_xs[NPHIXSPOINTS-1]);
      printout("[warning]   element %d, ion %d, level %d, epsilon %g, ionpot %g\n",element,ion,level,epsilon(element,ion,level),elements[element].ions[ion].ionpot);
    }
  #endif

  return sigma_bf;
}


///***************************************************************************/
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