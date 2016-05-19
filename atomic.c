#include "sn3d.h"
#include "atomic.h"

///***************************************************************************/
/// All the following deals with photoionisation cross-sections
/// Old version with only Kramers-like cross-sections,
/// and different versions of readin from a lookup table
/// The current active version is split into subfunctions for profiling
/// issues.
/*
double photoionization_crosssection(double nu_edge, double nu)
/// Calculates the photoionisation cross-section at frequency nu out of the atomic data.
/// Input: - edge frequency nu_edge of the desired bf-continuum
///        - nu
///        - BE AWARE: the elements of the global structure variable mastate
///                    must fit to the bound state of the desired bf-continuum
{
  double sigma_bf;
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  ///use a parameterization of sigma_bf by the Kramers formula and the cross-section at the edge
  sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[0] * pow(nu_edge/nu,3);

  return sigma_bf;
}
*/



///***************************************************************************/
/*double photoionization_crosssection_new(double nu_edge, double nu)
/// Calculates the photoionisation cross-section at frequency nu out of the atomic data.
/// Input: - edge frequency nu_edge of the desired bf-continuum
///        - nu
///        - BE AWARE: the elements of the global structure variable mastate
///                    must fit to the bound state of the desired bf-continuum
/// The cross sections are stored temporarily, to save on interpolation.
{
  double sigma_bf_stored,sigma_bf;
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;by(
  int level = mastate[tid].level;

  double nu_left,nu_right,sigma_left,sigma_right;
  int match;

  int npoints = elements[element].ions[ion].levels[level].photoion_xs_ngridpoints;

  int left = 0;
  int right = npoints-1;
  int middle=1;
  //printout("nu %g\n",nu);
  //nu_edge = elements[element].ions[ion].levels[level].photoion_xs[0].nu;
  if (nu == nu_edge)
  {
    /// return the egde cross-section
    sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[0].xs;
    elements[element].ions[ion].levels[level].currentphotoion_xs = sigma_bf;
  }
  else
  {
    sigma_bf_stored = elements[element].ions[ion].levels[level].currentphotoion_xs;
    if (sigma_bf_stored > 0)
    {
      sigma_bf = sigma_bf_stored;
    }
    else if (nu <= elements[element].ions[ion].levels[level].photoion_xs[npoints-1].nu)
    {
    /// find position of current frequency nu relative to the available sigma_bf frequency grid
      while (left <= right)
      {
        middle = left + ((right-left) / 2);

        if (nu >= elements[element].ions[ion].levels[level].photoion_xs[middle].nu
            && nu < elements[element].ions[ion].levels[level].photoion_xs[middle+1].nu) break;

        if (nu >= elements[element].ions[ion].levels[level].photoion_xs[middle+1].nu)
          left = middle + 1;
        else
          right = middle - 1;
      }
      match = middle;
    /// interpolate the photoionisation cross-section for the current frequency
    /// out of the next  neighbours
      sigma_left = elements[element].ions[ion].levels[level].photoion_xs[middle].xs;
      sigma_right = elements[element].ions[ion].levels[level].photoion_xs[middle+1].xs;
      nu_left = elements[element].ions[ion].levels[level].photoion_xs[middle].nu;
      nu_right = elements[element].ions[ion].levels[level].photoion_xs[middle+1].nu;
    //printout("nu - nu_left %g, nu %g, nu_left %g, nu_right %g\n",nu-nu_left,nu,nu_left,nu_right);
    //printout("nu_left, nu_right %g %g\n",nu_left,nu_right);

      sigma_bf = sigma_left + (sigma_right-sigma_left)/(nu_right-nu_left) * (nu-nu_left);
      elements[element].ions[ion].levels[level].currentphotoion_xs = sigma_bf;
    //printout("sigma_bf %g\n",sigma_bf);
    }
    else
    {
    /// use a parameterization of sigma_bf by the Kramers formula
    /// which anchor point should we take ??? the cross-section at the edge or at the highest grid point ???
    /// so far the highest grid point, otherwise the cross-section is not continuous
      sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[npoints-1].xs * pow(elements[element].ions[ion].levels[level].photoion_xs[npoints-1].nu/nu,3);
      elements[element].ions[ion].levels[level].currentphotoion_xs = sigma_bf;
    }
  }

  return sigma_bf;
}
*/


///should this be inlined or not?
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