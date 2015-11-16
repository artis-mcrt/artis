#include "sn3d.h"


///****************************************************************************
int get_element(int element)
/// Returns the atomic number associated with a given elementindex.
{
  return elements[element].anumber;
}


///****************************************************************************
int get_elementindex(int Z)
/// Returns the elementindex associated with a given atomic number.
/// If there is no element with the given atomic number in the atomic data
/// a negative value is returned to flag this event.
{
  int i;
  
  for (i = 0; i < nelements; i++)
  {
    //printf("i %d, Z %d, elements[i].anumber %d\n",i,Z,elements[i].anumber);
    if (Z == elements[i].anumber) return i;
  }
  
  //printout("[debug] get_elementindex: element Z=%d was not found in atomic data ... skip readin of cross sections for this element\n",Z);
  //printout("[fatal] get_elementindex: element Z=%d was not found in atomic data ... abort\n");
  //exit(0);
  return -100;
}


///****************************************************************************
int get_nions(int element)
/// Returns the number of ions associated with a specific element given by
/// its elementindex.
{
  return elements[element].nions;
}



///****************************************************************************
int get_ionstage(int element, int ion)
/// Returns the ionisation stage of an ion specified by its elementindex and
/// ionindex.
{
  return elements[element].ions[ion].ionstage;
}


///****************************************************************************
int get_nlevels(int element, int ion)
/// Returns the number of levels associated with with a specific ion given 
/// ist elementindex and ionindex.
{
  return elements[element].ions[ion].nlevels; 
}

///****************************************************************************
int get_nlevels_nlte(int element, int ion)
/// Returns the number of levels associated with with a specific ion given 
/// ist elementindex and ionindex.
{
  return elements[element].ions[ion].nlevels_nlte; 
}


///****************************************************************************
int get_ionisinglevels(int element, int ion)
/// Returns the number of levels associated with ion ion of element element which 
/// have energies below the ionisation threshold.
{
  return elements[element].ions[ion].ionisinglevels; 
}


///****************************************************************************
int get_bfcontinua(int element, int ion)
/// Returns the number of bf-continua associated with ion ion of element element.
{
  int get_ionisinglevels(int element, int ion);
  int nionisinglevels;
  
  nionisinglevels = get_ionisinglevels(element,ion);
  
  if (nionisinglevels < max_bf_continua)
    return nionisinglevels;
  else
    return max_bf_continua; 
}


///***************************************************************************/
double epsilon(int element, int ion, int level)
/// Returns the energy of (element,ion,level).
{
  return elements[element].ions[ion].levels[level].epsilon;
}


///***************************************************************************/
double stat_weight(int element, int ion, int level)
/// Returns the statistical weight of (element,ion,level).
{
  return elements[element].ions[ion].levels[level].stat_weight;
}

///***************************************************************************/
short is_nlte(int element, int ion, int level)
/// Returns 1 if (element,ion,level) is to be treated in nlte.
{
  if (level < 201)
    {
      elements[element].ions[ion].levels[level].is_nlte = 1;
    }
  else
    {
      elements[element].ions[ion].levels[level].is_nlte = 0;
    }
  return elements[element].ions[ion].levels[level].is_nlte;
}



///***************************************************************************/
int get_continuumindex(int element, int ion, int level)
/// Returns the index of the continuum associated to the given level.
{
  return elements[element].ions[ion].levels[level].cont_index;
}

///***************************************************************************/
int get_nphixstargets(int element, int ion, int level)
/// Returns the number of target states for photoionization of (element,ion,level).
{
  int nions, nionisinglevels;
  
  nions = get_nions(element);
  nionisinglevels = get_ionisinglevels(element,ion);
  if ((ion < nions-1) && (level < nionisinglevels))
  {
    return elements[element].ions[ion].levels[level].nphixstargets;
  }
  else
  {
    return 0;
  }
}

///***************************************************************************/
int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex)
/// Returns the level index of a target state for photoionization of (element,ion,level).
{
  int nphixstargets;
  
  nphixstargets = get_nphixstargets(element,ion,level);
  if (phixstargetindex < nphixstargets)
  {
    return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].levelindex;
  }
  else
  {
    //should throw an error? or just silenty return the top phixstarget levelindex
    return elements[element].ions[ion].levels[level].phixstargets[nphixstargets-1].levelindex;
  }
}

///***************************************************************************/
float get_phixsprobability(int element, int ion, int level, int phixstargetindex)
/// Returns the probability of a target state for photoionization of (element,ion,level).
{
  int nphixstargets;
  
  nphixstargets = get_nphixstargets(element,ion,level);
  if (phixstargetindex < nphixstargets)
  {
    return elements[element].ions[ion].levels[level].phixstargets[phixstargetindex].probability;
  }
  else
  {
    //should throw an error? or just silenty return the top phixstarget probability
    return elements[element].ions[ion].levels[level].phixstargets[nphixstargets-1].probability;
  }
}

///***************************************************************************/
int transitioncheck(int upper, int lower)
/// reads A_ul from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
  int index;
  int flag;
  
  index = (upper-lower) - 1;
  flag = transitions[upper].to[index];

  return flag;
}


///***************************************************************************/
double einstein_spontaneous_emission(int lineindex)
//double einstein_spontaneous_emission(int element, int ion, int upper, int lower)
/// reads A_ul from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
  //int index;
  double A_ul;
  
/*  index = (upper-lower) - 1;
  A_ul = elements[element].ions[ion].levels[upper].transitions[index].einstein_A;*/
  A_ul = linelist[lineindex].einstein_A;

  return A_ul;
}


///***************************************************************************/
double osc_strength(int lineindex)
//double osc_strength(int element, int ion, int upper, int lower)
/// reads f_lu from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
  return linelist[lineindex].osc_strength;
}

///***************************************************************************/
double coll_str(int lineindex)
{
  return linelist[lineindex].coll_str;
}

///***************************************************************************/
double statw_up(int lineindex)
{
  return elements[linelist[lineindex].elementindex].ions[linelist[lineindex].ionindex].levels[linelist[lineindex].upperlevelindex].stat_weight;
}
///***************************************************************************/
double statw_down(int lineindex)
{
  return elements[linelist[lineindex].elementindex].ions[linelist[lineindex].ionindex].levels[linelist[lineindex].lowerlevelindex].stat_weight;
}


///***************************************************************************/
double osc_strength_old(int lineindex)
//double osc_strength(int element, int ion, int upper, int lower)
/// reads f_lu from levellist which consists of
/// (epsilon_upper; 0) | (g_upper; 0) | (A_upper,upper-1; f_upper,upper-1) | (A_uppper,upper-2; f_upper,upper-2) | ... | (A_upper,1; f_upper,1)
{
  //int index;
  double f_ul;
  
/*  index = (upper-lower) - 1;
  f_ul = elements[element].ions[ion].levels[upper].transitions[index].oscillator_strength;*/
  f_ul = linelist[lineindex].osc_strength;
  return f_ul;
}



///***************************************************************************/
/// Alll the following deals with photoionisation cross-sections
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
  
  

///***************************************************************************/
double photoionization_crosssection(double nu_edge, double nu)
/// Calculates the photoionisation cross-section at frequency nu out of the atomic data.
/// Input: - edge frequency nu_edge of the desired bf-continuum
///        - nu
///        - BE AWARE: the elements of the global structure variable mastate
///                    must fit to the bound state of the desired bf-continuum!!!
{
  double epsilon(int element, int ion, int level);
  double sigma_bf;
  
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;
  int i;
  
  int check = 0;
  
  if (nu < nu_edge) 
  {
    printout("[warning] photoionization_crosssection was called with nu=%g < nu_edge=%g\n",nu,nu_edge);
    printout("[warning]   element %d, ion %d, level %d, epsilon %g, ionpot %g\n",element,ion,level,epsilon(element,ion,level),elements[element].ions[ion].ionpot);
    printout("[warning]   element %d, ion+1 %d, level %d epsilon %g, ionpot %g\n",element,ion+1,0,epsilon(element,ion+1,0),elements[element].ions[ion].ionpot);
    check=1;
    //abort();
  }
  
  if (nu == nu_edge)
  {
    sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[0];
  }
  else if (nu < nu_edge*(1+0.1*NPHIXSPOINTS))
  {
    i = floor(nu/(0.1*nu_edge)) - 10;
    sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[i];
  }
  else
  {
    /// use a parameterization of sigma_bf by the Kramers formula 
    /// which anchor point should we take ??? the cross-section at the edge or at the highest grid point ???
    /// so far the highest grid point, otherwise the cross-section is not continuous
    sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[NPHIXSPOINTS-1] * pow(nu_edge*(1+0.1*NPHIXSPOINTS)/nu,3);
  }
  
  if (check == 1) printout("[warning]   photoionization_crosssection %g\n",sigma_bf);
  if (sigma_bf < 0)
  {
    printout("[warning] photoionization_crosssection returns negative cross-section %g\n",sigma_bf);
    printout("[warning]   nu=%g,  nu_edge=%g, xs@edge=%g, xs@maxfreq\n",nu,nu_edge,elements[element].ions[ion].levels[level].photoion_xs[0],elements[element].ions[ion].levels[level].photoion_xs[NPHIXSPOINTS-1]);
    printout("[warning]   element %d, ion %d, level %d, epsilon %g, ionpot %g\n",element,ion,level,epsilon(element,ion,level),elements[element].ions[ion].ionpot);
    printout("[warning]   element %d, ion+1 %d, level %d epsilon %g, ionpot %g\n",element,ion+1,0,epsilon(element,ion+1,0),elements[element].ions[ion].ionpot);
  }
  
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
 
  if (nu == nu_edge)
  {
    sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[0];
  }
  else if (nu < nu_edge*(1+0.1*NPHIXSPOINTS))
  {
    double nu_step = 0.1*nu_edge;
    int lowerindex = floor((nu-nu_edge)/nu_step);
    int upperindex = lowerindex + 1;
    double f_upper = elements[element].ions[ion].levels[level].photoion_xs[upperindex];
    double f_lower = elements[element].ions[ion].levels[level].photoion_xs[lowerindex];
    double nu_upper =  nu_edge + upperindex*nu_step;
    double nu_lower =  nu_edge + lowerindex*nu_step;
    sigma_bf = f_lower + (f_upper-f_lower)/(nu_upper-nu_lower) * (nu-nu_lower);
  }
  else
  {
    /// use a parameterization of sigma_bf by the Kramers formula 
    /// which anchor point should we take ??? the cross-section at the edge or at the highest grid point ???
    /// so far the highest grid point, otherwise the cross-section is not continuous
    sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[NPHIXSPOINTS-1] * pow(nu_edge*(1+0.1*NPHIXSPOINTS)/nu,3);
  }
  
  return sigma_bf;
}
*/