#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "nuclear.h"


double arr_nucdecayenergygamma[RADIONUCLIDE_COUNT] = {0};


double nucdecayenergygamma(enum radionuclides nuclide_type)
// average energy (erg) per decay in the form of gamma rays
{
  return arr_nucdecayenergygamma[nuclide_type];
}


void set_nucdecayenergygamma(enum radionuclides nuclide_type, double value)
{
  arr_nucdecayenergygamma[nuclide_type] = value;
}


double nucdecayenergypositrons(enum radionuclides nuclide_type)
// average energy (erg) per decay in the form of positrons
{
  switch (nuclide_type)
  {
    case NUCLIDE_NI57:
      return 0.354 * MEV * 0.436;
    case NUCLIDE_CO56:
      return 0.63 * MEV * 0.19;
    case NUCLIDE_V48:
      return 0.290 * 0.499 * MEV;
    default:
      return 0.;
  }
}


double nucdecayenergy(enum radionuclides nuclide_type)
{
  return nucdecayenergygamma(nuclide_type) + nucdecayenergypositrons(nuclide_type);
}


double meanlife(enum radionuclides nuclide_type)
{
  switch (nuclide_type)
  {
    case NUCLIDE_NI56:
      return 8.80 * DAY;
    case NUCLIDE_NI57:
      return 51.36 * 60;
    case NUCLIDE_CO56:
      return 113.7 * DAY;
    case NUCLIDE_CR48:
      return 1.29602 * DAY;
    case NUCLIDE_V48:
      return 23.0442 * DAY;
    case NUCLIDE_CO57:
      return 392.03 * DAY;
    case NUCLIDE_FE52:
      return 0.497429 * DAY;
    case NUCLIDE_MN52:
      return 0.0211395 * DAY;
    case FAKE_GAM_LINE_ID:
    case RADIONUCLIDE_COUNT:
      assert(false);
  }
}


double nucmass(enum radionuclides nuclide_type)
{
  switch (nuclide_type)
  {
    case NUCLIDE_NI56:
      return 56 * MH;
    case NUCLIDE_NI57:
      return 57 * MH;
    case NUCLIDE_CO56:
      return 56 * MH;
    case NUCLIDE_CR48:
      return 48 * MH;
    case NUCLIDE_V48:
      return 48 * MH;
    case NUCLIDE_CO57:
      return 57 * MH;
    case NUCLIDE_FE52:
      return 52 * MH;
    case NUCLIDE_MN52:
      return 52 * MH;
    case FAKE_GAM_LINE_ID:
    case RADIONUCLIDE_COUNT:
      assert(false);
  }
}


enum radionuclides decaynuc1(enum decaypathways selected_chain)
{
  switch (selected_chain)
  {
    case DECAY_NI56:
      return NUCLIDE_NI56;

    case DECAY_NI56_CO56:
      return NUCLIDE_NI56;

    case DECAY_FE52:
      return NUCLIDE_FE52;

    case DECAY_FE52_MN52:
      return NUCLIDE_FE52;

    case DECAY_CR48:
      return NUCLIDE_CR48;

    case DECAY_CR48_V48:
      return NUCLIDE_CR48;

    case DECAY_CO56:
      return NUCLIDE_CO56;

    case DECAY_NI57:
      return NUCLIDE_NI57;

    case DECAY_NI57_CO57:
      return NUCLIDE_NI57;

    case DECAY_CO57:
      return NUCLIDE_CO57;

    case DECAYPATH_COUNT:
      assert(false);
  }
}


enum radionuclides decaynuc2(enum decaypathways selected_chain)
{
  switch (selected_chain)
  {
    case DECAY_NI56:
    case DECAY_FE52:
    case DECAY_CR48:
    case DECAY_CO56:
    case DECAY_NI57:
    case DECAY_CO57:
    case DECAYPATH_COUNT:
      // value for no second radioactive nuclide (decay only from initial abund and not from radioactive parent)
      return RADIONUCLIDE_COUNT;

    case DECAY_NI56_CO56:
      return NUCLIDE_CO56;

    case DECAY_FE52_MN52:
      return NUCLIDE_MN52;

    case DECAY_CR48_V48:
      return NUCLIDE_V48;

    case DECAY_NI57_CO57:
      return NUCLIDE_CO57;
  }
}


static bool decaypath_is_chain(enum decaypathways selected_chain)
// the second nucleus is radioactive
{
  switch (selected_chain)
  {
    case DECAY_NI56:
    case DECAY_FE52:
    case DECAY_CR48:
    case DECAY_CO56:
    case DECAY_NI57:
    case DECAY_CO57:
      return false;

    case DECAY_NI56_CO56:
    case DECAY_FE52_MN52:
    case DECAY_CR48_V48:
    case DECAY_NI57_CO57:
      return true;

    case DECAYPATH_COUNT:
      assert(false)
  }
}


static enum radionuclides find_nucparent(enum radionuclides nuclide)
// get the parent nuclide, or if it doesn't have one then the value RADIONUCLIDE_COUNT is used
{
  for (enum decaypathways d = 0; d < DECAYPATH_COUNT; d++)
  {
    if (decaypath_is_chain(d))
    {
      enum radionuclides nuc2 = decaynuc2(d);
      if (nuc2 == nuclide)
      {
        return decaynuc1(d);
      }
    }
  }
  return RADIONUCLIDE_COUNT;
}


static void calculate_double_decay_chain(
  const double initabund1, const double meanlife1,
  const double initabund2, const double meanlife2,
  const double t_current,
  double *abund1, double *abund2, double *abund3)
{
  // calculate abundances from double decay, e.g., Ni56 -> Co56 -> Fe56
  // initabund3 is assumed to be zero, so the abundance of species 3 is only from decays of species 2

  const double lambda1 = 1 / meanlife1;
  const double lambda2 = 1 / meanlife2;

  const double newabund1 = initabund1 * exp(-lambda1 * t_current);

  const double newabund2 = (
    initabund2 * exp(-lambda2 * t_current) +
    initabund1 * lambda1 / (lambda1 - lambda2) * (exp(-lambda2 * t_current) - exp(-lambda1 * t_current)));

  const double newabund3 = (
    (initabund2 + initabund1) * (lambda1 - lambda2) -
    initabund2 * lambda1 * exp(-lambda2 * t_current) +
    initabund2 * lambda2 * exp(-lambda2 * t_current) -
    initabund1 * lambda1 * exp(-lambda2 * t_current) +
    initabund1 * lambda2 * exp(-lambda1 * t_current)) / (lambda1 - lambda2);

  // printout("calculate_double_decay_chain: initabund1 %g, initabund2 %g\n", initabund1, initabund2);

  // update the pointed to values after calculating all threedimensional
  // to ensure no aliasing problems, e.g. if abund1 = &initabund1
  *abund1 = newabund1;
  *abund2 = newabund2;
  *abund3 = newabund3;

  // printout("calculate_double_decay_chain: abund1 %g, abund2 %g abund3 %g\n", abund1, abund2, abund3);

  // ensure that the decays haven't altered the total abundance of all three species
  assert(fabs((initabund1 + initabund2) - (*abund1 + *abund2 + *abund3)) < 0.001);
}


static void calculate_doubledecay_modelabund(
  const int modelgridindex,
  enum radionuclides nuclide1,
  enum radionuclides nuclide2,
  const double t_current,
  double *abund1, double *abund2, double *abund3)
{
  const double initabund1 = get_modelinitradioabund(modelgridindex, nuclide1);
  const double meanlife1 = meanlife(nuclide1);
  const double initabund2 = get_modelinitradioabund(modelgridindex, nuclide2);
  const double meanlife2 = meanlife(nuclide2);

  calculate_double_decay_chain(initabund1, meanlife1, initabund2, meanlife2, t_current, abund1, abund2, abund3);
}


static double get_modelradioabund_at_time(
  const int modelgridindex, const enum radionuclides nuclide_type, const double time)
// get the mass fraction of a nuclide accounting for all decays including those of its parent
// e.g. Co56 abundance may first increase with time due to Ni56 decays, then decease due to Co56 decay
{
  if (time == 0)
  {
    return get_modelinitradioabund(modelgridindex, nuclide_type);
  }
  assert(nuclide_type != FAKE_GAM_LINE_ID);
  assert(nuclide_type < RADIONUCLIDE_COUNT);

  enum radionuclides nucparent = find_nucparent(nuclide_type);
  if (nucparent == RADIONUCLIDE_COUNT) // no parent exists, so use simple decay formula
  {
    return get_modelinitradioabund(modelgridindex, nuclide_type) * exp(- time / meanlife(nuclide_type));
  }
  else
  {
    // nuclide is part of a double-decay chain, e.g., Co56 in the chain: Ni56 -> Co56 -> Fe56
    double abund1 = 0.;
    double abund2 = 0.;
    double abund3 = 0.;
    calculate_doubledecay_modelabund(modelgridindex, nucparent, nuclide_type, time, &abund1, &abund2, &abund3);
    return abund2;
  }
}


static double get_modelinitradioabund_decayed(
  const int modelgridindex, const enum radionuclides nuclide_type, const double time)
// only allow decays to decrease the nuclide abundance (i.e. don't count increases due to decay of parent)
{
  return get_modelinitradioabund(modelgridindex, nuclide_type) * exp(- time / meanlife(nuclide_type));
}


static double get_endecay_per_ejectamass_at_time(const int mgi, enum decaypathways decaypath, double time)
// decay energy that would be released from time tstart to time infinity from each decaypath
{
  const enum radionuclides nuc1 = decaynuc1(decaypath);
  if (decaypath_is_chain(decaypath))
  {
    // double decay, e.g. DECAY_NI56_CO56, the decay of Co56 nuclei that were produced from Ni56 decays
    // decays from the second nuclide (e.g. Co56) due to the initial abundance are not counted here

    const enum radionuclides nuc2 = decaynuc2(decaypath);
    assert(nuc2 < RADIONUCLIDE_COUNT);
    const double initabund1 = get_modelinitradioabund(mgi, nuc1);
    const double initabund2 = 0.; // don't count initial abundance
    double abund1;
    double abund2;
    double abund3;
    calculate_double_decay_chain(initabund1, meanlife(nuc1), initabund2, meanlife(nuc2), time, &abund1, &abund2, &abund3);
    return (abund1 + abund2) / nucmass(nuc1) * nucdecayenergy(nuc2);
  }
  else
  {
    // simple decay from initial abundance , e.g. DECAY_NI56 or DECAY_CO56
    return get_modelinitradioabund_decayed(mgi, nuc1, time) / nucmass(nuc1) * nucdecayenergy(nuc1);
  }
}


static double get_endecay_per_ejectamass_between_times(
  const int mgi, enum decaypathways decaypath, double tlow, double thigh)
// energy per mass released by a decaypath between two times (in seconds)
{
  assert(tlow <= thigh);
  const double energy_tlow = get_endecay_per_ejectamass_at_time(mgi, decaypath, tlow);
  const double energy_thigh = get_endecay_per_ejectamass_at_time(mgi, decaypath, thigh);
  assert(energy_tlow >= energy_thigh);
  return energy_tlow - energy_thigh;
}


double get_simtime_endecay_per_ejectamass(const int mgi, enum decaypathways decaypath)
{
#ifdef NO_INITIAL_PACKETS
  return get_endecay_per_ejectamass_between_times(mgi, decaypath, tmin, tmax);
#else
  // allow decays from time zero
  return get_endecay_per_ejectamass_between_times(mgi, decaypath, 0., tmax);
#endif
}


double get_positroninjection_rate_density(const int modelgridindex, const double t)
// in erg / s / cm^3
{
  const double rho = get_rho(modelgridindex);

  const double co56frac = get_modelradioabund_at_time(modelgridindex, NUCLIDE_CO56, t);

  const double co56_positron_dep = nucdecayenergypositrons(NUCLIDE_CO56) * co56frac / meanlife(NUCLIDE_CO56) / nucmass(NUCLIDE_CO56) * rho;

  const double ni57frac = get_modelradioabund_at_time(modelgridindex, NUCLIDE_NI57, t);
  const double ni57_positron_dep = nucdecayenergypositrons(NUCLIDE_NI57) * ni57frac / meanlife(NUCLIDE_NI57) / nucmass(NUCLIDE_NI57) * rho;

  const double v48frac = get_modelradioabund_at_time(modelgridindex, NUCLIDE_V48, t);
  const double v48_positron_dep = nucdecayenergypositrons(NUCLIDE_V48)  * v48frac * rho;

  const double pos_dep_sum = co56_positron_dep + ni57_positron_dep + v48_positron_dep;
  printout("positroninjection_rate_density(mgi %d time %g): %g erg/s/cm3 = co56 %g + ni57 %g + v48 %g\n",
          modelgridindex, t, pos_dep_sum, co56_positron_dep, ni57_positron_dep, v48_positron_dep);

  return pos_dep_sum;
}


void update_abundances(const int modelgridindex, const int timestep, const double t_current)
/// Updates the mass fractions of elements associated with the decay sequence
/// Parameters: - modelgridindex: the grid cell for which to update the abundances
///             - t_current: current time (here mid of current timestep)
{
  assert(!homogeneous_abundances); // no longer supported

  const double timediff = t_current; // subtract t_model?

  // Ni56 -> Co56 -> Fe56
  // abundances from the input model
  double ni56frac = 0.;
  double co56frac = 0.;
  double fe56frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_NI56, NUCLIDE_CO56, timediff, &ni56frac, &co56frac, &fe56frac_fromdecay);

  // Ni57 -> Co57 -> Fe57
  double ni57frac = 0.;
  double co57frac = 0.;
  double fe57frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_NI57, NUCLIDE_CO57, timediff, &ni57frac, &co57frac, &fe57frac_fromdecay);

  // Fe52 -> Mn52 -> Cr52
  double fe52frac = 0.;
  double mn52frac = 0.;
  double cr52frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_FE52, NUCLIDE_MN52, timediff, &fe52frac, &mn52frac, &cr52frac_fromdecay);

  // Cr48 -> V48 -> Ti48
  double cr48frac = 0.;
  double v48frac = 0.;
  double ti48frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_CR48, NUCLIDE_V48, timediff, &cr48frac, &v48frac, &ti48frac_fromdecay);

  // printout("model cell %d, has input radioactive ni56_init %g, co56_init %g, fe52_init %g\n",modelgridindex,ni56_init,co56_init,fe52_in);

  for (int element = nelements - 1; element >= 0; element--)
  {
    const int atomic_number = get_element(element);
    if (atomic_number == 28)
    {
      const double nifrac = get_stable_abund(modelgridindex, atomic_number) + ni56frac + ni57frac;
      modelgrid[modelgridindex].composition[element].abundance = nifrac;
    }
    else if (atomic_number == 27)
    {
      const double cofrac = get_stable_abund(modelgridindex, atomic_number) + co56frac + co57frac;
      modelgrid[modelgridindex].composition[element].abundance = cofrac;
    }
    else if (atomic_number == 26)
    {
      const double fefrac = get_stable_abund(modelgridindex, atomic_number) + fe52frac + fe56frac_fromdecay + fe57frac_fromdecay;
      modelgrid[modelgridindex].composition[element].abundance = fefrac;
    }
    else if (atomic_number == 25)
    {
      const double mnfrac = get_stable_abund(modelgridindex, atomic_number) + mn52frac;
      modelgrid[modelgridindex].composition[element].abundance = mnfrac;
    }
    else if (atomic_number == 24)
    {
      const double crfrac = get_stable_abund(modelgridindex, atomic_number) + cr48frac + cr52frac_fromdecay;
      modelgrid[modelgridindex].composition[element].abundance = crfrac;
    }
    else if (atomic_number == 23)
    {
      const double vfrac = get_stable_abund(modelgridindex, atomic_number) + v48frac;
      modelgrid[modelgridindex].composition[element].abundance = vfrac;
    }
    else if (atomic_number == 22)
    {
      const double tifrac = get_stable_abund(modelgridindex, atomic_number) + ti48frac_fromdecay;
      modelgrid[modelgridindex].composition[element].abundance = tifrac;
    }
  }
  // printout("model cell %d at t_current %g has frac: Ni %g Co %g Fe %g, stable: Ni %g Co %g Fe %g\n",
  //          modelgridindex, t_current,
  //          modelgrid[modelgridindex].composition[get_elementindex(28)].abundance,
  //          modelgrid[modelgridindex].composition[get_elementindex(27)].abundance,
  //          modelgrid[modelgridindex].composition[get_elementindex(26)].abundance,
  //          get_fnistabel(modelgridindex), get_fcostable(modelgridindex), get_ffestable(modelgridindex));
}


double get_decayedenergy_per_ejectamass(const int modelgridindex, const double tstart)
{
  double endecaytot = 0.;
  for (int i = 0; i < DECAYPATH_COUNT; i++)
  {
    endecaytot += get_endecay_per_ejectamass_between_times(modelgridindex, i, 0., tstart);
  }
  return endecaytot;
}
