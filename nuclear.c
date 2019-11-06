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


// work in progress
// the idea is to only measure decays that take place before tmax (to avoid packet reset)
static double get_modelinitradioabund_decaysintime(const int mgi, enum radionuclides nuclide_type)
{
  switch (nuclide_type)
  {
    case DECAY_NI56:
      return (1 - exp(-tmax / meanlife(NUCLIDE_NI56))) * get_modelinitradioabund(mgi, NUCLIDE_NI56);
    case DECAY_NI56_CO56:
      return (1 - (exp(-tmax / meanlife(NUCLIDE_CO56)) - exp(-tmax / meanlife(NUCLIDE_NI56))) / (meanlife(NUCLIDE_CO56) - meanlife(NUCLIDE_NI56))) * get_modelinitradioabund(mgi, NUCLIDE_NI56);
    default:
      return get_modelinitradioabund(mgi, nuclide_type);
  }
}


double get_decay_energy_per_ejectamass(const int mgi, enum decaypathways decaypath)
{
  switch (decaypath)
  {
    case DECAY_NI56:
      return get_modelinitradioabund(mgi, NUCLIDE_NI56) / nucmass(NUCLIDE_NI56) * nucdecayenergy(NUCLIDE_NI56);
    case DECAY_NI56_CO56:
      return get_modelinitradioabund(mgi, NUCLIDE_NI56) / nucmass(NUCLIDE_NI56) * nucdecayenergy(NUCLIDE_CO56);
    case DECAY_CO56:
      return get_modelinitradioabund(mgi, NUCLIDE_CO56) / nucmass(NUCLIDE_CO56) * nucdecayenergy(NUCLIDE_CO56);
    case DECAY_FE52:
      return get_modelinitradioabund(mgi, NUCLIDE_FE52) / nucmass(NUCLIDE_FE52) * nucdecayenergy(NUCLIDE_FE52);
    case DECAY_FE52_MN52:
      return get_modelinitradioabund(mgi, NUCLIDE_FE52) / nucmass(NUCLIDE_FE52) * nucdecayenergy(NUCLIDE_MN52);
    case DECAY_CR48:
      return get_modelinitradioabund(mgi, NUCLIDE_CR48) / nucmass(NUCLIDE_CR48) * nucdecayenergy(NUCLIDE_CR48);
    case DECAY_CR48_V48:
      return get_modelinitradioabund(mgi, NUCLIDE_CR48) / nucmass(NUCLIDE_CR48) * nucdecayenergy(NUCLIDE_V48);
    case DECAY_NI57:
      return get_modelinitradioabund(mgi, NUCLIDE_NI57) / nucmass(NUCLIDE_NI57) * nucdecayenergy(NUCLIDE_NI57);
    case DECAY_NI57_CO57:
      return get_modelinitradioabund(mgi, NUCLIDE_NI57) / nucmass(NUCLIDE_NI57) * nucdecayenergy(NUCLIDE_CO57);
    case DECAY_CO57:
      return get_modelinitradioabund(mgi, NUCLIDE_CO57) / nucmass(NUCLIDE_CO57) * nucdecayenergy(NUCLIDE_CO57);
    case DECAYPATH_COUNT:
      return 0.;
  }
}


double get_positroninjection_rate_density(const int modelgridindex, const double t)
// in erg / s / cm^3
{
  const double rho = get_rho(modelgridindex);

  // Ni56 -> Co56 -> Fe56
  // abundances from the input model
  double ni56frac = 0.;
  double co56frac = 0.;
  double fe56frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_NI56, NUCLIDE_CO56, t, &ni56frac, &co56frac, &fe56frac_fromdecay);

  const double co56_positron_dep = nucdecayenergypositrons(NUCLIDE_CO56) * co56frac / meanlife(NUCLIDE_CO56) / nucmass(NUCLIDE_CO56) * rho;

  const double ni57frac = get_modelinitradioabund(modelgridindex, NUCLIDE_NI57)  * exp(-t / meanlife(NUCLIDE_NI57));
  const double ni57_positron_dep = nucdecayenergypositrons(NUCLIDE_NI57) * ni57frac / meanlife(NUCLIDE_NI57) / nucmass(NUCLIDE_NI57) * rho;

  const double v48_positron_dep = nucdecayenergypositrons(NUCLIDE_V48)  *
    (exp(-t / meanlife(NUCLIDE_V48)) - exp(-t / meanlife(NUCLIDE_CR48))) /
    (meanlife(NUCLIDE_V48) - meanlife(NUCLIDE_CR48)) * get_modelinitradioabund(modelgridindex, NUCLIDE_CR48) / nucmass(NUCLIDE_CR48) * rho;

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

  const double timediff = t_current - t_model;

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
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_NI57, NUCLIDE_CO57, timediff, &ni56frac, &co56frac, &fe56frac_fromdecay);

  // Fe52 -> Mn52 -> Cr52
  double fe52frac = 0.;
  double mn52frac = 0.;
  double cr52frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_FE52, NUCLIDE_MN52, timediff, &ni56frac, &co56frac, &fe56frac_fromdecay);

  // Cr48 -> V48 -> Ti48
  double cr48frac = 0.;
  double v48frac = 0.;
  double ti48frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_CR48, NUCLIDE_V48, timediff, &ni56frac, &co56frac, &fe56frac_fromdecay);

  // printout("model cell %d, has input radioactive ni56_init %g, co56_init %g, fe52_init %g\n",modelgridindex,ni56_init,co56_init,fe52_in);

  for (int element = nelements-1; element >= 0; element--)
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
  const double factor56ni = get_modelinitradioabund(modelgridindex, NUCLIDE_NI56) / 56 / MH * (-1. / (tstart * (- meanlife(NUCLIDE_CO56) + meanlife(NUCLIDE_NI56))))
    * (- nucdecayenergy(NUCLIDE_NI56) * exp(- tstart / meanlife(NUCLIDE_NI56)) * tstart * meanlife(NUCLIDE_CO56) - nucdecayenergy(NUCLIDE_NI56) * exp(- tstart / meanlife(NUCLIDE_NI56)) * meanlife(NUCLIDE_NI56) * meanlife(NUCLIDE_CO56)
       + nucdecayenergy(NUCLIDE_NI56) * exp(- tstart / meanlife(NUCLIDE_NI56)) * tstart * meanlife(NUCLIDE_NI56) + pow(meanlife(NUCLIDE_NI56), 2) * nucdecayenergy(NUCLIDE_NI56) * exp(- tstart / meanlife(NUCLIDE_NI56))
       - meanlife(NUCLIDE_CO56) * tstart * nucdecayenergy(NUCLIDE_CO56) * exp(- tstart / meanlife(NUCLIDE_CO56)) - pow(meanlife(NUCLIDE_CO56), 2) * nucdecayenergy(NUCLIDE_CO56) * exp(- tstart / meanlife(NUCLIDE_CO56))
       + nucdecayenergy(NUCLIDE_CO56) * tstart * meanlife(NUCLIDE_NI56) * exp(- tstart / meanlife(NUCLIDE_NI56)) + pow(meanlife(NUCLIDE_NI56), 2) * nucdecayenergy(NUCLIDE_CO56) * exp(- tstart / meanlife(NUCLIDE_NI56))
       + nucdecayenergy(NUCLIDE_NI56) * meanlife(NUCLIDE_CO56) * meanlife(NUCLIDE_NI56) - nucdecayenergy(NUCLIDE_NI56) * pow(meanlife(NUCLIDE_NI56), 2) - pow(meanlife(NUCLIDE_NI56), 2) * nucdecayenergy(NUCLIDE_CO56) + nucdecayenergy(NUCLIDE_CO56) * pow(meanlife(NUCLIDE_CO56), 2));

  const double factor56co = get_modelinitradioabund(modelgridindex, NUCLIDE_CO56) / 56 / MH * (1. / (tstart * meanlife(NUCLIDE_CO56)))
    * (meanlife(NUCLIDE_CO56) * tstart * nucdecayenergy(NUCLIDE_CO56) * exp(- tstart / meanlife(NUCLIDE_CO56)) + pow(meanlife(NUCLIDE_CO56), 2) * nucdecayenergy(NUCLIDE_CO56) * exp(- tstart / meanlife(NUCLIDE_CO56)));

  // Ni56 -> Co56 -> Fe56
  // abundances from the input model
  // const double ni56_init = get_modelinitradioabund(modelgridindex, NUCLIDE_NI56);
  // double ni56frac = 0.;
  // double co56frac = 0.;
  // double fe56frac_fromdecay = 0.;
  // calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_NI56, NUCLIDE_CO56, tstart, &ni56frac, &co56frac, &fe56frac_fromdecay);
  // const double factor56ni_new = (ni56_init - ni56frac) * nucdecayenergy(NUCLIDE_NI56) / nucmass(NUCLIDE_NI56);
  // const double factor56co_new = fe56frac_fromdecay * nucdecayenergy(NUCLIDE_CO56) / nucmass(NUCLIDE_CO56);

  const double factor57ni = get_modelinitradioabund(modelgridindex, NUCLIDE_NI57) / 57 / MH * (-1. / (tstart * (- meanlife(NUCLIDE_CO57) + meanlife(NUCLIDE_NI57))))
    * (- nucdecayenergy(NUCLIDE_NI57) * exp(- tstart / meanlife(NUCLIDE_NI57)) * tstart * meanlife(NUCLIDE_CO57) - nucdecayenergy(NUCLIDE_NI57) * exp(- tstart / meanlife(NUCLIDE_NI57)) * meanlife(NUCLIDE_NI57) * meanlife(NUCLIDE_CO57)
       + nucdecayenergy(NUCLIDE_NI57) * exp(- tstart / meanlife(NUCLIDE_NI57)) * tstart * meanlife(NUCLIDE_NI57) + pow(meanlife(NUCLIDE_NI57), 2) * nucdecayenergy(NUCLIDE_NI57) * exp(- tstart / meanlife(NUCLIDE_NI57))
       - meanlife(NUCLIDE_CO57) * tstart * nucdecayenergy(NUCLIDE_CO57) * exp(- tstart / meanlife(NUCLIDE_CO57)) - pow(meanlife(NUCLIDE_CO57), 2) * nucdecayenergy(NUCLIDE_CO57) * exp(- tstart / meanlife(NUCLIDE_CO57))
       + nucdecayenergy(NUCLIDE_CO57) * tstart * meanlife(NUCLIDE_NI57) * exp(- tstart / meanlife(NUCLIDE_NI57)) + pow(meanlife(NUCLIDE_NI57), 2) * nucdecayenergy(NUCLIDE_CO57) * exp(- tstart / meanlife(NUCLIDE_NI57))
       + nucdecayenergy(NUCLIDE_NI57) * meanlife(NUCLIDE_CO57) * meanlife(NUCLIDE_NI57) - nucdecayenergy(NUCLIDE_NI57) * pow(meanlife(NUCLIDE_NI57), 2) - pow(meanlife(NUCLIDE_NI57), 2) * nucdecayenergy(NUCLIDE_CO57) + nucdecayenergy(NUCLIDE_CO57) * pow(meanlife(NUCLIDE_CO57), 2));

  const double factor52fe = get_modelinitradioabund(modelgridindex, NUCLIDE_FE52) / 52 / MH * (-1. / (tstart * (- meanlife(NUCLIDE_MN52) + meanlife(NUCLIDE_FE52))))
    * (- nucdecayenergy(NUCLIDE_FE52) * exp(- tstart / meanlife(NUCLIDE_FE52)) * tstart * meanlife(NUCLIDE_MN52) - nucdecayenergy(NUCLIDE_FE52) * exp(- tstart / meanlife(NUCLIDE_FE52)) * meanlife(NUCLIDE_FE52) * meanlife(NUCLIDE_MN52)
       + nucdecayenergy(NUCLIDE_FE52) * exp(- tstart / meanlife(NUCLIDE_FE52)) * tstart * meanlife(NUCLIDE_FE52) + pow(meanlife(NUCLIDE_FE52), 2) * nucdecayenergy(NUCLIDE_FE52) * exp(- tstart / meanlife(NUCLIDE_FE52))
       - meanlife(NUCLIDE_MN52) * tstart * nucdecayenergy(NUCLIDE_MN52) * exp(- tstart / meanlife(NUCLIDE_MN52)) - pow(meanlife(NUCLIDE_MN52), 2) * nucdecayenergy(NUCLIDE_MN52) * exp(- tstart / meanlife(NUCLIDE_MN52))
       + nucdecayenergy(NUCLIDE_MN52) * tstart * meanlife(NUCLIDE_FE52) * exp(- tstart / meanlife(NUCLIDE_FE52)) + pow(meanlife(NUCLIDE_FE52), 2) * nucdecayenergy(NUCLIDE_MN52) * exp(- tstart / meanlife(NUCLIDE_FE52))
       + nucdecayenergy(NUCLIDE_FE52) * meanlife(NUCLIDE_MN52) * meanlife(NUCLIDE_FE52) - nucdecayenergy(NUCLIDE_FE52) * pow(meanlife(NUCLIDE_FE52), 2) - pow(meanlife(NUCLIDE_FE52), 2) * nucdecayenergy(NUCLIDE_MN52) + nucdecayenergy(NUCLIDE_MN52) * pow(meanlife(NUCLIDE_MN52), 2));

  const double factor48cr = get_modelinitradioabund(modelgridindex, NUCLIDE_CR48) / 48 / MH * (-1. / (tstart * (- meanlife(NUCLIDE_V48) + meanlife(NUCLIDE_CR48))))
    * (- nucdecayenergy(NUCLIDE_CR48) * exp(- tstart / meanlife(NUCLIDE_CR48)) * tstart * meanlife(NUCLIDE_V48) - nucdecayenergy(NUCLIDE_CR48) * exp(- tstart / meanlife(NUCLIDE_CR48)) * meanlife(NUCLIDE_CR48) * meanlife(NUCLIDE_V48)
       + nucdecayenergy(NUCLIDE_CR48) * exp(- tstart / meanlife(NUCLIDE_CR48)) * tstart * meanlife(NUCLIDE_CR48) + pow(meanlife(NUCLIDE_CR48), 2) * nucdecayenergy(NUCLIDE_CR48) * exp(- tstart / meanlife(NUCLIDE_CR48))
       - meanlife(NUCLIDE_V48) * tstart * nucdecayenergy(NUCLIDE_V48) * exp(- tstart / meanlife(NUCLIDE_V48)) - pow(meanlife(NUCLIDE_V48), 2) * nucdecayenergy(NUCLIDE_V48) * exp(- tstart / meanlife(NUCLIDE_V48))
       + nucdecayenergy(NUCLIDE_V48) * tstart * meanlife(NUCLIDE_CR48) * exp(- tstart / meanlife(NUCLIDE_CR48)) + pow(meanlife(NUCLIDE_CR48), 2) * nucdecayenergy(NUCLIDE_V48) * exp(- tstart / meanlife(NUCLIDE_CR48))
       + nucdecayenergy(NUCLIDE_CR48) * meanlife(NUCLIDE_V48) * meanlife(NUCLIDE_CR48) - nucdecayenergy(NUCLIDE_CR48) * pow(meanlife(NUCLIDE_CR48), 2) - pow(meanlife(NUCLIDE_CR48), 2) * nucdecayenergy(NUCLIDE_V48) + nucdecayenergy(NUCLIDE_V48) * pow(meanlife(NUCLIDE_V48), 2));


  double entot = 0.;
  entot += factor56ni;
  entot += factor56co;
  entot += factor57ni;
  entot += factor52fe;
  entot += factor48cr;

  return entot;
}