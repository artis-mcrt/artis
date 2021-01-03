#include "sn3d.h"
#include "atomic.h"
#include "grid.h"
#include "decay.h"

namespace decay
{

__managed__ double arr_nucdecayenergygamma[RADIONUCLIDE_COUNT] = {0};


__host__ __device__
double nucdecayenergygamma(enum radionuclides nuclide_type)
// average energy (erg) per decay in the form of gamma rays
{
  return arr_nucdecayenergygamma[nuclide_type];
}


__host__ __device__
void set_nucdecayenergygamma(enum radionuclides nuclide_type, double value)
{
  arr_nucdecayenergygamma[nuclide_type] = value;
}


__host__ __device__
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


__host__ __device__
double nucdecayenergy(enum radionuclides nuclide_type)
{
  return nucdecayenergygamma(nuclide_type) + nucdecayenergypositrons(nuclide_type);
}


__host__ __device__
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
      ;
  }
  assert(false);
  return -1;
}


__host__ __device__
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
      ;
  }
  assert(false);
  return -1;
}


__host__ __device__
enum radionuclides decayparent(enum decaypathways decaypath)
{
  switch (decaypath)
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
      ;
  }
  assert(false);
  return RADIONUCLIDE_COUNT;
}


__host__ __device__
enum radionuclides decaydaughter(enum decaypathways decaypath)
{
  switch (decaypath)
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
  assert(false);
  return RADIONUCLIDE_COUNT;
}


__host__ __device__
static bool decaypath_is_chain(enum decaypathways decaypath)
// return true if the second nuclide in the decay path is also radioactive
{
  switch (decaypath)
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
      ;
  }
  assert(false);
  return false;
}


__host__ __device__
double sample_decaytime(enum decaypathways decaypath, const double tdecaymin, const double tdecaymax)
{
  double tdecay = -1;
  const bool ischain = decaypath_is_chain(decaypath);
  const double meanlife1 = meanlife(decayparent(decaypath));

  if (!ischain)
  {
    // simple decay from initial abundances, e.g. Ni56 -> Co56
    while (tdecay <= tdecaymin || tdecay >= tdecaymax)
    {
      const double zrand = gsl_rng_uniform_pos(rng);
      tdecay = -meanlife1 * log(zrand);
    }
  }
  else
  {
    // decay of daughter nuclei produced by decay of parent, e.g. Co56 decay from in Ni56 -> Co56 -> Fe56 (no initial contribution)
    const double meanlife2 = meanlife(decaydaughter(decaypath));
    while (tdecay <= tdecaymin || tdecay >= tdecaymax)
    {
      const double zrand = gsl_rng_uniform_pos(rng);
      const double zrand2 = gsl_rng_uniform_pos(rng);
      tdecay = (-meanlife1 * log(zrand)) + (-meanlife2 * log(zrand2));
    }
  }
  return tdecay;
}


__host__ __device__
enum packet_type get_decay_pellet_type(enum decaypathways decaypath, bool *originated_from_positron)
{
  *originated_from_positron = false; // will be changed if necessary before returning
  switch (decaypath)
  {
    case DECAY_NI56:  // Ni56 pellet
      return TYPE_56NI_PELLET;

    case DECAY_NI56_CO56: // Ni56 -> Co56 pellet
    {
      const double zrand = gsl_rng_uniform(rng);
      if (zrand < nucdecayenergygamma(NUCLIDE_CO56) / nucdecayenergy(NUCLIDE_CO56))
      {
        return TYPE_56CO_PELLET;
      }
      else
      {
        *originated_from_positron = true;
        return TYPE_56CO_POSITRON_PELLET;
      }
    }

    case DECAY_FE52:
      return TYPE_52FE_PELLET;

    case DECAY_FE52_MN52:
      return TYPE_52MN_PELLET;

    case DECAY_CR48:
      return TYPE_48CR_PELLET;

    case DECAY_CR48_V48:
      return TYPE_48V_PELLET;

    case DECAY_CO56:
    {
      // Now it is a 56Co pellet, choose whether it becomes a positron
      const double zrand = gsl_rng_uniform(rng);
      if (zrand < nucdecayenergygamma(NUCLIDE_CO56) / nucdecayenergy(NUCLIDE_CO56))
      {
        return TYPE_56CO_PELLET;
      }
      else
      {
        *originated_from_positron = true;
        return TYPE_56CO_POSITRON_PELLET;
      }
    }

    case DECAY_NI57: // Ni57 pellet
    {
      const double zrand = gsl_rng_uniform(rng);
      if (zrand < nucdecayenergygamma(NUCLIDE_NI57) / nucdecayenergy(NUCLIDE_NI57))
      {
        return TYPE_57NI_PELLET;
      }
      else
      {
        *originated_from_positron = true;
        return TYPE_57NI_POSITRON_PELLET;
      }
    }

    case DECAY_NI57_CO57: // Ni57 -> Co57 pellet
      return TYPE_57CO_PELLET;

    case DECAY_CO57: // Co57 pellet
      return TYPE_57CO_PELLET;

    case DECAYPATH_COUNT:
    {
      printout("Problem selecting pellet type\n");
      abort();
    }
  }
  assert(false);
  return TYPE_ESCAPE; // will never reach here, but gcc needs a return value
}


__host__ __device__
static enum radionuclides find_nucparent(enum radionuclides nuclide)
// get the parent nuclide, or if it doesn't have one then the value RADIONUCLIDE_COUNT is used
{
  for (int d = 0; d < DECAYPATH_COUNT; d++)
  {
    enum decaypathways decaypath = (enum decaypathways)(d);
    if (decaypath_is_chain(decaypath))
    {
      enum radionuclides nuc2 = decaydaughter(decaypath);
      if (nuc2 == nuclide)
      {
        return decayparent(decaypath);
      }
    }
  }
  return RADIONUCLIDE_COUNT;
}


__host__ __device__
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


__host__ __device__
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

  const double tdiff = t_current - globals::t_model;
  calculate_double_decay_chain(initabund1, meanlife1, initabund2, meanlife2, tdiff, abund1, abund2, abund3);
}


__host__ __device__
static double get_modelinitradioabund_decayed(
  const int modelgridindex, const enum radionuclides nuclide_type, const double time)
// only allow decays to decrease the nuclide abundance (i.e. don't count increases due to decay of parent)
{
  assert(time >= 0.);
  const double tdiff = time - globals::t_model;
  return get_modelinitradioabund(modelgridindex, nuclide_type) * exp(- tdiff / meanlife(nuclide_type));
}


__host__ __device__
static double get_modelradioabund_at_time(
  const int modelgridindex, const enum radionuclides nuclide_type, const double time)
// get the mass fraction of a nuclide accounting for all decays including those of its parent
// e.g. Co56 abundance may first increase with time due to Ni56 decays, then decease due to Co56 decay
{
  if (time == 0)
  {
    return get_modelinitradioabund(modelgridindex, nuclide_type);
  }
  else if (nuclide_type == FAKE_GAM_LINE_ID)
  {
    return 0.;
  }
  assert(nuclide_type < RADIONUCLIDE_COUNT);
  assert(time >= 0.);

  enum radionuclides nucparent = find_nucparent(nuclide_type);
  if (nucparent == RADIONUCLIDE_COUNT) // no parent exists, so use simple decay formula
  {
    return get_modelinitradioabund_decayed(modelgridindex, nuclide_type, time);
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


__host__ __device__
static double get_endecay_per_ejectamass_at_time(const int mgi, enum decaypathways decaypath, double time)
// decay energy that would be released from time tstart to time infinity from each decaypath
{
  const enum radionuclides nuc1 = decayparent(decaypath);
  if (decaypath_is_chain(decaypath))
  {
    // double decay, e.g. DECAY_NI56_CO56, the decay of Co56 nuclei that were produced from Ni56 decays
    // decays from the second nuclide (e.g. Co56) due to the initial abundance are not counted here

    const enum radionuclides nuc2 = decaydaughter(decaypath);
    assert(nuc2 < RADIONUCLIDE_COUNT);
    const double initabund1 = get_modelinitradioabund(mgi, nuc1);
    const double initabund2 = 0.; // don't count initial abundance
    double abund1;
    double abund2;
    double abund3;
    calculate_double_decay_chain(initabund1, meanlife(nuc1), initabund2, meanlife(nuc2), time, &abund1, &abund2, &abund3);
    return (abund1 / nucmass(nuc1) + abund2 / nucmass(nuc2)) * nucdecayenergy(nuc2);
  }
  else
  {
    // simple decay from initial abundance , e.g. DECAY_NI56 or DECAY_CO56
    return get_modelinitradioabund_decayed(mgi, nuc1, time) / nucmass(nuc1) * nucdecayenergy(nuc1);
  }
}


__host__ __device__
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


__host__ __device__
double get_simtime_endecay_per_ejectamass(const int mgi, enum decaypathways decaypath)
// get the decay energy released during the simulation time
{
#ifdef NO_INITIAL_PACKETS
  return get_endecay_per_ejectamass_between_times(mgi, decaypath, globals::tmin, globals::tmax);
#else
  // allow decays from time zero
  return get_endecay_per_ejectamass_between_times(mgi, decaypath, 0., globals::tmax);
#endif
}


__host__ __device__
double get_decay_power_per_ejectamass(enum decaypathways decaypath, const int modelgridindex, const double time)
// total decay energy injection rate in erg / s / kg
{
  double decaypower = 0.;

  const enum radionuclides nuc1 = decayparent(decaypath);
  if (decaypath_is_chain(decaypath))
  {
    // double decay, e.g. DECAY_NI56_CO56, the decay of Co56 nuclei that were produced from Ni56 decays
    // decays from the second nuclide (e.g. Co56) due to the initial abundance are not counted here

    const enum radionuclides nuc2 = decaydaughter(decaypath);
    assert(nuc2 < RADIONUCLIDE_COUNT);
    const double initabund1 = get_modelinitradioabund(modelgridindex, nuc1);
    const double initabund2 = 0.; // don't count initial abundance
    double abund1;
    double abund2;
    double abund3;
    calculate_double_decay_chain(initabund1, meanlife(nuc1), initabund2, meanlife(nuc2), time, &abund1, &abund2, &abund3);
    decaypower = abund2 / nucmass(nuc2) / meanlife(nuc2) * nucdecayenergy(nuc2);
  }
  else
  {
    // simple decay from initial abundance , e.g. DECAY_NI56 or DECAY_CO56
    decaypower = get_modelinitradioabund_decayed(modelgridindex, nuc1, time) / nucmass(nuc1) / meanlife(nuc1) * nucdecayenergy(nuc1);
  }
  // const double time2 = time * 1.001;
  // const double decaypower2 = get_endecay_per_ejectamass_between_times(modelgridindex, decaypath, time, time2) / (time2 - time);
  // printout("compare decaypath %d answer %g and %g\n", decaypath, decaypower, decaypower2);
  return decaypower;
}


__host__ __device__
double get_positroninjection_rate_density(const int modelgridindex, const double t)
// energy injection rate from positrons in erg / s / cm^3
{
  const double rho = get_rho(modelgridindex);

  double pos_dep_sum = 0.;
  for (int n = 0; n < RADIONUCLIDE_COUNT; n++)
  {
    enum radionuclides nuclide = (enum radionuclides)(n);
    if (nucdecayenergypositrons(nuclide) > 0. && get_modelradioabund_at_time(modelgridindex, nuclide, t) > 0.)
    {
      // printout("positrons coming from nuclide %d en %g abund %g\n", nuclide, nucdecayenergypositrons(nuclide), get_modelradioabund_at_time(modelgridindex, nuclide, t));
      const double decayratefactor = get_modelradioabund_at_time(modelgridindex, nuclide, t) / meanlife(nuclide);
      pos_dep_sum += decayratefactor * nucdecayenergypositrons(nuclide) * rho / nucmass(nuclide);
    }
  }

  return pos_dep_sum;
}


__host__ __device__
void update_abundances(const int modelgridindex, const int timestep, const double t_current)
/// Updates the mass fractions of elements associated with the decay sequence
/// Parameters: - modelgridindex: the grid cell for which to update the abundances
///             - t_current: current time (here mid of current timestep)
{
  assert(!globals::homogeneous_abundances); // no longer supported

  // Ni56 -> Co56 -> Fe56
  // abundances from the input model
  double ni56frac = 0.;
  double co56frac = 0.;
  double fe56frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_NI56, NUCLIDE_CO56, t_current, &ni56frac, &co56frac, &fe56frac_fromdecay);

  // Ni57 -> Co57 -> Fe57
  double ni57frac = 0.;
  double co57frac = 0.;
  double fe57frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_NI57, NUCLIDE_CO57, t_current, &ni57frac, &co57frac, &fe57frac_fromdecay);

  // Fe52 -> Mn52 -> Cr52
  double fe52frac = 0.;
  double mn52frac = 0.;
  double cr52frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_FE52, NUCLIDE_MN52, t_current, &fe52frac, &mn52frac, &cr52frac_fromdecay);

  // Cr48 -> V48 -> Ti48
  double cr48frac = 0.;
  double v48frac = 0.;
  double ti48frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(modelgridindex, NUCLIDE_CR48, NUCLIDE_V48, t_current, &cr48frac, &v48frac, &ti48frac_fromdecay);

  // printout("model cell %d, has input radioactive ni56_init %g, co56_init %g, fe52_init %g\n",modelgridindex,ni56_init,co56_init,fe52_in);

  for (int element = get_nelements() - 1; element >= 0; element--)
  {
    const int atomic_number = get_element(element);
    if (atomic_number == 28)
    {
      const double nifrac = get_stable_abund(modelgridindex, atomic_number) + ni56frac + ni57frac;
      set_elem_abundance(modelgridindex, element, nifrac);
    }
    else if (atomic_number == 27)
    {
      const double cofrac = get_stable_abund(modelgridindex, atomic_number) + co56frac + co57frac;
      set_elem_abundance(modelgridindex, element, cofrac);
    }
    else if (atomic_number == 26)
    {
      const double fefrac = get_stable_abund(modelgridindex, atomic_number) + fe52frac + fe56frac_fromdecay + fe57frac_fromdecay;
      set_elem_abundance(modelgridindex, element, fefrac);
    }
    else if (atomic_number == 25)
    {
      const double mnfrac = get_stable_abund(modelgridindex, atomic_number) + mn52frac;
      set_elem_abundance(modelgridindex, element, mnfrac);
    }
    else if (atomic_number == 24)
    {
      const double crfrac = get_stable_abund(modelgridindex, atomic_number) + cr48frac + cr52frac_fromdecay;
      set_elem_abundance(modelgridindex, element, crfrac);
    }
    else if (atomic_number == 23)
    {
      const double vfrac = get_stable_abund(modelgridindex, atomic_number) + v48frac;
      set_elem_abundance(modelgridindex, element, vfrac);
    }
    else if (atomic_number == 22)
    {
      const double tifrac = get_stable_abund(modelgridindex, atomic_number) + ti48frac_fromdecay;
      set_elem_abundance(modelgridindex, element, tifrac);
    }
  }
  // printout("model cell %d at t_current %g has frac: Ni %g Co %g Fe %g, stable: Ni %g Co %g Fe %g\n",
  //          modelgridindex, t_current,
  //          get_elem_abundance(modelgridinex, get_elementindex(28)),
  //          get_elem_abundance(modelgridinex, get_elementindex(27)),
  //          get_elem_abundance(modelgridinex, get_elementindex(26)),
  //          get_fnistabel(modelgridindex), get_fcostable(modelgridindex), get_ffestable(modelgridindex));
}


__host__ __device__
double get_decayedenergy_per_ejectamass(const int modelgridindex, const double tstart)
{
  double endecaytot = 0.;
  for (int i = 0; i < DECAYPATH_COUNT; i++)
  {
    endecaytot += get_endecay_per_ejectamass_between_times(modelgridindex, (enum decaypathways)(i), 0., tstart);
  }
  return endecaytot;
}

}