#include "sn3d.h"
#include "atomic.h"
#include "gamma.h"
#include "grid.h"
#include "decay.h"

namespace decay
{

const int NUCLIDE_NI57 = 0;
const int NUCLIDE_NI56 = 1;
const int NUCLIDE_CO56 = 2;
const int NUCLIDE_CR48 = 4;
const int NUCLIDE_V48 = 5;
const int NUCLIDE_CO57 = 6;
const int NUCLIDE_FE52 = 7;
const int NUCLIDE_MN52 = 8;

enum decaypathways {
  DECAY_NI56 = 0,
  DECAY_NI56_CO56 = 1,
  DECAY_FE52 = 2,
  DECAY_FE52_MN52 = 3,
  DECAY_CR48 = 4,
  DECAY_CR48_V48 = 5,
  DECAY_CO56 = 6,
  DECAY_NI57 = 7,
  DECAY_NI57_CO57 = 8,
  DECAY_CO57 = 9,
  DECAYPATH_COUNT = 10,
};


struct nuclide {
  int z;
  int a;
  double meanlife;
  double endecay_positrons;
  double endecay_gamma;
};

struct nuclide *nuclides = NULL;
int num_nuclides = 0;


__host__ __device__
int get_num_nuclides(void)
{
  assert_always(num_nuclides > 0);
  return num_nuclides;
}


__host__ __device__
int get_nuc_z(int nucindex)
{
  assert_always(nucindex >= 0);
  assert_always(nucindex < get_num_nuclides());
  return nuclides[nucindex].z;
}


__host__ __device__
int get_nuc_a(int nucindex)
{
  assert_always(nucindex >= 0);
  assert_always(nucindex < get_num_nuclides());
  return nuclides[nucindex].a;
}


__host__ __device__
int get_nuc_index(int z, int a)
// get the nuclide array index from the atomic number and mass number
{
  assert_always(get_num_nuclides() > 0);
  assert_always(nuclides != NULL);

  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
  {
    if (nuclides[nucindex].z == z && nuclides[nucindex].a == a)
    {
      return nucindex;
    }
  }
  printout("Could not find nuclide Z=%d A=%d\n", z, a);
  assert_always(false); // nuclide not found
  return -1;
}


__host__ __device__
static bool nuc_exists(int z, int a)
// check if nuclide exists in the simulation
{
  assert_always(get_num_nuclides() > 0);
  assert_always(nuclides != NULL);

  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
  {
    if (nuclides[nucindex].z == z && nuclides[nucindex].a == a)
    {
      return true;
    }
  }
  return false;
}


__host__ __device__
void init_nuclides(void)
{
  num_nuclides = 9;
  nuclides = (struct nuclide *) calloc(num_nuclides, sizeof(struct nuclide));
  assert_always(nuclides != NULL);

  nuclides[0].z = 28; // Ni
  nuclides[0].a = 57;
  nuclides[0].meanlife = 51.36 * 60;
  nuclides[0].endecay_positrons = 0.354 * MEV * 0.436;

  nuclides[1].z = 28;
  nuclides[1].a = 56;
  nuclides[1].meanlife = 8.80 * DAY;

  nuclides[2].z = 27; // Co
  nuclides[2].a = 56;
  nuclides[2].meanlife = 113.7 * DAY;
  nuclides[2].endecay_positrons = 0.63 * MEV * 0.19;

  nuclides[3].z = -1;  // FAKE_GAM_LINE_ID
  nuclides[3].a = -1;

  nuclides[4].z = 24; // Cr
  nuclides[4].a = 48;
  nuclides[4].meanlife = 1.29602 * DAY;

  nuclides[5].z = 23; // V
  nuclides[5].a = 48;
  nuclides[5].meanlife = 23.0442 * DAY;
  nuclides[5].endecay_positrons = 0.290 * 0.499 * MEV;

  nuclides[6].z = 27; // Co
  nuclides[6].a = 57;
  nuclides[6].meanlife = 392.03 * DAY;

  nuclides[7].z = 26; // Fe
  nuclides[7].a = 52;
  nuclides[7].meanlife = 0.497429 * DAY;

  nuclides[8].z = 25; // Mn
  nuclides[8].a = 52;
  nuclides[8].meanlife = 0.0211395 * DAY;

  printout("init_nuclides() done. num_nuclides %d\n", get_num_nuclides());

  /// Read in data for gamma ray lines and make a list of them in energy order.
  init_gamma_linelist();
}


__host__ __device__
double nucdecayenergygamma(int z, int a)
// average energy (erg) per decay in the form of gamma rays
{
  return nuclides[get_nuc_index(z, a)].endecay_gamma;
}


__host__ __device__
void set_nucdecayenergygamma(int z, int a, double value)
{
  nuclides[get_nuc_index(z, a)].endecay_gamma = value;
}


__host__ __device__
double nucdecayenergypositrons(int z, int a)
// average energy (erg) per decay in the form of positrons
{
  return nuclides[get_nuc_index(z, a)].endecay_positrons;
}


__host__ __device__
double nucdecayenergy(int z, int a)
{
  return nucdecayenergygamma(z, a) + nucdecayenergypositrons(z, a);
}


__host__ __device__
double get_meanlife(int z, int a)
{
  assert_always(z > 0);
  assert_always(a >= z);
  const int nucindex = get_nuc_index(z, a);
  return nuclides[nucindex].meanlife;
}


__host__ __device__
double nucmass(int z, int a)
{
  assert_always(z > 0);
  assert_always(a >= z);

  return a * MH;

  // const int nucindex = get_nuc_index(z, a);
  // return nuclides[nucindex].amass;
}


__host__ __device__
static double sample_decaytime(bool from_parent_abund, int z, int a, const double tdecaymin, const double tdecaymax)
{
  double tdecay = -1;
  const double meanlife = get_meanlife(z, a);

  if (!from_parent_abund)
  {
    // simple decay from initial abundances, e.g. Ni56 -> Co56
    while (tdecay <= tdecaymin || tdecay >= tdecaymax)
    {
      const double zrand = gsl_rng_uniform_pos(rng);
      tdecay = -meanlife * log(zrand);
    }
  }
  else
  {
    // decay nuclei that were produced by decay of parent, e.g. Co56 decay from Ni56 -> Co56 -> Fe56
    // (no initial Co56 contribution)

    assert_always(nuc_exists(z + 1, a));
    const double meanlife_parent = get_meanlife(z + 1, a);
    while (tdecay <= tdecaymin || tdecay >= tdecaymax)
    {
      const double zrand = gsl_rng_uniform_pos(rng);
      const double zrand2 = gsl_rng_uniform_pos(rng);
      tdecay = (-meanlife_parent * log(zrand)) + (-meanlife * log(zrand2));
    }
  }
  return tdecay;
}


__host__ __device__
static enum packet_type get_decay_pellet_type(const int z, const int a, bool *originated_from_positron)
{
  *originated_from_positron = false; // will be changed if necessary before returning
  if (z == 28 && a == 56)
  {
    return TYPE_56NI_PELLET;
  }
  else if (z == 27 && a == 56)
  {
    const double zrand = gsl_rng_uniform(rng);
    if (zrand < nucdecayenergygamma(27, 56) / nucdecayenergy(27, 56))
    {
      return TYPE_56CO_PELLET;
    }
    else
    {
      *originated_from_positron = true;
      return TYPE_56CO_POSITRON_PELLET;
    }
  }
  else if (z == 26 && a == 52)
  {
    return TYPE_52FE_PELLET;
  }
  else if (z == 25 && a == 52)
  {
    return TYPE_52MN_PELLET;
  }
  else if (z == 24 && a == 48)
  {
    return TYPE_48CR_PELLET;
  }
  else if (z == 23 && a == 48)
  {
    return TYPE_48V_PELLET;
  }
  else if (z == 28 && a == 57)
  {
    const double zrand = gsl_rng_uniform(rng);
    if (zrand < nucdecayenergygamma(28, 57) / nucdecayenergy(28, 57))
    {
      return TYPE_57NI_PELLET;
    }
    else
    {
      *originated_from_positron = true;
      return TYPE_57NI_POSITRON_PELLET;
    }
  }

  assert_always(false);
  return TYPE_ESCAPE; // will never reach here, but gcc needs a return value
}


__host__ __device__
static void calculate_double_decay_chain(
  const double initabund1, const double meanlife1,
  const double initabund2, const double meanlife2,
  const double t_current,
  double *abund1, double *abund2, double *abund3)
{
  // calculate abundances from double decay, e.g., Ni56 -> Co56 -> Fe56 (nuc1 -> nuc2 -> nuc3)
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
  assert_always(fabs((initabund1 + initabund2) - (*abund1 + *abund2 + *abund3)) < 0.001);
}


__host__ __device__
static void calculate_doubledecay_modelabund(
  const int modelgridindex,
  const int nuclide1,
  const int nuclide2,
  const double t_current,
  double *abund1, double *abund2, double *abund3)
{
  const int z1 = get_nuc_z(nuclide1);
  const int a = get_nuc_a(nuclide1);
  const int z2 = get_nuc_z(nuclide2);
  const int a2 = get_nuc_a(nuclide2);
  assert_always(a2 == a); // mass number is conserved
  assert_always(z2 == (z1 - 1)); // beta decay p -> n

  const double initabund1 = get_modelinitradioabund(modelgridindex, z1, a);
  const double meanlife1 = get_meanlife(z1, a);
  const double initabund2 = get_modelinitradioabund(modelgridindex, z2, a);
  const double meanlife2 = get_meanlife(z2, a);

  const double tdiff = t_current - globals::t_model;
  calculate_double_decay_chain(initabund1, meanlife1, initabund2, meanlife2, tdiff, abund1, abund2, abund3);
}


__host__ __device__
static double get_modelinitradioabund_decayed(
  const int modelgridindex, const int z, const int a, const double time)
// only allow decays to decrease the nuclide abundance (i.e. don't count increases due to decay of parent)
{
  assert_always(time >= 0.);
  const double tdiff = time - globals::t_model;
  return get_modelinitradioabund(modelgridindex, z, a) * exp(- tdiff / get_meanlife(z, a));
}


__host__ __device__
static double get_modelradioabund_at_time(
  const int modelgridindex, const int z, const int a, const double time)
// get the mass fraction of a nuclide accounting for all decays including those of its parent
// e.g. Co56 abundance may first increase with time due to Ni56 decays, then decease due to Co56 decay
{
  if (time == 0)
  {
    return get_modelinitradioabund(modelgridindex, z, a);
  }
  else if (z < 1) // skip FAKE_GAM_LINE_ID
  {
    return 0.;
  }
  assert_always(time >= 0.);

  const int nucindex = get_nuc_index(z, a);
  if (!nuc_exists(z + 1, a)) // no parent exists, so use simple decay formula
  {
    return get_modelinitradioabund_decayed(modelgridindex, z, a, time);
  }
  else
  {
    // nuclide is part of a double-decay chain, e.g., Co56 in the chain: Ni56 -> Co56 -> Fe56
    const int nucparent = get_nuc_index(z + 1, a);
    assert_always(!nuc_exists(z + 2, a)); // only three-nuclide chains work for now
    double abund1 = 0.;
    double abund2 = 0.;
    double abund3 = 0.;
    calculate_doubledecay_modelabund(modelgridindex, nucparent, nucindex, time, &abund1, &abund2, &abund3);
    return abund2;
  }
}


__host__ __device__
static double get_endecay_per_ejectamass_at_time(
  const int mgi, const bool from_parent_abund, const int z, const int a, const double time)
// returns decay energy [erg] that would be released from time tstart [s] to time infinity a given decaypath
{
  if (from_parent_abund)
  {
    // double decaypath, e.g. DECAY_NI56_CO56, represents the decay of Co56 nuclei
    // that were produced from decays of Ni56 in the initial abundance.
    // Decays from Co56 due to the initial abundance of Co56 are not counted here,
    // nor is the energy from decays of Ni56

    if (!nuc_exists(z + 1, a)) // parent is not included (e.g. Ni56 has no parent)
    {
      return 0.;
    }
    assert(!nuc_exists(z + 2, a)); // no more than three-nuclide chains for now

    const double initabund1 = get_modelinitradioabund(mgi, z + 1, a);
    const double initabund2 = 0.; // don't count initial abundance
    double abund1;
    double abund2;
    double abund3;
    calculate_double_decay_chain(initabund1, get_meanlife(z + 1, a), initabund2, get_meanlife(z, a), time, &abund1, &abund2, &abund3);
    return (abund1 / nucmass(z + 1, a) + abund2 / nucmass(z, a)) * nucdecayenergy(z, a);
  }
  else
  {
    // simple decay from initial abundance , e.g. DECAY_NI56 or DECAY_CO56
    return get_modelinitradioabund_decayed(mgi, z, a, time) / nucmass(z, a) * nucdecayenergy(z, a);
  }
}


__host__ __device__
static double get_endecay_per_ejectamass_between_times(
  const int mgi, const bool from_parent_abund, const int z, const int a, double tlow, double thigh)
// energy per mass [erg/g] released by a decaypath between two times [s]
{
  assert_always(tlow <= thigh);
  const double energy_tlow = get_endecay_per_ejectamass_at_time(mgi, from_parent_abund, a, z, tlow);
  const double energy_thigh = get_endecay_per_ejectamass_at_time(mgi, from_parent_abund, a, z, thigh);
  assert_always(energy_tlow >= energy_thigh);
  return energy_tlow - energy_thigh;
}


__host__ __device__
double get_simtime_endecay_per_ejectamass(const int mgi, const bool from_parent_abund, const int z, const int a)
// get the decay energy released during the simulation time
{
#ifdef NO_INITIAL_PACKETS
  // get decay energy released from t=tmin to tmax
  return get_endecay_per_ejectamass_between_times(mgi, from_parent_abund, a, z, globals::tmin, globals::tmax);
#else
  // get decay energy released from t=0 to tmax
  return get_endecay_per_ejectamass_between_times(mgi, from_parent_abund, a, z, 0., globals::tmax);
#endif
}


__host__ __device__
static double get_decay_power_per_ejectamass(
  const bool from_parent_abund, const int z, const int a, const int modelgridindex, const double time)
// total decay power per mass [erg / s / kg] for a given decaypath
{
  double nucabund = 0.;
  if (from_parent_abund)
  {
    // double decay, e.g. DECAY_NI56_CO56, the decay of Co56 nuclei that were produced from Ni56 decays
    // decays from the second nuclide (e.g. Co56) due to the initial abundance are not counted here
    if (!nuc_exists(z + 1, a))
    {
      return 0.;
    }
    assert(!nuc_exists(z + 2, a)); // no more than three-nuclide chains for now

    const double initabund1 = get_modelinitradioabund(modelgridindex, z + 1, a);
    const double initabund2 = 0.; // don't count initial abundance
    double abund1;
    double abund3;
    calculate_double_decay_chain(initabund1, get_meanlife(z + 1, a), initabund2, get_meanlife(z + 1, a), time, &abund1, &nucabund, &abund3);
  }
  else
  {
    // simple decay from initial abundance , e.g. DECAY_NI56 or DECAY_CO56
    nucabund = get_modelinitradioabund_decayed(modelgridindex, z, a, time);
  }
  const double decaypower = nucabund / nucmass(z, a) / get_meanlife(z, a) * nucdecayenergy(z, a);
  // const double time2 = time * 1.001;
  // const double decaypower2 = get_endecay_per_ejectamass_between_times(modelgridindex, decaypath, time, time2) / (time2 - time);
  // printout("compare decaypath %d answer %g and %g\n", decaypath, decaypower, decaypower2);
  return decaypower;
}


__host__ __device__
double get_modelcell_decay_energy_density(const int mgi)
{
  double modelcell_decay_energy_density = 0.;
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
  {
    const int z = decay::get_nuc_z(nucindex);
    if (z > 0) // skip FAKE_GAM_LINE_ID
    {
      const int a = decay::get_nuc_a(nucindex);
      for (int from_parent_abund = 0; from_parent_abund < 2; from_parent_abund++)
      {
        modelcell_decay_energy_density += (
          get_rhoinit(mgi) * get_simtime_endecay_per_ejectamass(mgi, from_parent_abund, z, a) * MH);
      }
    }
  }
  return modelcell_decay_energy_density;
}


__host__ __device__
double get_positroninjection_rate_density(const int modelgridindex, const double t)
// energy release rate from positrons in erg / s / cm^3
{
  const double rho = get_rho(modelgridindex);

  double pos_dep_sum = 0.;
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
  {
    const int z = decay::get_nuc_z(nucindex);
    const int a = decay::get_nuc_a(nucindex);
    if (nucdecayenergypositrons(z, a) > 0. && get_modelradioabund_at_time(modelgridindex, z, a, t) > 0.)
    {
      // printout("positrons coming from nuclide %d en %g abund %g\n", nuclide, nucdecayenergypositrons(nuclide), get_modelradioabund_at_time(modelgridindex, nuclide, t));
      const double decayratefactor = get_modelradioabund_at_time(modelgridindex, z, a, t) / get_meanlife(z, a);
      pos_dep_sum += decayratefactor * nucdecayenergypositrons(z, a) * rho / nucmass(z, a);
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
  assert_always(!globals::homogeneous_abundances); // no longer supported

  // Ni56 -> Co56 -> Fe56
  // abundances from the input model
  double ni56frac = 0.;
  double co56frac = 0.;
  double fe56frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(
    modelgridindex, NUCLIDE_NI56, NUCLIDE_CO56, t_current, &ni56frac, &co56frac, &fe56frac_fromdecay);

  // Ni57 -> Co57 -> Fe57
  double ni57frac = 0.;
  double co57frac = 0.;
  double fe57frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(
    modelgridindex, NUCLIDE_NI57, NUCLIDE_CO57, t_current, &ni57frac, &co57frac, &fe57frac_fromdecay);

  // Fe52 -> Mn52 -> Cr52
  double fe52frac = 0.;
  double mn52frac = 0.;
  double cr52frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(
    modelgridindex, NUCLIDE_FE52, NUCLIDE_MN52, t_current, &fe52frac, &mn52frac, &cr52frac_fromdecay);

  // Cr48 -> V48 -> Ti48
  double cr48frac = 0.;
  double v48frac = 0.;
  double ti48frac_fromdecay = 0.;
  calculate_doubledecay_modelabund(
    modelgridindex, NUCLIDE_CR48, NUCLIDE_V48, t_current, &cr48frac, &v48frac, &ti48frac_fromdecay);

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


void setup_radioactive_pellet(const double e0, const int mgi, PKT *pkt_ptr)
{
  double cumulative_decay_energy_per_mass[2 * get_num_nuclides()];
  for (int i = 0; i < 2 * get_num_nuclides(); i++)
  {
    // visit each radioactive nuclide twice - once for decays from the initial abundance,
    // and again for decays from the abundance coming from decays of a parent nucleus
    const int nucindex = i / 2;
    const bool from_parent_abund = i % 2;
    const int z = decay::get_nuc_z(nucindex);
    const int a = decay::get_nuc_a(nucindex);
    double simtime_endecay_thispath = 0.;

    if (z > 0) // skip for the FAKE_GAM_LINE_ID
    {
      simtime_endecay_thispath = decay::get_simtime_endecay_per_ejectamass(mgi, from_parent_abund, z, a);
    }

    const double lower_sum = ((i > 0) ? cumulative_decay_energy_per_mass[i - 1] : 0);
    cumulative_decay_energy_per_mass[i] = lower_sum + simtime_endecay_thispath;
  }

  const double zrand_chain = gsl_rng_uniform(rng) * cumulative_decay_energy_per_mass[2 * get_num_nuclides() - 1];

  int nucindex = -1;
  bool from_parent_abund;
  for (int i = 0; i < 2 * get_num_nuclides(); i++)
  {
    if (zrand_chain <= cumulative_decay_energy_per_mass[i])
    {
      nucindex = i / 2;
      from_parent_abund = i % 2;
      break;
    }
  }
  assert_always(nucindex >= 0); // Failed to select pellet

  #ifdef NO_INITIAL_PACKETS
  const double tdecaymin = globals::tmin;
  #else
  const double tdecaymin = 0.; // allow decays before the first timestep
  #endif

  const int z = get_nuc_z(nucindex);
  const int a = get_nuc_a(nucindex);
  if (UNIFORM_PELLET_ENERGIES)
  {
    pkt_ptr->tdecay = decay::sample_decaytime(from_parent_abund, z, a, tdecaymin, globals::tmax);
    pkt_ptr->e_cmf = e0;
  }
  else
  {
    // use uniform decay time distribution (scale the packet energies instead)
    // keeping the pellet decay rate constant will give better statistics at very late times when very little
    // energy is released
    const double zrand = gsl_rng_uniform(rng);
    pkt_ptr->tdecay = zrand * tdecaymin + (1. - zrand) * globals::tmax;

    // we need to scale the packet energy up or down according to decay rate at the randomly selected time.
    // e0 is the average energy per packet for this cell and decaypath, so we scale this up or down
    // according to: decay power at this time relative to the average decay power
    const double avgpower = decay::get_simtime_endecay_per_ejectamass(mgi, from_parent_abund, z, a) / (globals::tmax - tdecaymin);
    pkt_ptr->e_cmf = e0 * decay::get_decay_power_per_ejectamass(from_parent_abund, z, a, mgi, pkt_ptr->tdecay) / avgpower;
    // assert_always(pkt_ptr->e_cmf >= 0);
  }

  bool from_positron;
  pkt_ptr->type = decay::get_decay_pellet_type(z, a, &from_positron); // set the packet tdecay and type
  pkt_ptr->originated_from_positron = from_positron;
}

}