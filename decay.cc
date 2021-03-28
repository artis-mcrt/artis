#include <algorithm> // std::max
#include <vector>

#include "sn3d.h"
#include "atomic.h"
#include "gamma.h"
#include "grid.h"
#include "nonthermal.h"

#include "decay.h"

namespace decay
{

struct nuclide {
  int z;                     // atomic number
  int a;                     // mass number
  double meanlife;           // mean lifetime before decay [s]
  double endecay_positrons;  // average energy per decay in kinetic energy of emitted positrons [erg]
  double endecay_gamma;      // average energy per decay in gamma rays [erg]
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


static bool nuc_is_parent(const int z_parent, const int a_parent, const int z, const int a)
// check if (z_parent, a_parent) is a parent of (z, a)
{
  // electron capture/beta plus decay only for now
  return nuc_exists(z_parent, a_parent) && (a_parent == a) && z_parent == z + 1;
}

static void print_chain(std::vector<int> &z_list, std::vector<int> &a_list)
{
  assert_always(z_list.size() == a_list.size());
  if (z_list.size() > 0)
  {
    printout(" chain (len %lu): (%d,%d)", z_list.size(), z_list[0], a_list[0]);
  }

  for (size_t i = 1; i < z_list.size(); i++)
  {
    printout(" -> (%d,%d)", z_list[i], a_list[i]);
  }
  printout("\n");
}


static void find_chains(void)
{
  int maxchainlength = 0;
  int numchains = 0;
  for (int endnuc = 0; endnuc < get_num_nuclides(); endnuc++)
  {
    if (get_nuc_z(endnuc) < 1) // FAKE_GAM_LINE_ID
    {
      continue;
    }

    std::vector<int> z_list = {get_nuc_z(endnuc)};
    std::vector<int> a_list = {get_nuc_a(endnuc)};
    maxchainlength = std::max(maxchainlength, (int) z_list.size());
    numchains++;

    print_chain(z_list, a_list);

    bool endofchain = false;
    while (!endofchain)
    {
      endofchain = true; // start true and set to false if we find a parent
      for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
      {
        if (nuc_is_parent(get_nuc_z(nucindex), get_nuc_a(nucindex), z_list[0], a_list[0]))
        {
          z_list.insert(z_list.begin(), get_nuc_z(nucindex));
          a_list.insert(a_list.begin(), get_nuc_a(nucindex));
          maxchainlength = std::max(maxchainlength, (int) z_list.size());
          numchains++;

          print_chain(z_list, a_list);
          endofchain = false;
        }
      }
    }
  }
  printout("Number of chains: %d (max length %d)\n", numchains, maxchainlength);
}


__host__ __device__
void init_nuclides(void)
{
  // all decays are currently assumed to be electron-capture or beta+ (Z -> Z - 1)

  num_nuclides = 9;
  nuclides = (struct nuclide *) calloc(num_nuclides, sizeof(struct nuclide));
  assert_always(nuclides != NULL);

  for (int nucindex = 0; nucindex < num_nuclides; nucindex++)
  {
    nuclides[0].z = -1;
    nuclides[0].a = -1;
    nuclides[0].meanlife = -1;
    nuclides[0].endecay_positrons = 0.;
    nuclides[0].endecay_gamma = 0.;
  }

  nuclides[0].z = 28; // Ni57
  nuclides[0].a = 57;
  nuclides[0].meanlife = 51.36 * 60;
  nuclides[0].endecay_positrons = 0.354 * MEV * 0.436;

  nuclides[1].z = 28; // Ni56
  nuclides[1].a = 56;
  nuclides[1].meanlife = 8.80 * DAY;

  nuclides[2].z = 27; // Co56
  nuclides[2].a = 56;
  nuclides[2].meanlife = 113.7 * DAY;
  nuclides[2].endecay_positrons = 0.63 * MEV * 0.19;

  nuclides[3].z = -1;  // FAKE_GAM_LINE_ID
  nuclides[3].a = -1;

  nuclides[4].z = 24; // Cr48
  nuclides[4].a = 48;
  nuclides[4].meanlife = 1.29602 * DAY;

  nuclides[5].z = 23; // V48
  nuclides[5].a = 48;
  nuclides[5].meanlife = 23.0442 * DAY;
  nuclides[5].endecay_positrons = 0.290 * 0.499 * MEV;

  nuclides[6].z = 27; // Co57
  nuclides[6].a = 57;
  nuclides[6].meanlife = 392.03 * DAY;

  nuclides[7].z = 26; // Fe52
  nuclides[7].a = 52;
  nuclides[7].meanlife = 0.497429 * DAY;

  nuclides[8].z = 25; // Mn52
  nuclides[8].a = 52;
  nuclides[8].meanlife = 0.0211395 * DAY;

  printout("init_nuclides() done. num_nuclides %d\n", get_num_nuclides());

  /// Read in data for gamma ray lines and make a list of them in energy order.
  init_gamma_linelist();

  find_chains();
}


__host__ __device__
double nucdecayenergygamma(int z, int a)
// average energy (erg) per decay in the form of gamma rays
{
  return nuclides[get_nuc_index(z, a)].endecay_gamma;
}


__host__ __device__
void set_nucdecayenergygamma(int z, int a, double value)
// average energy per decay in the form of gamma rays [erg]
{
  nuclides[get_nuc_index(z, a)].endecay_gamma = value;
}


__host__ __device__
double nucdecayenergypositrons(int z, int a)
// average energy (erg) per decay in the form of positron kinetic energy [erg]
{
  return nuclides[get_nuc_index(z, a)].endecay_positrons;
}


__host__ __device__
double nucdecayenergy(int z, int a)
// energy release per decay [erg]
{
  return nucdecayenergygamma(z, a) + nucdecayenergypositrons(z, a);
}


__host__ __device__
double get_meanlife(int z, int a)
{
  assert_always(z > 0);
  assert_always(a >= z);
  if (!nuc_exists(z, a))
  {
    assert_always(nuc_exists(z + 1, a));
    return -1;
  }
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
static double calculate_decaychain_abund(
  const double firstinitabund, const double *meanlifetimes, const int num_nuclides, const double time)
{
  // calculate final abundance from multiple decays, e.g., Ni56 -> Co56 -> Fe56 (nuc[0] -> nuc[1] -> nuc[2])
  // the top nuclide abundance is given and the end abundance is returned (all intermediates are assumed to start
  // with zero abundance)
  // note: first and last can be nuclide can be the same if num_nuclides==1, giving simple decay formula
  // meanlifetimes: array of mean lifetimes from nuc[0]
  // assuming intermediate nuclides start with no abundance, return the abundance at the end of the chain

  assert_always(num_nuclides >= 1);

  // if the meanlife is zero or negative, that indicates a stable nuclide

  // last nuclide might be stable (meanlife <= 0.)
  double lambdas[num_nuclides];
  for (int i = 0; i < num_nuclides; i++)
  {
    assert_always(meanlifetimes[i] > 0. || (i == num_nuclides - 1)); // only the last nuclide can be stable
    lambdas[i] = (meanlifetimes[i] > 0.) ? 1. / meanlifetimes[i] : 0.;
  }

  double lambdaproduct = 1.;
  for (int j = 0; j < num_nuclides - 1; j++)
  {
    lambdaproduct *= lambdas[j];
  }

  double sum = 0;
  for (int j = 0; j < num_nuclides; j++)
  {
    double denominator = 1.;
    for (int p = 0; p < num_nuclides; p++)
    {
      if (p != j)
      {
        denominator *= (lambdas[p] - lambdas[j]);
      }
    }
    sum += exp(-lambdas[j] * time) / denominator;
  }

  const double lastabund = firstinitabund * lambdaproduct * sum;

  return lastabund;
}


__host__ __device__
static void calculate_double_decay_chain(
  const double initabund1, const double meanlife1,
  const double initabund2, const double meanlife2,
  const double t_afterinit,
  double *abund1, double *abund2, double *abund3)
{
  // calculate abundances from double decay, e.g., Ni56 -> Co56 -> Fe56 (nuc1 -> nuc2 -> nuc3)
  // initabund3 is assumed to be zero, so the abundance of species 3 is only from decays of species 2

  double meanlifetimes[3];
  meanlifetimes[0] = meanlife1;
  meanlifetimes[1] = meanlife2;
  meanlifetimes[2] = -1.;  // stable nuclide

  *abund1 = calculate_decaychain_abund(initabund1, meanlifetimes, 1, t_afterinit);

  *abund2 = (
    calculate_decaychain_abund(initabund1, meanlifetimes, 2, t_afterinit) +
    calculate_decaychain_abund(initabund2, &meanlifetimes[1], 1, t_afterinit));

  *abund3 = (
    calculate_decaychain_abund(initabund1, meanlifetimes, 3, t_afterinit) +
    calculate_decaychain_abund(initabund2, &meanlifetimes[1], 2, t_afterinit));

  // ensure that the decays haven't altered the total abundance of all three species
  assert_always(fabs((initabund1 + initabund2) - (*abund1 + *abund2 + *abund3)) < 0.001);
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


static double get_ancestor_abundcontrib(
  const double modelgridindex, const int *const z_list, const int *const a_list, const int chainlength, const double time)
// get the contribution to the abundance of (z_list[-1] a_list[-1]) from the initial abundance of (z_list[0], a_list[0])
// with chained decays (recursive, with items appended to the start of the z_list and a_list)
// chainlength can be one, in which case we get the deca
{
  assert_always(chainlength > 0);
  double meanlifetimes[chainlength];
  for (int i = 0; i < chainlength; i++)
  {
    meanlifetimes[i] = get_meanlife(z_list[i], a_list[i]);
  }

  const double t_afterinit = time - globals::t_model;

  double abundcontrib = 0.;

  const int z_top = z_list[0]; // current top of chain considered (but might be a partial chain)
  const int a_top = a_list[0];

  // add the decayed contribution from the current top of the chain (if it's in the radioactive list)
  // we're only counting ancestor contributions, so don't add for the single-nuclide 'chain'
  if (nuc_exists(z_top, a_top) && chainlength > 1)
  {
    const double top_abund = get_modelinitradioabund(modelgridindex, z_top, a_top);
    abundcontrib += calculate_decaychain_abund(top_abund, meanlifetimes, chainlength, t_afterinit);
  }

  // find parent nulides to append to the top of the chain
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
  {
    const int z_check = get_nuc_z(nucindex);
    const int a_check = get_nuc_a(nucindex);
    if (nuc_is_parent(z_check, a_check, z_top, a_top))
    {
      // found a parent nuclide of the current top nuclide, so add it to the chain
      int zlist_plusone[chainlength + 1];
      int alist_plusone[chainlength + 1];
      for (int i = 0; i < chainlength; i++)
      {
        zlist_plusone[i + 1] = z_list[i];
        alist_plusone[i + 1] = a_list[i];
      }
      zlist_plusone[0] = z_check;
      alist_plusone[0] = a_check;

      abundcontrib += get_ancestor_abundcontrib(modelgridindex, zlist_plusone, alist_plusone, chainlength + 1, time);
    }
  }

  return abundcontrib;
}


__host__ __device__
static double get_modelradioabund_at_time(
  const int modelgridindex, const int z, const int a, const double time)
// Get the mass fraction of a nuclide accounting for all decays including those of its parent and grandparent.
// e.g., Co56 abundance may first increase with time due to Ni56 decays, then decease due to Co56 decay
// Can be called for stable nuclides that are one daughters of the radioactive nuclide list e.g., Fe56
// For stable nuclides, abundance returned only comes from other decays (some could be included in init model elem frac)
{
  if (z < 1) // skip FAKE_GAM_LINE_ID
  {
    return 0.;
  }
  assert_always(time >= 0.);

  // start with contribution of decays from ancestor nuclides and add decayed initial abundance (if present)
  double abund = get_ancestor_abundcontrib(modelgridindex, &z, &a, 1, time);

  if (nuc_exists(z, a))  // stable nuclide e.g. Fe56 will not exist in the list or have an init abundance
  {
    abund += get_modelinitradioabund_decayed(modelgridindex, z, a, time);
  }

  return abund;
}


__host__ __device__
static double get_endecay_per_ejectamass_at_time(
  const int mgi, const bool from_parent_abund, const int z, const int a, const double time)
// returns decay energy [erg] that would be released from time tstart [s] to time infinity a given decaypath
{
  if (from_parent_abund)
  {
    // e.g. DECAY_NI56_CO56, represents the decay of Co56 nuclei
    // that were produced from Ni56 in the initial abundance.
    // Decays from Co56 due to the initial abundance of Co56 are not counted here,
    // nor is the energy from Ni56 decays

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
    // double decay, e.g. Ni56 -> Co56 -> Fe56, the decay of Co56 nuclei that were produced from Ni56 decays
    // decays of the second nuclide (e.g. Co56) due to its initial abundance are not counted here
    if (!nuc_exists(z + 1, a))
    {
      return 0.;
    }
    assert(!nuc_exists(z + 2, a)); // no more than three-nuclide chains for now

    const double initabund1 = get_modelinitradioabund(modelgridindex, z + 1, a);
    const double initabund2 = 0.; // don't count initial abundance
    double abund1;
    double abund3;
    calculate_double_decay_chain(initabund1, get_meanlife(z + 1, a), initabund2, get_meanlife(z, a), time, &abund1, &nucabund, &abund3);
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

  printout("update_abundances for cell %d timestep %d\n", modelgridindex, timestep);

  for (int element = get_nelements() - 1; element >= 0; element--)
  {
    const int atomic_number = get_element(element);
    double isofracsum = 0.; // mass fraction sum of radioactive isotopes, and stable nuclei coming from other decays
    for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
    {
      const int a = get_nuc_a(nucindex);
      if (get_nuc_z(nucindex) == atomic_number)
      {
        // radioactive isotope of this element
        isofracsum += get_modelradioabund_at_time(modelgridindex, atomic_number, a, t_current);
      }
      else if (!nuc_exists(atomic_number, a) && get_nuc_z(nucindex) == atomic_number + 1)
      {
        // nuclide is one step off the network (only includes radioactive nuclides), e.g. Fe56
        // note: there could also be Fe56 included in stable_initabund(z), but
        // here we only count the contribution from decays
        isofracsum += get_modelradioabund_at_time(modelgridindex, atomic_number, a, t_current);
      }
    }

    const double elmassfrac = get_stable_initabund(modelgridindex, element) + isofracsum;
    set_elem_abundance(modelgridindex, element, elmassfrac);
  }

  nonthermal::calculate_deposition_rate_density(modelgridindex, timestep);
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

  pkt_ptr->type = TYPE_RADIOACTIVE_PELLET;
  pkt_ptr->pellet_nucindex = nucindex;

  const double zrand = gsl_rng_uniform(rng);
  pkt_ptr->originated_from_positron = (zrand >= nucdecayenergygamma(z, a) / nucdecayenergy(z, a));
}

}