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

const char *elsymbols[119] = {
  "n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
  "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo"};

struct nuclide {
  int z;                     // atomic number
  int a;                     // mass number
  double meanlife;           // mean lifetime before decay [s]
  double endecay_positrons;  // average energy per decay in kinetic energy of emitted positrons [erg]
  double endecay_gamma;      // average energy per decay in gamma rays [erg]
};

struct nuclide *nuclides = NULL;
int num_nuclides = 0;

std::vector<std::vector<int>> decaychains_z;
std::vector<std::vector<int>> decaychains_a;

enum nucmodes {
  MODE_ABUND = 0,  // get the abundance
  MODE_DECAYRATE = 1,  // get the decayrate
  MODE_ABUNDEXPANSION = 2,  // get the abundance with expansion factor
};

__host__ __device__
int get_num_nuclides(void)
{
  assert_always(num_nuclides > 0);
  return num_nuclides;
}


static void printout_nuclidename(const int z, const int a)
{
  printout("%s-%d", elsymbols[z], a);
}


static const char *get_elname(const int z)
{
  return elsymbols[z];
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


static int decay_daughter_z(const int z_parent, const int a_parent)
// check if (z_parent, a_parent) is a parent of (z, a)
{
  assert_always(nuc_exists(z_parent, a_parent));
  // electron capture/beta plus decay only for now
  return z_parent - 1;
}


static int decay_daughter_a(const int z_parent, const int a_parent)
// check if (z_parent, a_parent) is a parent of (z, a)
{
  assert_always(nuc_exists(z_parent, a_parent));
  // electron capture/beta +/- keep mass number constant (no alpha decay or fission yet)
  return a_parent;
}


static bool nuc_is_parent(const int z_parent, const int a_parent, const int z, const int a)
// check if (z_parent, a_parent) is a parent of (z, a)
{
  assert_always(nuc_exists(z_parent, a_parent));
  // each radioactive nuclide is limited to one daughter nuclide
  return (decay_daughter_z(z_parent, a_parent) == z && decay_daughter_a(z_parent, a_parent) == a);
}


static void print_chain(const int chainindex)
{
  assert_always(decaychains_z[chainindex].size() == decaychains_a[chainindex].size());
  if (decaychains_z.size() > 0)
  {
    printout(" decay chain %d: ", chainindex);
    printout_nuclidename(decaychains_z[chainindex][0], decaychains_a[chainindex][0]);
  }

  for (size_t i = 1; i < decaychains_z[chainindex].size(); i++)
  {
    printout(" -> ");
    printout_nuclidename(decaychains_z[chainindex][i], decaychains_a[chainindex][i]);
  }
  printout("\n");
}


static void add_ancestorchains(const int z, const int a, const int startchainindex)
{
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
  {
    const int z_parent = get_nuc_z(nucindex); // possible parent, will check if true
    const int a_parent = get_nuc_a(nucindex);
    // printout("z_parent %d a_parent %d isparent(%d, %d) %d\n", z_parent, a_parent, z_list[0], a_list[0], nuc_is_parent(z_parent, a_parent, z_list[0], a_list[0]));
    if (nuc_is_parent(z_parent, a_parent, z, a))
    {
      std::vector<int> new_z_list(decaychains_z[startchainindex]);
      std::vector<int> new_a_list(decaychains_a[startchainindex]);

      // check for repeated nuclides, which would indicate a loop in the decay chain
      for (size_t i = 1; i < new_z_list.size(); i++)
      {
        if (new_z_list[i] == z_parent && new_a_list[i] == a_parent)
        {
          printout("\nERROR: Loop found in nuclear decay chain.\n");
          abort();
        }
      }
      new_z_list.insert(new_z_list.begin(), z_parent);
      new_a_list.insert(new_a_list.begin(), a_parent);
      decaychains_z.push_back(new_z_list);
      decaychains_a.push_back(new_a_list);

      add_ancestorchains(z_parent, a_parent, decaychains_z.size() - 1);
    }
  }
}


static void find_chains(void)
{
  for (int endnuc = 0; endnuc < get_num_nuclides(); endnuc++)
  {
    if (get_nuc_z(endnuc) < 1) // FAKE_GAM_LINE_ID
    {
      continue;
    }

    std::vector<int> z_list = {get_nuc_z(endnuc)};
    std::vector<int> a_list = {get_nuc_a(endnuc)};
    decaychains_z.push_back(z_list);
    decaychains_a.push_back(a_list);

    add_ancestorchains(get_nuc_z(endnuc), get_nuc_a(endnuc), decaychains_z.size() - 1);
  }
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
    nuclides[nucindex].z = -1;
    nuclides[nucindex].a = -1;
    nuclides[nucindex].meanlife = -1;
    nuclides[nucindex].endecay_positrons = 0.;
    nuclides[nucindex].endecay_gamma = 0.;
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

  printout("init_nuclides: num_nuclides %d\n", get_num_nuclides());

  /// Read in data for gamma ray lines and make a list of them in energy order.
  init_gamma_linelist();

  find_chains();

  int maxchainlength = 0;
  for (size_t chainindex = 0; chainindex < decaychains_z.size(); chainindex++)
  {
    print_chain(chainindex);
    maxchainlength = std::max(maxchainlength, (int) decaychains_a[chainindex].size());
  }
  printout("Number of chains: %d (max length %d)\n", (int) decaychains_z.size(), maxchainlength);

  printout("decayenergy(NI56), decayenergy(CO56), decayenergy_gamma(CO56): %g, %g, %g\n",
           nucdecayenergy(28, 56) / MEV, nucdecayenergy(27, 56) / MEV,
           nucdecayenergygamma(27, 56) / MEV);
  printout("decayenergy(NI57), decayenergy_gamma(NI57), nucdecayenergy(CO57): %g, %g, %g\n",
           nucdecayenergy(28, 57) / MEV, nucdecayenergygamma(28, 57) / MEV,
           nucdecayenergy(27, 57) / MEV);
  printout("decayenergy(CR48), decayenergy(V48): %g %g\n",
           nucdecayenergy(24, 48) / MEV, nucdecayenergy(23, 48) / MEV);
  printout("decayenergy(FE52), decayenergy(MN52): %g %g\n",
           nucdecayenergy(26, 52) / MEV, nucdecayenergy(25, 52) / MEV);
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
static double sample_decaytime(const int decaychainindex, const double tdecaymin, const double tdecaymax)
{
  double tdecay = -1;
  while (tdecay <= tdecaymin || tdecay >= tdecaymax)
  {
    tdecay = 0.;

    for (size_t i = 0; i < decaychains_z[decaychainindex].size(); i++)
    {
      const int z = decaychains_z[decaychainindex][i];
      const int a = decaychains_a[decaychainindex][i];
      const double zrand = gsl_rng_uniform_pos(rng);
      tdecay += -get_meanlife(z, a) * log(zrand);
    }
  }
  return tdecay;
}


__host__ __device__
static double calculate_decaychain(
  const double firstinitabund, const double *meanlifetimes, const int num_nuclides, const double timediff, const enum nucmodes mode)
{
  // calculate final abundance or decayrate from multiple decays, e.g., Ni56 -> Co56 -> Fe56 (nuc[0] -> nuc[1] -> nuc[2])
  // the top nuclide initial abundance is set and the chain-end abundance is returned (all intermediates nuclides
  // are assumed to start with zero abundance)
  // note: first and last can be nuclide can be the same if num_nuclides==1, reducing to simple decay formula
  //
  // meanlifetimes: array of mean lifetimes for nuc[0]..nuc[num_nuclides-1]
  // mode:
  //   MODE_ABUND:          assuming intermediate nuclides start with no abundance,
  //                          return the abundance at the end of the chain
  //   MODE_DECAYRATE:      will return the last nuclide's decay rate instead of its abundance
  //   MODE_ABUNDEXPANSION: return a modified abundance at the end of the chain, with a weighting factor2
  //                          accounting from phton energy loss from expansion since the prior decays occured
  //                          (This is needed to get the initial temperature)

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

    if (mode == MODE_ABUND) // abundance output
    {
      sum += exp(-lambdas[j] * timediff) / denominator;
    }
    else if (mode == MODE_DECAYRATE)  // decay rate: derivative of abundance with respect to time
    {
      sum += lambdas[j] * exp(-lambdas[j] * timediff) / denominator;
    }
    else if (mode == MODE_ABUNDEXPANSION)
    {
      if (lambdas[j] > 0.)
      {
        const double sumtermtop = (1 + meanlifetimes[j] / timediff) * exp(-timediff / meanlifetimes[j]) - meanlifetimes[j] / timediff;
        sum += sumtermtop / denominator;
      }
    }
  }

  const double lastabund = firstinitabund * lambdaproduct * sum;

  return lastabund;
}


__host__ __device__
static double get_nuc_abund_or_rate(
  const int modelgridindex, const int z, const int a, const double time, const enum nucmodes mode)
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

  assert_always(mode == MODE_ABUND || mode == MODE_DECAYRATE);
  assert_always(mode != MODE_DECAYRATE || (nuc_exists(z, a) && get_meanlife(z, a))); // decay rate for radioactive nuclides only

  const double t_afterinit = time - get_t_model();

  // decay chains include all paths from radionuclides to other radionuclides (including trivial size-one chains)

  double nuctotal = 0.;  // abundance or decay rate, depending on mode parameter
  for (size_t decaychainindex = 0; decaychainindex < decaychains_a.size(); decaychainindex++)
  {
    const int z_top = decaychains_z[decaychainindex][0];
    const int a_top = decaychains_a[decaychainindex][0];
    const int z_end = decaychains_z[decaychainindex].back();
    const int a_end = decaychains_a[decaychainindex].back();

    if (nuc_exists(z, a) && !(z_end == z && a_end == a)) // requested nuclide is radioactive, so match last nuc in chain
    {
      continue;
    }

    if (!nuc_exists(z, a) && !nuc_is_parent(z_end, a_end, z, a)) // requested nuclide is stable, so add match daughter of last nuc in chain
    {
      continue;
    }

    const double top_initabund = get_modelinitradioabund(modelgridindex, z_top, a_top);
    assert_always(top_initabund >= 0.)
    if (top_initabund <= 0.)
    {
      return 0.;
    }

    int chainlength = decaychains_z[decaychainindex].size();
    double meanlifetimes[chainlength + 1];
    for (int i = 0; i < chainlength; i++)
    {
      meanlifetimes[i] = get_meanlife(decaychains_z[decaychainindex][i], decaychains_a[decaychainindex][i]);
    }

    int fullchainlength = chainlength;
    if (!nuc_exists(z, a))
    {
      // the nuclide is past the end of the chain, in case requested (Z, A) is stable and not in the radionuclides
      meanlifetimes[chainlength] = -1.;
      fullchainlength = chainlength + 1;
    }

    nuctotal += calculate_decaychain(top_initabund, meanlifetimes, fullchainlength, t_afterinit, mode);
   }

  return nuctotal;
}


static double get_nuc_abund(
  const int modelgridindex, const int z, const int a, const double time)
{
  return get_nuc_abund_or_rate(modelgridindex, z, a, time, MODE_ABUND);
}


static double get_nuc_decayrate(
  const int modelgridindex, const int z, const int a, const double time)
{
  return get_nuc_abund_or_rate(modelgridindex, z, a, time, MODE_DECAYRATE);
}


__host__ __device__
static double get_endecay_to_tinf_per_ejectamass_at_time(
  const int modelgridindex, const int decaychainindex, const double time)
// returns decay energy [erg] that would be released from time tstart [s] to time infinity by a given decaypath
{
  // e.g. NI56 -> CO56, represents the decay of Co56 nuclei
  // that were produced from Ni56 in the initial abundance.
  // Decays from Co56 due to the initial abundance of Co56 are not counted here,
  // nor is the energy from Ni56 decays
  // decaying nucleus at the end of the chain

  const int z_top = decaychains_z[decaychainindex][0];
  const int a_top = decaychains_a[decaychainindex][0];
  // if we're a single-nuclide decay chain, then contribution the initial abundance, otherwise contribute
  // all ancestors
  const int z_end = decaychains_z[decaychainindex].back();
  const int a_end = decaychains_a[decaychainindex].back();

  const int chainlength = decaychains_z[decaychainindex].size();
  double meanlifetimes[chainlength + 1];
  for (int i = 0; i < chainlength; i++)
  {
    meanlifetimes[i] = get_meanlife(decaychains_z[decaychainindex][i], decaychains_a[decaychainindex][i]);
  }
  // the nuclide past the end of the chain radionuclide
  meanlifetimes[chainlength] = -1.; // nuclide at the end is a sink, so treat it as stable (even if it's not)

  const double top_initabund = get_modelinitradioabund(modelgridindex, z_top, a_top) / nucmass(z_top, a_top);
  assert_always(top_initabund >= 0.)
  if (top_initabund <= 0.)
  {
    return 0.;
  }
  const double t_afterinit = time - get_t_model();

  // count the number of chain-top nuclei that haven't decayed past the end of the chain

  const double abund_endplusone = calculate_decaychain(top_initabund, meanlifetimes, chainlength + 1, t_afterinit, MODE_ABUND);
  const double ndecays_remaining = top_initabund - abund_endplusone;

  // // alternative: add up the ancestor abundances that will eventually cause decays at the end of chain
  // double ndecays_remaining = 0.;
  // for (int c = 1; c <= chainlength; c++)
  // {
  //   ndecays_remaining += calculate_decaychain(top_initabund, meanlifetimes, c, t_afterinit);
  // }

  const double endecay = ndecays_remaining * nucdecayenergy(z_end, a_end);

  return endecay;
}


__host__ __device__
double get_endecay_per_ejectamass_t0_to_time_withexpansion_chain_numerical(
  const int modelgridindex, const int chainindex, const double tstart)
// just here as as check on the analytic result from get_endecay_per_ejectamass_t0_to_time_withexpansion()
// this version does an Euler integration
{
  double min_meanlife = -1;
  for (size_t i = 0; i < decaychains_a[chainindex].size(); i++)
  {
    const double meanlife = get_meanlife(decaychains_z[chainindex][i], decaychains_a[chainindex][i]);
    if (min_meanlife < 0. or meanlife < min_meanlife)
    {
      min_meanlife = meanlife;
    }
  }

  const int nsteps = ceil((tstart - get_t_model()) / min_meanlife) * 100000; // min steps across the meanlifetime
  double chain_endecay = 0.;
  double last_chain_endecay = -1.;
  double last_t = -1.;
  for (int i = 0; i < nsteps; i++)
  {
    const double t = get_t_model() + (tstart - get_t_model()) * i / nsteps;
    const double chain_endecay_t = get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, chainindex, t);
    if (last_chain_endecay >= 0)
    {
      const double chain_step_endecay_diff = last_chain_endecay - chain_endecay_t;
      const double expansionfactor = 0.5 * (t + last_t) / tstart; // photons lose energy as 1/t for homologous expansion
      chain_endecay += chain_step_endecay_diff * expansionfactor;
    }
    last_chain_endecay = chain_endecay_t;
    last_t = t;
  }

  const double chain_endecay_noexpansion = (
    get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, chainindex, get_t_model()) - get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, chainindex, tstart));

  printout("  chain_endecay:              %g\n", chain_endecay);
  printout("  chain_endecay_noexpansion:  %g\n", chain_endecay_noexpansion);
  printout("  expansion energy factor:    %g\n", chain_endecay / chain_endecay_noexpansion);

  return chain_endecay;
}


__host__ __device__
double get_endecay_per_ejectamass_t0_to_time_withexpansion(const int modelgridindex, const double tstart)
// calculate the decay energy per unit mass [erg/g] released from time zero to tstart, accounting for
// the photon energy loss due to expansion between time of decays and tstart (equation 18 of Lucy 2005)
{
  double tot_endecay = 0.;
  for (size_t chainindex = 0; chainindex < decaychains_z.size(); chainindex++)
  {
    if (get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, chainindex, get_t_model()) <= 0.)
    {
      // skip unused chains
      continue;
    }
    // print_chain(chainindex);
    // get_endecay_per_ejectamass_t0_to_time_withexpansion_chain_numerical(modelgridindex, chainindex, tstart);

    const int chainlength = decaychains_z[chainindex].size();
    double meanlifetimes[chainlength + 1];
    for (int i = 0; i < chainlength; i++)
    {
      meanlifetimes[i] = get_meanlife(decaychains_z[chainindex][i], decaychains_a[chainindex][i]);
    }
    // the nuclide past the end of the chain radionuclide
    meanlifetimes[chainlength] = -1.; // nuclide at the end is a sink, so treat it as stable (even if it's not)

    // const double numerator = calculate_decaychain(1., meanlifetimes, chainlength + 1, tdiff, MODE_ABUNDEXPANSION);
    // const double factor = numerator / calculate_decaychain(1., meanlifetimes, chainlength + 1, tdiff, MODE_ABUND);
    // printout("  Analytical expansion factor: %g\n", factor);

    const int z_top = decaychains_z[chainindex][0];
    const int a_top = decaychains_a[chainindex][0];
    const int z_end = decaychains_z[chainindex].back();
    const int a_end = decaychains_a[chainindex].back();
    const double top_initabund = get_modelinitradioabund(modelgridindex, z_top, a_top) / nucmass(z_top, a_top);
    const double chain_endecay = calculate_decaychain(top_initabund, meanlifetimes, chainlength + 1, tstart - get_t_model(), MODE_ABUNDEXPANSION) * nucdecayenergy(z_end, a_end);
    // printout("  Analytical chain_endecay: %g\n", chain_endecay);
    tot_endecay += chain_endecay;
  }

  return tot_endecay;
}


__host__ __device__
static double get_endecay_per_ejectamass_between_times(
  const int mgi, const int decaychainindex, double tlow, double thigh)
// energy per mass [erg/g] released by a decaypath between two times [s]
{
  assert_always(tlow <= thigh);
  const double energy_tlow = get_endecay_to_tinf_per_ejectamass_at_time(mgi, decaychainindex, tlow);
  const double energy_thigh = get_endecay_to_tinf_per_ejectamass_at_time(mgi, decaychainindex, thigh);
  assert_always(energy_tlow >= energy_thigh);
  return energy_tlow - energy_thigh;
}


__host__ __device__
double get_simtime_endecay_per_ejectamass(const int mgi, const int decaychainindex)
// get the decay energy released during the simulation time
{
#ifdef NO_INITIAL_PACKETS
  // get decay energy released from t=tmin to tmax
  return get_endecay_per_ejectamass_between_times(mgi, decaychainindex, globals::tmin, globals::tmax);
#else
  // get decay energy released from t=0 to tmax
  return get_endecay_per_ejectamass_between_times(mgi, decaychainindex, get_t_model(), globals::tmax);
#endif
}


__host__ __device__
static double get_chain_decay_power_per_ejectamass(
  const int decaychainindex, const int modelgridindex, const double time)
// total decay power per mass [erg / s / kg] for a given decaypath
{
  // only decays at the end of the chain contributed from the initial abundance of the top of the chain are counted
  // (these can be can be same for a chain of length one)

  const int z_top = decaychains_z[decaychainindex][0];
  const int a_top = decaychains_a[decaychainindex][0];
  const int z_end = decaychains_z[decaychainindex].back();
  const int a_end = decaychains_a[decaychainindex].back();

  const double top_initabund = get_modelinitradioabund(modelgridindex, z_top, a_top) / nucmass(z_top, a_top);
  assert_always(top_initabund >= 0.)
  if (top_initabund <= 0.)
  {
    return 0.;
  }

  const int chainlength = decaychains_z[decaychainindex].size();
  double meanlifetimes[chainlength];
  for (int i = 0; i < chainlength; i++)
  {
    meanlifetimes[i] = get_meanlife(decaychains_z[decaychainindex][i], decaychains_a[decaychainindex][i]);
  }

  const double t_afterinit = time - get_t_model();

  const double decaypower = (
    calculate_decaychain(top_initabund, meanlifetimes, chainlength, t_afterinit, MODE_DECAYRATE)
     * nucdecayenergy(z_end, a_end));

   // const double decaypower = (
   //   get_nuc_abund(modelgridindex, z_end, a_end, time) / nucmass(z_top, a_top)
   //   / get_meanlife(z_end, a_end) * nucdecayenergy(z_end, a_end));

  // const double time2 = time * 1.001;
  // const double decaypower2 = get_endecay_per_ejectamass_between_times(modelgridindex, decaychainindex, time, time2) / (time2 - time);
  // printout("compare decaychainindex %d answer %g and %g\n", decaychainindex, decaypower, decaypower2);

  return decaypower;
}


__host__ __device__
double get_modelcell_decay_energy_density(const int mgi)
// get the density at time tmin of decay energy that will
// be released during the simulation time range [erg/cm3]
{
  double modelcell_decay_energy_density = 0.;
  for (size_t decaychainindex = 0; decaychainindex < decaychains_z.size(); decaychainindex++)
  {
    modelcell_decay_energy_density += (
      get_rhoinit(mgi) * get_simtime_endecay_per_ejectamass(mgi, decaychainindex));
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
    const int z = get_nuc_z(nucindex);
    const int a = get_nuc_a(nucindex);
    if (nucdecayenergypositrons(z, a) > 0. && get_nuc_abund(modelgridindex, z, a, t) > 0.)
    {
      // printout("positrons coming from nuclide %d en %g abund %g\n", nuclide, nucdecayenergypositrons(nuclide), get_nuc_abund(modelgridindex, nuclide, t));
      const double decayratefactor = get_nuc_decayrate(modelgridindex, z, a, t);
      pos_dep_sum += decayratefactor * nucdecayenergypositrons(z, a) * rho / nucmass(z, a);
    }
  }

  return pos_dep_sum;
}


double get_global_etot_t0_tinf(void)
{
  double etot_tinf = 0.;
  for (size_t decaychainindex = 0; decaychainindex < decaychains_z.size(); decaychainindex++)
  {
    const int z_top = decaychains_z[decaychainindex][0];
    const int a_top = decaychains_a[decaychainindex][0];
    // if we're a single-nuclide decay chain, then contribution the initial abundance, otherwise contribute
    // all ancestors
    const int z_end = decaychains_z[decaychainindex].back();
    const int a_end = decaychains_a[decaychainindex].back();
    etot_tinf += (
      get_totmassradionuclide(z_top, a_top) / nucmass(z_top, a_top) * nucdecayenergy(z_end, a_end));
  }
  return etot_tinf;
}


__host__ __device__
void update_abundances(const int modelgridindex, const int timestep, const double t_current)
/// Updates the mass fractions of elements associated with the decay sequence
/// Parameters: - modelgridindex: the grid cell for which to update the abundances
///             - t_current: current time (here mid of current timestep)
{
  assert_always(!globals::homogeneous_abundances); // no longer supported

  printout("update_abundances for cell %d timestep %d\n", modelgridindex, timestep);

  double nucfracsum = 0.;
  for (int element = get_nelements() - 1; element >= 0; element--)
  {
    const int atomic_number = get_element(element);
    double isofracsum = 0.; // mass fraction sum of radioactive isotopes, and stable nuclei coming from other decays
    for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
    {
      const int nuc_z = get_nuc_z(nucindex);
      const int a = get_nuc_a(nucindex);
      if (nuc_z == atomic_number)
      {
        // radioactive isotope of the element
        isofracsum += get_nuc_abund(modelgridindex, atomic_number, a, t_current);
      }
      else if (!nuc_exists(decay_daughter_z(nuc_z, a), decay_daughter_a(nuc_z, a)) && decay_daughter_z(nuc_z, a) == atomic_number)
      {
        // nuclide decays into correct atomic number but outside of the radionuclide list
        // note: there could also be stable isotopes of this element included in stable_initabund(z), but
        // here we only count the contribution from decays
        isofracsum += get_nuc_abund(modelgridindex, decay_daughter_z(nuc_z, a), decay_daughter_a(nuc_z, a), t_current);
      }
    }

    const double elmassfrac = get_stable_initabund(modelgridindex, element) + isofracsum;
    set_elem_abundance(modelgridindex, element, elmassfrac);

    nucfracsum += isofracsum;
  }

  double initnucfracsum = 0.;
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
  {
    const int z = get_nuc_z(nucindex);
    if (z < 1) // FAKE_GAM_LINE_ID
      continue;
    initnucfracsum += get_modelinitradioabund(modelgridindex, z, get_nuc_a(nucindex));
  }

  assert_always(fabs(nucfracsum - initnucfracsum) < 0.001); // decays shouldn't change nuclear mass fraction sum

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
  double lower_sum = 0.;
  double cumulative_decay_energy_per_mass[decaychains_z.size()];
  for (size_t decaychainindex = 0; decaychainindex < decaychains_z.size(); decaychainindex++)
  {
    // visit each radioactive nuclide and any chains of ancestors
    // the ancestor chains need to be treated separately so that the decay time can be randomly sampled
    double simtime_endecay_thispath = 0.;

    simtime_endecay_thispath = get_simtime_endecay_per_ejectamass(mgi, decaychainindex);
    cumulative_decay_energy_per_mass[decaychainindex] = lower_sum + simtime_endecay_thispath;
    lower_sum += simtime_endecay_thispath;
  }

  const double zrand_chain = gsl_rng_uniform(rng) * cumulative_decay_energy_per_mass[decaychains_z.size() - 1];

  int decaychainindex = -1;
  for (size_t i = 0; i < decaychains_z.size(); i++)
  {
    if (zrand_chain <= cumulative_decay_energy_per_mass[i])
    {
      decaychainindex = i;
      break;
    }
  }
  assert_always(decaychainindex >= 0); // Failed to select pellet

  #ifdef NO_INITIAL_PACKETS
  const double tdecaymin = globals::tmin;
  #else
  const double tdecaymin = 0.; // allow decays before the first timestep
  #endif

  if (UNIFORM_PELLET_ENERGIES)
  {
    pkt_ptr->tdecay = sample_decaytime(decaychainindex, tdecaymin, globals::tmax);
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
    const double avgpower = get_simtime_endecay_per_ejectamass(mgi, decaychainindex) / (globals::tmax - tdecaymin);
    pkt_ptr->e_cmf = e0 * get_chain_decay_power_per_ejectamass(decaychainindex, mgi, pkt_ptr->tdecay) / avgpower;
    // assert_always(pkt_ptr->e_cmf >= 0);
  }

  // final decaying nuclide at the end of the chain
  const int z = decaychains_z[decaychainindex].back();
  const int a = decaychains_a[decaychainindex].back();

  pkt_ptr->type = TYPE_RADIOACTIVE_PELLET;
  pkt_ptr->pellet_nucindex = get_nuc_index(z, a);

  const double zrand = gsl_rng_uniform(rng);
  pkt_ptr->originated_from_positron = (zrand >= nucdecayenergygamma(z, a) / nucdecayenergy(z, a));
}

}