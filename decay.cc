#include "sn3d.h"
#include "atomic.h"
#include "gamma.h"
#include "grid.h"
#include "input.h"
#include "nonthermal.h"

#include "decay.h"

#include <algorithm> // std::max
#include <vector>
#include <string>
#include <regex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace decay
{

const char *elsymbols[119] = {
  "n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
  "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo"};

enum decaytypes {
  DECAYTYPE_ELECTRONCAPTURE,
  DECAYTYPE_BETAPLUS,
  DECAYTYPE_BETAMINUS,
  DECAYTYPE_NONE,
};

struct nuclide {
  int z;                     // atomic number
  int a;                     // mass number
  double meanlife;           // mean lifetime before decay [s]
  double endecay_positrons;  // average energy per decay in kinetic energy of emitted positrons [erg]
  double endecay_gamma;      // average energy per decay in gamma rays [erg]
  enum decaytypes decaytype;
};

struct nuclide *nuclides = NULL;
int num_nuclides = 0;

struct decaypath {
  int pathlength;
  int *z;                     // atomic number
  int *a;                     // mass number
};

std::vector<struct decaypath> decaychains;

// cumulative_chain_energy_per_mass point to an array of length npts_model * num_decaypaths
// the index [mgi * num_decaypaths + i] will hold the decay energy released by chains [0, i] in cell mgi
double *cumulative_chain_energy_per_mass = NULL;


__host__ __device__
int get_num_nuclides(void)
{
  assert_always(num_nuclides > 0);
  return num_nuclides;
}


static void printout_nuclidename(const int z, const int a)
{
  printout("(Z=%d)%s%d", z, elsymbols[z], a);
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
  const enum decaytypes decaytype = nuclides[get_nuc_index(z_parent, a_parent)].decaytype;
  switch (decaytype)
  {
    case DECAYTYPE_BETAPLUS:
    case DECAYTYPE_ELECTRONCAPTURE:
    {
      return z_parent - 1; // lose a proton, gain a neutron
    }
    case DECAYTYPE_BETAMINUS:
    {
      return z_parent + 1;  // lose a neutron, gain a proton
    }
    case DECAYTYPE_NONE:
    {
      return -1; // no daughter
    }
  }
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


static int get_num_decaypaths(void)
{
  return decaychains.size();
}


static int get_decaypathlength(int decaypathindex)
{
  return decaychains[decaypathindex].pathlength;
}


static void printout_chain(const int decaypathindex)
{
  if (decaychains.size() > 0)
  {
    printout(" decay chain %d: ", decaypathindex);
    printout_nuclidename(decaychains[decaypathindex].z[0], decaychains[decaypathindex].a[0]);
  }

  for (int i = 1; i < get_decaypathlength(decaypathindex); i++)
  {
    printout(" -> ");
    printout_nuclidename(decaychains[decaypathindex].z[i], decaychains[decaypathindex].a[i]);
  }
  printout("\n");
}


static void add_ancestorchains(const int z, const int a, const int startdecaypathindex)
{
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
  {
    const int z_parent = get_nuc_z(nucindex); // possible parent, will check if true
    const int a_parent = get_nuc_a(nucindex);
    // printout("z_parent %d a_parent %d isparent(%d, %d) %d\n", z_parent, a_parent, z_list[0], a_list[0], nuc_is_parent(z_parent, a_parent, z_list[0], a_list[0]));
    if (nuc_is_parent(z_parent, a_parent, z, a))
    {
      struct decaypath newdecaypath;
      newdecaypath.pathlength = get_decaypathlength(startdecaypathindex) + 1;
      newdecaypath.z = (int *) malloc((get_decaypathlength(startdecaypathindex) + 1) * sizeof(int));
      newdecaypath.a = (int *) malloc((get_decaypathlength(startdecaypathindex) + 1) * sizeof(int));

      // check for repeated nuclides, which would indicate a loop in the decay chain
      for (int i = 1; i < newdecaypath.pathlength; i++)
      {
        newdecaypath.z[i] = decaychains[startdecaypathindex].z[i - 1];
        newdecaypath.a[i] = decaychains[startdecaypathindex].a[i - 1];
        if (newdecaypath.z[i] == z_parent && newdecaypath.a[i] == a_parent)
        {
          printout("\nERROR: Loop found in nuclear decay chain.\n");
          abort();
        }
      }
      newdecaypath.z[0] = z_parent;
      newdecaypath.a[0] = a_parent;
      decaychains.push_back(newdecaypath);

      add_ancestorchains(z_parent, a_parent, decaychains.size() - 1);
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

    struct decaypath newdecaypath;
    newdecaypath.pathlength = 1;
    newdecaypath.z = (int *) malloc(sizeof(int));
    newdecaypath.a = (int *) malloc(sizeof(int));

    newdecaypath.z[0] = get_nuc_z(endnuc);
    newdecaypath.a[0] = get_nuc_a(endnuc);
    decaychains.push_back(newdecaypath);

    add_ancestorchains(get_nuc_z(endnuc), get_nuc_a(endnuc), decaychains.size() - 1);
  }
}


int get_nucstring_z(const char *strnuc)
// convert something like Ni56 to integer 28
{
  std::string elcode = strnuc;
  elcode.erase(std::remove_if(elcode.begin(), elcode.end(), &isdigit), elcode.end());

  for (int z = 1; z < 110; z++)
  {
    if (strcmp(elcode.c_str(), elsymbols[z]) == 0)  // first to letters match el symbol
    {
      return z;
    }
  }
  printout("Could not get atomic number of '%s' '%s'\n", strnuc, elcode.c_str());
  return -1;
}


int get_nucstring_a(const char *strnuc)
// convert something like Ni56 to integer 56
{
  std::string strmassnum = std::regex_replace(strnuc, std::regex("[^0-9]*([0-9]+).*"), std::string("$1"));
  const int a = std::stoi(strmassnum);
  assert(a > 0);
  return a;
}


__host__ __device__
void init_nuclides(std::vector<int> custom_zlist, std::vector<int> custom_alist) // std::vector<std::string> nuclides
{
  // all decays are currently assumed to be electron-capture or beta+ (Z -> Z - 1)

  assert_always(custom_zlist.size() == custom_alist.size());

  num_nuclides = 9 + custom_zlist.size();
  nuclides = (struct nuclide *) malloc(num_nuclides * sizeof(struct nuclide));
  assert_always(nuclides != NULL);

  for (int nucindex = 0; nucindex < num_nuclides; nucindex++)
  {
    nuclides[nucindex].z = -1;
    nuclides[nucindex].a = -1;
    nuclides[nucindex].meanlife = -1;
    nuclides[nucindex].endecay_positrons = 0.;
    nuclides[nucindex].endecay_gamma = 0.;
    nuclides[nucindex].decaytype = DECAYTYPE_NONE;
  }

  int nucindex = 0;

  nuclides[nucindex].z = 28; // Ni57
  nuclides[nucindex].a = 57;
  nuclides[nucindex].meanlife = 51.36 * 60;
  nuclides[nucindex].endecay_positrons = 0.354 * MEV * 0.436;
  nuclides[nucindex].decaytype = DECAYTYPE_BETAPLUS;
  nucindex++;

  nuclides[nucindex].z = 28; // Ni56
  nuclides[nucindex].a = 56;
  nuclides[nucindex].meanlife = 8.80 * DAY;
  nuclides[nucindex].decaytype = DECAYTYPE_ELECTRONCAPTURE;
  nucindex++;

  nuclides[nucindex].z = 27; // Co56
  nuclides[nucindex].a = 56;
  nuclides[nucindex].meanlife = 113.7 * DAY;
  nuclides[nucindex].endecay_positrons = 0.63 * MEV * 0.19;
  nuclides[nucindex].decaytype = DECAYTYPE_BETAPLUS;
  nucindex++;

  nuclides[nucindex].z = -1;  // FAKE_GAM_LINE_ID
  nuclides[nucindex].a = -1;
  nuclides[nucindex].decaytype = DECAYTYPE_NONE;
  nucindex++;

  nuclides[nucindex].z = 24; // Cr48
  nuclides[nucindex].a = 48;
  nuclides[nucindex].meanlife = 1.29602 * DAY;
  nuclides[nucindex].decaytype = DECAYTYPE_ELECTRONCAPTURE;
  nucindex++;

  nuclides[nucindex].z = 23; // V48
  nuclides[nucindex].a = 48;
  nuclides[nucindex].meanlife = 23.0442 * DAY;
  nuclides[nucindex].endecay_positrons = 0.290 * 0.499 * MEV;
  nuclides[nucindex].decaytype = DECAYTYPE_BETAPLUS;
  nucindex++;

  nuclides[nucindex].z = 27; // Co57
  nuclides[nucindex].a = 57;
  nuclides[nucindex].meanlife = 392.03 * DAY;
  nuclides[nucindex].decaytype = DECAYTYPE_ELECTRONCAPTURE;
  nucindex++;

  nuclides[nucindex].z = 26; // Fe52
  nuclides[nucindex].a = 52;
  nuclides[nucindex].meanlife = 0.497429 * DAY;
  nuclides[nucindex].decaytype = DECAYTYPE_ELECTRONCAPTURE;
  nucindex++;

  nuclides[nucindex].z = 25; // Mn52
  nuclides[nucindex].a = 52;
  nuclides[nucindex].meanlife = 0.0211395 * DAY;
  nuclides[nucindex].decaytype = DECAYTYPE_ELECTRONCAPTURE;
  nucindex++;

  for (int i = 0; i < (int) custom_zlist.size(); i++)
  {
    nuclides[nucindex].z = custom_zlist[i];
    nuclides[nucindex].a = custom_alist[i];
    nuclides[nucindex].decaytype = DECAYTYPE_BETAMINUS;

    // todo: read Hotokezaka files for beta minus
    // file path 'data/betaminus/' + A + '.txt'
    // columns: A, Z, Q[MeV], Egamma[MeV], Eelec[MeV], Eneutrino[MeV], tau[s]

    nucindex++;
  }


  printout("init_nuclides: num_nuclides %d\n", get_num_nuclides());

  /// Read in data for gamma ray lines and make a list of them in energy order.
  init_gamma_linelist();

  find_chains();

  int maxdecaypathlength = 0;
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++)
  {
    printout_chain(decaypathindex);
    maxdecaypathlength = std::max(maxdecaypathlength, get_decaypathlength(decaypathindex));
  }
  printout("Number of decay paths: %d (max length %d)\n", (int) get_num_decaypaths(), maxdecaypathlength);

  // TODO: generalise this to all included nuclides
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
static double sample_decaytime(const int decaypathindex, const double tdecaymin, const double tdecaymax)
{
  double tdecay = -1;
  const double t_model = get_t_model();
  while (tdecay <= tdecaymin || tdecay >= tdecaymax)
  {
    tdecay = t_model; // can't decay before initial model snapshot time

    for (int i = 0; i < get_decaypathlength(decaypathindex); i++)
    {
      const int z = decaychains[decaypathindex].z[i];
      const int a = decaychains[decaypathindex].a[i];
      const double zrand = gsl_rng_uniform_pos(rng);
      tdecay += -get_meanlife(z, a) * log(zrand);
    }
  }
  return tdecay;
}


__host__ __device__
static double calculate_decaychain(
  const double firstinitabund, const double *meanlifetimes, const int num_nuclides, const double timediff, bool useexpansionfactor)
{
  // calculate final abundance or decayrate from multiple decays, e.g., Ni56 -> Co56 -> Fe56 (nuc[0] -> nuc[1] -> nuc[2])
  // the top nuclide initial abundance is set and the chain-end abundance is returned (all intermediates nuclides
  // are assumed to start with zero abundance)
  // note: first and last can be nuclide can be the same if num_nuclides==1, reducing to simple decay formula
  //
  // meanlifetimes:      array of mean lifetimes for nuc[0]..nuc[num_nuclides-1]
  // useexpansionfactor: if true, return a modified abundance at the end of the chain, with a weighting factor
  //                          accounting for photon energy loss from expansion since the decays occured
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

    if (!useexpansionfactor) // abundance output
    {
      sum += exp(-lambdas[j] * timediff) / denominator;
    }
    else
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
static double get_nuc_abund(
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

  const double t_afterinit = time - get_t_model();

  // decay chains include all paths from radionuclides to other radionuclides (including trivial size-one chains)

  double nuctotal = 0.;  // abundance or decay rate, depending on mode parameter
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++)
  {
    const int z_top = decaychains[decaypathindex].z[0];
    const int a_top = decaychains[decaypathindex].a[0];
    const int z_end = decaychains[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
    const int a_end = decaychains[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];

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
      continue;
    }

    int decaypathlength = get_decaypathlength(decaypathindex);
    double meanlifetimes[decaypathlength + 1];
    for (int i = 0; i < decaypathlength; i++)
    {
      meanlifetimes[i] = get_meanlife(decaychains[decaypathindex].z[i], decaychains[decaypathindex].a[i]);
    }

    int fulldecaypathlength = decaypathlength;
    if (!nuc_exists(z, a))
    {
      // the nuclide is past the end of the chain, in case requested (Z, A) is stable and not in the radionuclides
      meanlifetimes[decaypathlength] = -1.;
      fulldecaypathlength = decaypathlength + 1;
    }

    nuctotal += calculate_decaychain(top_initabund, meanlifetimes, fulldecaypathlength, t_afterinit, false);
   }

  return nuctotal;
}


__host__ __device__
static double get_endecay_to_tinf_per_ejectamass_at_time(
  const int modelgridindex, const int decaypathindex, const double time)
// returns decay energy [erg] that would be released from time tstart [s] to time infinity by a given decaypath
{
  // e.g. NI56 -> CO56, represents the decay of Co56 nuclei
  // that were produced from Ni56 in the initial abundance.
  // Decays from Co56 due to the initial abundance of Co56 are not counted here,
  // nor is the energy from Ni56 decays
  // decaying nucleus at the end of the chain

  assert_always(decaypathindex >= 0);
  assert_always(decaypathindex < get_num_decaypaths());

  const int z_top = decaychains[decaypathindex].z[0];
  const int a_top = decaychains[decaypathindex].a[0];
  // if we're a single-nuclide decay chain, then contribution the initial abundance, otherwise contribute
  // all ancestors
  const int z_end = decaychains[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
  const int a_end = decaychains[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];

  const int decaypathlength = get_decaypathlength(decaypathindex);
  double meanlifetimes[decaypathlength + 1];
  for (int i = 0; i < decaypathlength; i++)
  {
    meanlifetimes[i] = get_meanlife(decaychains[decaypathindex].z[i], decaychains[decaypathindex].a[i]);
  }
  // the nuclide past the end of the chain radionuclide
  meanlifetimes[decaypathlength] = -1.; // nuclide at the end is a sink, so treat it as stable (even if it's not)

  const double top_initabund = get_modelinitradioabund(modelgridindex, z_top, a_top) / nucmass(z_top, a_top);
  assert_always(top_initabund >= 0.)
  if (top_initabund <= 0.)
  {
    return 0.;
  }
  const double t_afterinit = time - get_t_model();

  // count the number of chain-top nuclei that haven't decayed past the end of the chain

  const double abund_endplusone = calculate_decaychain(top_initabund, meanlifetimes, decaypathlength + 1, t_afterinit, false);
  const double ndecays_remaining = top_initabund - abund_endplusone;

  // // alternative: add up the ancestor abundances that will eventually cause decays at the end of chain
  // double ndecays_remaining = 0.;
  // for (int c = 1; c <= decaypathlength; c++)
  // {
  //   ndecays_remaining += calculate_decaychain(top_initabund, meanlifetimes, c, t_afterinit);
  // }

  const double endecay = ndecays_remaining * nucdecayenergy(z_end, a_end);

  return endecay;
}


__host__ __device__
double get_endecay_per_ejectamass_t0_to_time_withexpansion_chain_numerical(
  const int modelgridindex, const int decaypathindex, const double tstart)
// just here as as check on the analytic result from get_endecay_per_ejectamass_t0_to_time_withexpansion()
// this version does an Euler integration
{
  double min_meanlife = -1;
  for (int i = 0; i < get_decaypathlength(decaypathindex); i++)
  {
    const double meanlife = get_meanlife(decaychains[decaypathindex].z[i], decaychains[decaypathindex].a[i]);
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
    const double chain_endecay_t = get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, decaypathindex, t);
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
    get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, decaypathindex, get_t_model()) - get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, decaypathindex, tstart));

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
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++)
  {
    if (get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, decaypathindex, get_t_model()) <= 0.)
    {
      // skip unused chains
      continue;
    }
    // printout_chain(decaypathindex);
    // get_endecay_per_ejectamass_t0_to_time_withexpansion_chain_numerical(modelgridindex, decaypathindex, tstart);

    const int decaypathlength = get_decaypathlength(decaypathindex);
    double meanlifetimes[decaypathlength + 1];
    for (int i = 0; i < decaypathlength; i++)
    {
      meanlifetimes[i] = get_meanlife(decaychains[decaypathindex].z[i], decaychains[decaypathindex].a[i]);
    }
    // the nuclide past the end of the chain radionuclide
    meanlifetimes[decaypathlength] = -1.; // nuclide at the end is a sink, so treat it as stable (even if it's not)

    // const double numerator = calculate_decaychain(1., meanlifetimes, decaypathlength + 1, tdiff, true);
    // const double factor = numerator / calculate_decaychain(1., meanlifetimes, decaypathlength + 1, tdiff, MODE_ABUND);
    // printout("  Analytical expansion factor: %g\n", factor);

    const int z_top = decaychains[decaypathindex].z[0];
    const int a_top = decaychains[decaypathindex].a[0];
    const int z_end = decaychains[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
    const int a_end = decaychains[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];
    const double top_initabund = get_modelinitradioabund(modelgridindex, z_top, a_top) / nucmass(z_top, a_top);
    const double chain_endecay = calculate_decaychain(top_initabund, meanlifetimes, decaypathlength + 1, tstart - get_t_model(), true) * nucdecayenergy(z_end, a_end);
    // printout("  Analytical chain_endecay: %g\n", chain_endecay);
    tot_endecay += chain_endecay;
  }

  return tot_endecay;
}


__host__ __device__
static double get_endecay_per_ejectamass_between_times(
  const int mgi, const int decaypathindex, double tlow, double thigh)
// energy per mass [erg/g] released by a decaypath between two times [s]
{
  assert_always(tlow <= thigh);
  const double energy_tlow = get_endecay_to_tinf_per_ejectamass_at_time(mgi, decaypathindex, tlow);
  const double energy_thigh = get_endecay_to_tinf_per_ejectamass_at_time(mgi, decaypathindex, thigh);
  assert_always(energy_tlow >= energy_thigh);
  const double endiff = energy_tlow - energy_thigh;
  assert_always(std::isfinite(endiff));
  return endiff;
}


__host__ __device__
double get_simtime_endecay_per_ejectamass(const int mgi, const int decaypathindex)
// get the decay energy released during the simulation time
{
#ifdef NO_INITIAL_PACKETS
  // get decay energy released from t=tmin to tmax
  return get_endecay_per_ejectamass_between_times(mgi, decaypathindex, globals::tmin, globals::tmax);
#else
  // get decay energy released from t=0 to tmax
  return get_endecay_per_ejectamass_between_times(mgi, decaypathindex, get_t_model(), globals::tmax);
#endif
}


__host__ __device__
static double get_chain_decay_power_per_ejectamass(
  const int decaypathindex, const int modelgridindex, const double time)
// total decay power per mass [erg / s / kg] for a given decaypath
{
  // only decays at the end of the chain contributed from the initial abundance of the top of the chain are counted
  // (these can be can be same for a chain of length one)

  const int z_top = decaychains[decaypathindex].z[0];
  const int a_top = decaychains[decaypathindex].a[0];
  const int z_end = decaychains[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
  const int a_end = decaychains[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];;

  const double top_initabund = get_modelinitradioabund(modelgridindex, z_top, a_top);
  assert_always(top_initabund >= 0.)
  if (top_initabund <= 0.)
  {
    return 0.;
  }

  const double t_afterinit = time - get_t_model();

  int decaypathlength = get_decaypathlength(decaypathindex);
  double meanlifetimes[decaypathlength];
  for (int i = 0; i < decaypathlength; i++)
  {
    meanlifetimes[i] = get_meanlife(decaychains[decaypathindex].z[i], decaychains[decaypathindex].a[i]);
  }

  // contribution to the end nuclide abundance from the top of chain (could be a length-one chain Z,A_top = Z,A_end
  // so contribution would be from init abundance only)
  const double endnucabund = calculate_decaychain(top_initabund, meanlifetimes, decaypathlength, t_afterinit, false);

  const double decaypower = endnucabund / get_meanlife(z_end, a_end) / nucmass(z_top, a_top);

  assert_always(decaypower >= 0.);
  assert_always(std::isfinite(decaypower));

  return decaypower;
}


__host__ __device__
double get_modelcell_decay_energy_density(const int mgi)
// get the density at time tmin of decay energy that will
// be released during the simulation time range [erg/cm3]
{
  assert_always(cumulative_chain_energy_per_mass != NULL);
  // use the last chain from the cell's cumulative chain energies
  return get_rhoinit(mgi) *  cumulative_chain_energy_per_mass[mgi * get_num_decaypaths() + get_num_decaypaths() - 1];
}


void setup_cumulative_chain_energy_per_mass(void)
{
  assert_always(cumulative_chain_energy_per_mass == NULL); // ensure not allocated yet
  cumulative_chain_energy_per_mass = (double *) malloc(get_npts_model() * get_num_decaypaths() * sizeof(double));
  for (int mgi = 0; mgi < get_npts_model(); mgi++)
  {
    double lower_sum = 0.;
    for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++)
    {
      // visit each radioactive nuclide and any chains of ancestors
      // the ancestor chains need to be treated separately so that the decay time can be randomly sampled
      double simtime_endecay_thispath = 0.;

      simtime_endecay_thispath = get_simtime_endecay_per_ejectamass(mgi, decaypathindex);
      cumulative_chain_energy_per_mass[mgi * get_num_decaypaths() + decaypathindex] = lower_sum + simtime_endecay_thispath;
      lower_sum += simtime_endecay_thispath;
    }
  }
}


void free_cumulative_chain_energy_per_mass(void)
{
  free(cumulative_chain_energy_per_mass);
  cumulative_chain_energy_per_mass = NULL;
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
      const double nucdecayrate = get_nuc_abund(modelgridindex, z, a, t) / get_meanlife(z, a);
      // printout("positrons coming from %s-%d en_e+ %g abund %g nucdecayrate %g\n", get_elname(z), a, nucdecayenergypositrons(z, a), get_nuc_abund(modelgridindex, z, a, t), nucdecayrate);
      assert_always(nucdecayrate >= 0);
      pos_dep_sum += nucdecayrate * nucdecayenergypositrons(z, a) * rho / nucmass(z, a);
    }
  }

  return pos_dep_sum;
}


double get_global_etot_t0_tinf(void)
{
  double etot_tinf = 0.;
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++)
  {
    const int z_top = decaychains[decaypathindex].z[0];
    const int a_top = decaychains[decaypathindex].a[0];
    // if we're a single-nuclide decay chain, then contribution the initial abundance, otherwise contribute
    // all ancestors
    const int z_end = decaychains[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
    const int a_end = decaychains[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];
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
  }

  double initnucfracsum = 0.;
  double nucfracsum = 0.;
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
  {
    const int z = get_nuc_z(nucindex);
    const int a = get_nuc_a(nucindex);
    if (z < 1) // FAKE_GAM_LINE_ID
      continue;
    initnucfracsum += get_modelinitradioabund(modelgridindex, z, a);
    nucfracsum += get_nuc_abund(modelgridindex, z, a, t_current);

    // printout_nuclidename(z, a);
    // printout(" init: %g now: %g\n", get_modelinitradioabund(modelgridindex, z, a), get_nuc_abund(modelgridindex, z, a, t_current));

    if (!nuc_exists(decay_daughter_z(z, a), decay_daughter_a(z, a)))
    {
      // printout_nuclidename(decay_daughter_z(z, a), decay_daughter_a(z, a));
      // printout("(stable) init: 0 now: %g\n", get_nuc_abund(modelgridindex, decay_daughter_z(z, a), decay_daughter_a(z, a), t_current));
      // this decay steps off the nuclide list, so add its daughter abundance to the total
      nucfracsum += get_nuc_abund(modelgridindex, decay_daughter_z(z, a), decay_daughter_a(z, a), t_current);
    }
  }

  // printout("initnucfracsum %g\n", initnucfracsum);
  // printout("nucfracsum %g\n", nucfracsum);

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
  assert_testmodeonly(cumulative_chain_energy_per_mass != NULL);
  const double *cumulative_decay_energy_per_mass_thiscell = &cumulative_chain_energy_per_mass[mgi * get_num_decaypaths()];
  const double zrand_chain = gsl_rng_uniform(rng) * cumulative_decay_energy_per_mass_thiscell[get_num_decaypaths() - 1];

  int decaypathindex = -1;
  for (int i = 0; i < get_num_decaypaths(); i++)
  {
    if (cumulative_decay_energy_per_mass_thiscell[i] > zrand_chain)
    {
      decaypathindex = i;
      break;
    }
  }
  assert_always(decaypathindex >= 0); // Failed to select pellet

  #ifdef NO_INITIAL_PACKETS
  const double tdecaymin = globals::tmin;
  #else
  const double tdecaymin = 0.; // allow decays before the first timestep
  #endif

  if (UNIFORM_PELLET_ENERGIES)
  {
    pkt_ptr->tdecay = sample_decaytime(decaypathindex, tdecaymin, globals::tmax);
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
    const double avgpower = get_simtime_endecay_per_ejectamass(mgi, decaypathindex) / (globals::tmax - tdecaymin);
    assert_always(avgpower > 0.);
    assert_always(std::isfinite(avgpower));
    pkt_ptr->e_cmf = e0 * get_chain_decay_power_per_ejectamass(decaypathindex, mgi, pkt_ptr->tdecay) / avgpower;
    assert_always(pkt_ptr->e_cmf >= 0);
    assert_always(std::isfinite(pkt_ptr->e_cmf));
  }

  // final decaying nuclide at the end of the chain
  const int z = decaychains[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
  const int a = decaychains[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];

  pkt_ptr->type = TYPE_RADIOACTIVE_PELLET;
  pkt_ptr->pellet_nucindex = get_nuc_index(z, a);

  const double zrand = gsl_rng_uniform(rng);
  pkt_ptr->originated_from_positron = (zrand >= nucdecayenergygamma(z, a) / nucdecayenergy(z, a));
}

}