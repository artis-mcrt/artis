#include "decay.h"

#include <algorithm>  // std::max
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "atomic.h"
#include "gammapkt.h"
#include "grid.h"
#include "input.h"
#include "nonthermal.h"
#include "sn3d.h"

namespace decay {

constexpr int Z_MAX = 119;
const char *elsymbols[1 + Z_MAX] = {
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",  "Mg", "Al",  "Si", "P",   "S",
    "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",  "Cu", "Zn",  "Ga", "Ge",  "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",  "Pd", "Ag",  "Cd", "In",  "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu", "Gd",  "Tb", "Dy",  "Ho",
    "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",  "Hg", "Tl",  "Pb", "Bi",  "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",  "Bk", "Cf",  "Es", "Fm",  "Md",
    "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo"};

struct nuclide {
  int z;                                // atomic number
  int a;                                // mass number
  double meanlife;                      // mean lifetime before decay [s]
  double endecay_electron;              // average energy per decay in kinetic energy of emitted electons [erg]
  double endecay_positron;              // average energy per decay in kinetic energy of emitted positrons [erg]
  double endecay_gamma;                 // average energy per decay in gamma rays [erg]
  double endecay_alpha;                 // average energy per decay in kinetic energy of alpha particles [erg]
  double endecay_q[DECAYTYPE_COUNT];    // Q-value for decay (reactant minus product energy) of each decay type
  double branchprobs[DECAYTYPE_COUNT];  // branching probabilities of each decay type
};

std::vector<struct nuclide> nuclides;

// a decay path follows the contribution from an initial nuclear abundance
// to another (daughter of last nuclide in decaypath) via decays
// every different path within the network is considered, e.g. 56Ni -> 56Co -> 56Fe is separate to 56Ni -> 56Co
struct decaypath {
  int pathlength;
  std::unique_ptr<int[]> z = nullptr;  // atomic number
  std::unique_ptr<int[]> a = nullptr;  // mass number
  std::unique_ptr<int[]> decaytypes = nullptr;
};

std::vector<struct decaypath> decaypaths;

// decaypath_energy_per_mass points to an array of length npts_model * num_decaypaths
// the index [mgi * num_decaypaths + i] will hold the decay energy per mass [erg/g] released by chain i in cell mgi
// during the simulation time range
double *decaypath_energy_per_mass = NULL;
#ifdef MPI_ON
MPI_Win win_decaypath_energy_per_mass = MPI_WIN_NULL;
#endif

__host__ __device__ int get_num_nuclides(void) { return nuclides.size(); }

const char *get_elname(const int z) {
  assert_always(z <= Z_MAX);
  return elsymbols[z];
}

__host__ __device__ int get_nuc_z(int nucindex) {
  assert_always(nucindex >= 0);
  assert_always(nucindex < get_num_nuclides());
  return nuclides[nucindex].z;
}

__host__ __device__ int get_nuc_a(int nucindex) {
  assert_always(nucindex >= 0);
  assert_always(nucindex < get_num_nuclides());
  return nuclides[nucindex].a;
}

__host__ __device__ int get_nuc_index(int z, int a)
// get the nuclide array index from the atomic number and mass number
{
  assert_always(get_num_nuclides() > 0);

  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
    if (nuclides[nucindex].z == z && nuclides[nucindex].a == a) {
      return nucindex;
    }
  }
  printout("Could not find nuclide Z=%d A=%d\n", z, a);
  assert_always(false);  // nuclide not found
  return -1;
}

__host__ __device__ bool nuc_exists(int z, int a)
// check if nuclide exists in the simulation
{
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
    if (nuclides[nucindex].z == z && nuclides[nucindex].a == a) {
      return true;
    }
  }
  return false;
}

static void printout_nuclidename(const int z, const int a) { printout("(Z=%d)%s%d", z, get_elname(z), a); }

static void printout_nuclidemeanlife(const int z, const int a) {
  if (nuc_exists(z, a) && get_meanlife(z, a) > 0.) {
    printout("[tau %.1es]", get_meanlife(z, a));
  } else if (nuc_exists(z, a)) {
    printout("[stable,in_net]");
  } else {
    printout("[stable,offnet]");
  }
}

__host__ __device__ static double get_nuc_decaybranchprob(const int z_parent, const int a_parent, int decaytype) {
  assert_always(nuc_exists(z_parent, a_parent));
  assert_always(decaytype >= 0);
  assert_always(decaytype < DECAYTYPE_COUNT);
  return nuclides[get_nuc_index(z_parent, a_parent)].branchprobs[decaytype];
}

constexpr int decay_daughter_z(const int z_parent, const int a_parent, int decaytype)
// check if (z_parent, a_parent) is a parent of (z, a)
{
  assert_always(decaytype >= 0);
  assert_always(decaytype < DECAYTYPE_COUNT);

  switch ((enum decaytypes)decaytype) {
    case DECAYTYPE_ALPHA: {
      return z_parent - 2;  // lose two protons and two neutrons
    }
    case DECAYTYPE_BETAPLUS:
    case DECAYTYPE_ELECTRONCAPTURE: {
      return z_parent - 1;  // lose a proton, gain a neutron
    }
    case DECAYTYPE_BETAMINUS: {
      return z_parent + 1;  // lose a neutron, gain a proton
    }
    case DECAYTYPE_NONE: {
      return -1;  // no daughter
    }
    case DECAYTYPE_COUNT: {
      assert_always(false);
    }
  }
  return -1;  // no daughter
}

constexpr int decay_daughter_a(const int z_parent, const int a_parent, int decaytype)
// check if (z_parent, a_parent) is a parent of (z, a)
{
  switch ((enum decaytypes)decaytype) {
    case DECAYTYPE_ALPHA: {
      return a_parent - 4;  // lose two protons and two neutrons
    }
    case DECAYTYPE_BETAPLUS:
    case DECAYTYPE_ELECTRONCAPTURE:
    case DECAYTYPE_BETAMINUS: {
      return a_parent;  // swap a neutron to proton or vice-versa
    }
    case DECAYTYPE_NONE: {
      return -1;  // no daughter
    }
    case DECAYTYPE_COUNT: {
      assert_always(false);
    }
  }
  return -1;  // no daughter
}

static bool nuc_is_parent(const int z_parent, const int a_parent, const int z, const int a)
// check if (z_parent, a_parent) is a parent of (z, a)
{
  assert_always(nuc_exists(z_parent, a_parent));
  // each radioactive nuclide is limited to one daughter nuclide
  for (int dectypeindex = 0; dectypeindex < DECAYTYPE_COUNT; dectypeindex++) {
    if (decay_daughter_z(z_parent, a_parent, dectypeindex) == z &&
        decay_daughter_a(z_parent, a_parent, dectypeindex) == a &&
        get_nuc_decaybranchprob(z_parent, a_parent, dectypeindex) > 0.) {
      return true;
    }
  }
  return false;
}

__host__ __device__ double nucdecayenergygamma(int z, int a)
// average energy (erg) per decay in the form of gamma rays
{
  return nuclides[get_nuc_index(z, a)].endecay_gamma;
}

__host__ __device__ void set_nucdecayenergygamma(int z, int a, double value)
// average energy per decay in the form of gamma rays [erg]
{
  nuclides[get_nuc_index(z, a)].endecay_gamma = value;
}

double nucdecayenergyparticle(const int z_parent, const int a_parent, const int decaytype)
// decay energy in the form of kinetic energy of electrons, positrons, or alpha particles,
// depending on the relevant decay type (but not including neutrinos)
// important: the branching factor has been applied. so, e.g. energy in positrons is
// the average energy per all decays (including electron captures)
{
  assert_always(nuc_exists(z_parent, a_parent));
  assert_always(decaytype >= 0);
  assert_always(decaytype < DECAYTYPE_COUNT);

  switch (decaytype) {
    case DECAYTYPE_ALPHA: {
      return nuclides[get_nuc_index(z_parent, a_parent)].endecay_alpha;
    }
    case DECAYTYPE_BETAPLUS: {
      return nuclides[get_nuc_index(z_parent, a_parent)].endecay_positron;
    }
    case DECAYTYPE_ELECTRONCAPTURE: {
      return 0.;
    }
    case DECAYTYPE_BETAMINUS: {
      return nuclides[get_nuc_index(z_parent, a_parent)].endecay_electron;
    }
    case DECAYTYPE_NONE: {
      return 0.;
    }
  }
  return 0.;
}

__host__ __device__ static double nucdecayenergytotal(int z, int a)
// average energy (erg) per decay in the form of gammas and particles [erg]
{
  double endecay = 0.;
  endecay += nuclides[get_nuc_index(z, a)].endecay_gamma;
  for (int decaytype = 0; decaytype < DECAYTYPE_COUNT; decaytype++) {
    endecay += nucdecayenergyparticle(z, a, decaytype) * get_nuc_decaybranchprob(z, a, decaytype);
  }

  return endecay;
}

double nucdecayenergy(int z, int a, int decaytype)
// contributed energy release per decay [erg] for decaytype (e.g. DECAYTYPE_BETAPLUS)
// (excludes neutrinos!)
{
  assert_always(nuc_exists(z, a));
  const double endecay = nucdecayenergygamma(z, a) + nucdecayenergyparticle(z, a, decaytype);

  return endecay;
}

static double nucdecayenergyqval(int z, int a, int decaytype) {
  return nuclides[get_nuc_index(z, a)].endecay_q[decaytype];
}

__host__ __device__ double get_meanlife(int z, int a) {
  assert_always(z > 0);
  assert_always(a >= z);
  if (!nuc_exists(z, a)) {
    return -1;
  }
  const int nucindex = get_nuc_index(z, a);
  return nuclides[nucindex].meanlife;
}

__host__ __device__ double nucmass(int z, int a) {
  assert_always(z > 0);
  assert_always(a >= z);

  return a * MH;

  // const int nucindex = get_nuc_index(z, a);
  // return nuclides[nucindex].amass;
}

static int get_num_decaypaths(void) { return decaypaths.size(); }

static int get_decaypathlength(int decaypathindex) { return decaypaths[decaypathindex].pathlength; }

static double get_decaypath_branchproduct(int decaypathindex)
// return the product of all branching factors in the decay path
{
  double branchprod = 1.;
  for (int i = 0; i < get_decaypathlength(decaypathindex); i++) {
    branchprod = branchprod * get_nuc_decaybranchprob(decaypaths[decaypathindex].z[i], decaypaths[decaypathindex].a[i],
                                                      decaypaths[decaypathindex].decaytypes[i]);
  }
  return branchprod;
}

static double get_decaypath_lastnucdecayenergy(const int decaypathindex)
// a decaypath's energy is the decay energy of the last nuclide and decaytype in the chain
{
  const int z_end = decaypaths[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
  const int a_end = decaypaths[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];
  const int decaytype_end = decaypaths[decaypathindex].decaytypes[get_decaypathlength(decaypathindex) - 1];
  return nucdecayenergy(z_end, a_end, decaytype_end);
}

static void printout_decaytype(const int decaytype) {
  switch (decaytype) {
    case DECAYTYPE_ALPHA: {
      printout("alpha");
      break;
    }
    case DECAYTYPE_BETAPLUS: {
      printout("beta+");
      break;
    }
    case DECAYTYPE_ELECTRONCAPTURE: {
      printout("ec");
      break;
    }
    case DECAYTYPE_BETAMINUS: {
      printout("beta-");
      break;
    }
    case DECAYTYPE_NONE: {
      printout("none");
      break;
    }
  }
}

static void printout_decaypath(const int decaypathindex) {
  assert_always(decaypaths.size() > 0);
  printout(" decaypath %d: ", decaypathindex);
  printout_nuclidename(decaypaths[decaypathindex].z[0], decaypaths[decaypathindex].a[0]);
  printout_nuclidemeanlife(decaypaths[decaypathindex].z[0], decaypaths[decaypathindex].a[0]);

  for (int i = 1; i < get_decaypathlength(decaypathindex); i++) {
    printout(" -> ");
    printout_decaytype(decaypaths[decaypathindex].decaytypes[i - 1]);
    printout(" -> ");
    printout_nuclidename(decaypaths[decaypathindex].z[i], decaypaths[decaypathindex].a[i]);
    printout_nuclidemeanlife(decaypaths[decaypathindex].z[i], decaypaths[decaypathindex].a[i]);
  }
  const int last_z = decaypaths[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
  const int last_a = decaypaths[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];
  const int last_decaytype = decaypaths[decaypathindex].decaytypes[get_decaypathlength(decaypathindex) - 1];
  const int end_z = decay_daughter_z(last_z, last_a, last_decaytype);
  const int end_a = decay_daughter_a(last_z, last_a, last_decaytype);
  printout(" -> ");
  printout_decaytype(last_decaytype);
  printout(" -> ");
  printout_nuclidename(end_z, end_a);
  printout_nuclidemeanlife(end_z, end_a);

  printout("\n");
}

static void extend_lastdecaypath(void)
// follow decays at the ends of the current list of decaypaths,
// to get decaypaths from all descendants
{
  const int startdecaypathindex = decaypaths.size() - 1;
  const int last_z = decaypaths[startdecaypathindex].z[get_decaypathlength(startdecaypathindex) - 1];
  const int last_a = decaypaths[startdecaypathindex].a[get_decaypathlength(startdecaypathindex) - 1];
  const int dectypeindex = decaypaths[startdecaypathindex].decaytypes[get_decaypathlength(startdecaypathindex) - 1];

  const int daughter_z = decay_daughter_z(last_z, last_a, dectypeindex);
  const int daughter_a = decay_daughter_a(last_z, last_a, dectypeindex);
  if (nuc_exists(daughter_z, daughter_a)) {
    for (int dectypeindex2 = 0; dectypeindex2 < DECAYTYPE_COUNT; dectypeindex2++) {
      if (get_nuc_decaybranchprob(daughter_z, daughter_a, dectypeindex2) == 0.) {
        continue;
      }
      const int pathlength = get_decaypathlength(startdecaypathindex) + 1;
      decaypaths.push_back({.pathlength = pathlength,
                            .z = std::make_unique<int[]>(pathlength),
                            .a = std::make_unique<int[]>(pathlength),
                            .decaytypes = std::make_unique<int[]>(pathlength)});
      const int lastindex = decaypaths.size() - 1;

      // check for repeated nuclides, which would indicate a loop in the decay chain
      for (int i = 0; i < get_decaypathlength(startdecaypathindex); i++) {
        decaypaths[lastindex].z[i] = decaypaths[startdecaypathindex].z[i];
        decaypaths[lastindex].a[i] = decaypaths[startdecaypathindex].a[i];
        decaypaths[lastindex].decaytypes[i] = decaypaths[startdecaypathindex].decaytypes[i];
        if (decaypaths[lastindex].z[i] == daughter_z && decaypaths[lastindex].a[i] == daughter_a) {
          printout("\nERROR: Loop found in nuclear decay chain.\n");
          abort();
        }
      }
      decaypaths[lastindex].z[decaypaths[lastindex].pathlength - 1] = daughter_z;
      decaypaths[lastindex].a[decaypaths[lastindex].pathlength - 1] = daughter_a;
      decaypaths[lastindex].decaytypes[decaypaths[lastindex].pathlength - 1] = dectypeindex2;

      extend_lastdecaypath();
    }
  }
}

constexpr bool operator<(const struct decaypath &d1, const struct decaypath &d2)
// true if d1 < d2
// order the chains in the same way as when the search moved up from the descendant
// instead of down from the ancestor, for ease of test comparison
// chains are sorted by mass number of first, second, third, etc position in chain
{
  const int smallestpathlength = std::min(d1.pathlength, d2.pathlength);
  bool matchingoverlap = true;
  for (int i = 0; i < smallestpathlength; i++) {
    const int d1pos = d1.pathlength - 1 - i;
    // assert_always(d1pos >= 0);
    // assert_always(d1pos < d1.pathlength);
    const int d2pos = d2.pathlength - 1 - i;
    // assert_always(d2pos >= 0);
    // assert_always(d2pos < d2.pathlength);
    // if (get_nuc_index(d1.z[d1pos], d1.a[d1pos]) < get_nuc_index(d2.z[d2pos], d2.a[d2pos]))
    if (d1.a[d1pos] < d2.a[d2pos]) {
      return true;
    } else if (d1.a[d1pos] == d2.a[d2pos] && d1.z[d1pos] < d2.z[d2pos]) {
      return true;
    }
    if (d1.a[d1pos] != d2.a[d2pos] || d1.z[d1pos] != d2.z[d2pos]) {
      matchingoverlap = false;
    }
  }
  // one is an extension of the other
  if (matchingoverlap && d1.pathlength < d2.pathlength) {
    return true;
  }

  return false;
}

static void find_decaypaths(void) {
  for (int startnucindex = 0; startnucindex < get_num_nuclides(); startnucindex++) {
    if (get_nuc_z(startnucindex) < 1)  // FAKE_GAM_LINE_ID
    {
      continue;
    }
    const int z = get_nuc_z(startnucindex);
    const int a = get_nuc_a(startnucindex);

    for (int dectypeindex = 0; dectypeindex < DECAYTYPE_COUNT; dectypeindex++) {
      if (get_nuc_decaybranchprob(z, a, dectypeindex) == 0. || get_meanlife(z, a) <= 0.) {
        continue;
      }

      constexpr int pathlength = 1;
      decaypaths.push_back({.pathlength = pathlength,
                            .z = std::make_unique<int[]>(pathlength),
                            .a = std::make_unique<int[]>(pathlength),
                            .decaytypes = std::make_unique<int[]>(pathlength)});
      const int lastindex = decaypaths.size() - 1;

      decaypaths[lastindex].z[0] = z;
      decaypaths[lastindex].a[0] = a;
      decaypaths[lastindex].decaytypes[0] = dectypeindex;

      extend_lastdecaypath();  // take this single step chain and find all descendants
    }
  }

  std::sort(decaypaths.begin(), decaypaths.end());
}

int get_nucstring_z(const std::string &strnuc)
// convert something like Ni56 to integer 28
{
  std::string elcode = strnuc;
  elcode.erase(std::remove_if(elcode.begin(), elcode.end(), &isdigit), elcode.end());

  for (int z = 1; z <= Z_MAX; z++) {
    if (strcmp(elcode.c_str(), get_elname(z)) == 0)  // first to letters match el symbol
    {
      return z;
    }
  }
  printout("Could not get atomic number of '%s' '%s'\n", strnuc.c_str(), elcode.c_str());
  assert_always(false);  // could not match to an element
  return -1;
}

int get_nucstring_a(const std::string &strnuc)
// convert something like Ni56 to integer 56
{
  // find first digit character
  size_t i = 0;
  for (; i < strnuc.length(); i++) {
    if (isdigit(strnuc[i])) break;
  }

  // remove the non-digit charts
  const std::string strmassnum = strnuc.substr(i);

  const int a = std::atoi(strmassnum.c_str());
  assert_always(a > 0);
  return a;
}

__host__ __device__ void init_nuclides(std::vector<int> custom_zlist, std::vector<int> custom_alist) {
  assert_always(custom_zlist.size() == custom_alist.size());

  struct nuclide default_nuclide;
  default_nuclide.z = -1;
  default_nuclide.a = -1;
  default_nuclide.meanlife = -1;
  default_nuclide.endecay_electron = 0.;
  default_nuclide.endecay_positron = 0.;
  default_nuclide.endecay_gamma = 0.;
  default_nuclide.endecay_alpha = 0.;
  for (int dectypeindex = 0; dectypeindex < DECAYTYPE_COUNT; dectypeindex++) {
    default_nuclide.branchprobs[dectypeindex] = 0.;
    default_nuclide.endecay_q[dectypeindex] = 0.;
  }

  nuclides.push_back(default_nuclide);
  nuclides.back().z = 28;  // Ni57
  nuclides.back().a = 57;
  nuclides.back().meanlife = 51.36 * 60;
  nuclides.back().endecay_positron = 0.354 * MEV * 0.436;
  nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 1.;
  // nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 0.436;
  // nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1. - 0.436;

  nuclides.push_back(default_nuclide);
  nuclides.back().z = 28;  // Ni56
  nuclides.back().a = 56;
  nuclides.back().meanlife = 8.80 * DAY;
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;

  nuclides.push_back(default_nuclide);
  nuclides.back().z = 27;  // Co56
  nuclides.back().a = 56;
  nuclides.back().meanlife = 113.7 * DAY;
  // old method: move BETAPLUS branching factor into average positron energy
  nuclides.back().endecay_positron = 0.63 * MEV * 0.19;
  nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 1.;
  // should produce the same results if everything is done correctly
  // nuclides.back().endecay_positron = 0.63 * MEV;
  // nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 0.19;
  // nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 0.81;

  nuclides.push_back(default_nuclide);
  nuclides.back().z = -1;  // FAKE_GAM_LINE_ID
  nuclides.back().a = -1;
  nuclides.back().branchprobs[DECAYTYPE_NONE] = 1.;

  nuclides.push_back(default_nuclide);
  nuclides.back().z = 24;  // Cr48
  nuclides.back().a = 48;
  nuclides.back().meanlife = 1.29602 * DAY;
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;

  nuclides.push_back(default_nuclide);
  nuclides.back().z = 23;  // V48
  nuclides.back().a = 48;
  nuclides.back().meanlife = 23.0442 * DAY;
  nuclides.back().endecay_positron = 0.290 * MEV * 0.499;
  nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 1.;

  nuclides.push_back(default_nuclide);
  nuclides.back().z = 27;  // Co57
  nuclides.back().a = 57;
  nuclides.back().meanlife = 392.03 * DAY;
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;

  nuclides.push_back(default_nuclide);
  nuclides.back().z = 26;  // Fe52
  nuclides.back().a = 52;
  nuclides.back().meanlife = 0.497429 * DAY;
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;

  nuclides.push_back(default_nuclide);
  nuclides.back().z = 25;  // Mn52
  nuclides.back().a = 52;
  nuclides.back().meanlife = 0.0211395 * DAY;
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;

  if (custom_alist.size() > 0) {
    std::ifstream fbetaminus("betaminusdecays.txt");
    assert_always(fbetaminus.is_open());
    std::string line;
    while (get_noncommentline(fbetaminus, line)) {
      // energies are average per beta decay
      // columns: # A, Z, Q[MeV], E_gamma[MeV], E_elec[MeV], E_neutrino[MeV], meanlife[s]
      int a = -1;
      int z = -1;
      double q_mev = 0.;
      double e_gamma_mev = 0.;
      double e_elec_mev = 0.;
      double e_neutrino = 0.;
      double tau_sec = 0.;
      std::stringstream(line) >> a >> z >> q_mev >> e_gamma_mev >> e_elec_mev >> e_neutrino >> tau_sec;

      bool keeprow = false;  // keep if the mass number matches one of the input nuclides
      for (int i = 0; i < (int)custom_alist.size(); i++) {
        if (custom_alist[i] == a) {
          keeprow = true;
          break;
        }
      }
      if (keeprow) {
        assert_always(!nuc_exists(z, a));
        nuclides.push_back(default_nuclide);
        nuclides.back().z = z;
        nuclides.back().a = a;
        nuclides.back().meanlife = tau_sec;
        nuclides.back().branchprobs[DECAYTYPE_BETAMINUS] = 1.;
        nuclides.back().endecay_q[DECAYTYPE_BETAMINUS] = q_mev * MEV;
        nuclides.back().endecay_electron = e_elec_mev * MEV;
        nuclides.back().endecay_gamma = e_gamma_mev * MEV;
        // printout("betaminus file: Adding (Z=%d)%s-%d endecay_electron %g endecay_gamma %g tau_s %g\n",
        //          z, get_elname(z), a, e_elec_mev, e_gamma_mev, tau_sec);
        assert_always(e_elec_mev >= 0.);
      }
    }
    fbetaminus.close();

    std::ifstream falpha("alphadecays.txt");
    assert_always(falpha.is_open());
    while (get_noncommentline(falpha, line)) {
      // columns: # A, Z, branch_alpha, branch_beta, halflife[s], Q_total_alphadec[MeV], Q_total_betadec[MeV],
      // E_alpha[MeV], E_gamma[MeV], E_beta[MeV]
      int a = -1;
      int z = -1;
      double branch_alpha = 0.;
      double branch_beta = 0.;
      double halflife = 0.;
      double Q_total_alphadec = 0.;
      double Q_total_betadec = 0.;
      double e_alpha_mev = 0.;
      double e_gamma_mev = 0.;
      double e_beta_mev = 0.;
      std::stringstream(line) >> a >> z >> branch_alpha >> branch_beta >> halflife >> Q_total_alphadec >>
          Q_total_betadec >> e_alpha_mev >> e_gamma_mev >> e_beta_mev;

      bool keeprow = ((branch_alpha > 0. || branch_beta > 0.) && halflife > 0.);
      if (keeprow) {
        const double tau_sec = halflife / log(2);
        int alphanucindex = -1;
        if (nuc_exists(z, a)) {
          alphanucindex = get_nuc_index(z, a);
          // printout("compare z %d a %d e_gamma_mev1 %g e_gamma_mev2 %g\n", z, a, nucdecayenergygamma(z, a) / MEV,
          // e_gamma_mev); printout("compare z %d a %d tau1 %g tau2 %g\n", z, a, get_meanlife(z, a), tau_sec);
          // printout("compare z %d a %d e_beta_mev1 %g e_beta_mev2 %g\n", z, a, nuclides[get_nuc_index(z,
          // a)].endecay_positron / MEV, e_beta_mev);
        } else {
          nuclides.push_back(default_nuclide);
          nuclides.back().z = z;
          nuclides.back().a = a;
          nuclides.back().meanlife = tau_sec;
          nuclides.back().endecay_gamma = e_gamma_mev * MEV;
          alphanucindex = nuclides.size() - 1;
        }
        nuclides[alphanucindex].endecay_alpha = e_alpha_mev * MEV;
        nuclides[alphanucindex].branchprobs[DECAYTYPE_BETAMINUS] = branch_beta;
        nuclides[alphanucindex].endecay_q[DECAYTYPE_BETAMINUS] = Q_total_betadec * MEV;
        nuclides[alphanucindex].branchprobs[DECAYTYPE_ALPHA] = branch_alpha;
        nuclides[alphanucindex].endecay_q[DECAYTYPE_ALPHA] = Q_total_alphadec * MEV;

        // printout("alphadecay file: Adding (Z=%d)%s-%d endecay_alpha %g endecay_gamma %g tau_s %g\n",
        //          z, get_elname(z), a, e_alpha_mev, e_gamma_mev, tau_sec);
      }
    }
    falpha.close();
  }

  for (int i = 0; i < (int)custom_alist.size(); i++) {
    const int z = custom_zlist[i];
    const int a = custom_alist[i];
    if (!nuc_exists(z, a)) {
      // printout("Adding Z %d A %d with no decay data (assuming stable)\n", z, a);
      nuclides.push_back(default_nuclide);
      nuclides.back().z = z;
      nuclides.back().a = a;
      nuclides.back().meanlife = -1;
    }
  }

  // TESTING: REMOVE
  // for (int i = 0; i < get_num_nuclides(); i++)
  // {
  //   for (int dectypeindex = 0; dectypeindex < DECAYTYPE_COUNT; dectypeindex++)
  //   {
  //     nuclides[i].branchprobs[dectypeindex] *= 1e-4;
  //     nuclides[i].endecay_gamma = 0.;
  //   }
  // }

  printout("init_nuclides: num_nuclides %d\n", get_num_nuclides());

  /// Read in data for gamma ray lines and make a list of them in energy order.
  gammapkt::init_gamma_linelist();

  find_decaypaths();

  int maxdecaypathlength = 0;
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
    // printout_decaypath(decaypathindex);
    maxdecaypathlength = std::max(maxdecaypathlength, get_decaypathlength(decaypathindex));
  }
  printout("Number of decay paths: %d (max length %d)\n", (int)get_num_decaypaths(), maxdecaypathlength);

  // TODO: generalise this to all included nuclides
  printout("decayenergy(NI56), decayenergy(CO56), decayenergy_gamma(CO56): %g, %g, %g\n",
           nucdecayenergytotal(28, 56) / MEV, nucdecayenergytotal(27, 56) / MEV, nucdecayenergygamma(27, 56) / MEV);
  printout("decayenergy(NI57), decayenergy_gamma(NI57), nucdecayenergy(CO57): %g, %g, %g\n",
           nucdecayenergytotal(28, 57) / MEV, nucdecayenergygamma(28, 57) / MEV, nucdecayenergytotal(27, 57) / MEV);
  printout("decayenergy(CR48), decayenergy(V48): %g %g\n", nucdecayenergytotal(24, 48) / MEV,
           nucdecayenergytotal(23, 48) / MEV);
  printout("decayenergy(FE52), decayenergy(MN52): %g %g\n", nucdecayenergytotal(26, 52) / MEV,
           nucdecayenergytotal(25, 52) / MEV);
}

__host__ __device__ static double sample_decaytime(const int decaypathindex, const double tdecaymin,
                                                   const double tdecaymax) {
  double tdecay = -1;
  const double t_model = grid::get_t_model();
  // rejection method. draw random times until they are within the correct range.
  while (tdecay <= tdecaymin || tdecay >= tdecaymax) {
    tdecay = t_model;  // can't decay before initial model snapshot time

    for (int i = 0; i < get_decaypathlength(decaypathindex); i++) {
      const int z = decaypaths[decaypathindex].z[i];
      const int a = decaypaths[decaypathindex].a[i];
      const double zrand = gsl_rng_uniform_pos(rng);
      tdecay += -get_meanlife(z, a) * log(zrand);
    }
  }
  return tdecay;
}

__host__ __device__ static double calculate_decaychain(const double firstinitabund,
                                                       std::unique_ptr<double[]> &meanlifetimes, const int num_nuclides,
                                                       const double timediff, bool useexpansionfactor) {
  // calculate final number abundance from multiple decays, e.g., Ni56 -> Co56 -> Fe56 (nuc[0] -> nuc[1] -> nuc[2])
  // the top nuclide initial abundance is set and the chain-end abundance is returned (all intermediates nuclides
  // are assumed to start with zero abundance)
  // note: first and last can be nuclide can be the same if num_nuclides==1, reducing to simple decay formula
  //
  // timediff:           time elapsed since firstinitabund was true [seconds]
  // numnuclides:        number of items in meanlifetimes to use
  // meanlifetimes:      array of mean lifetimes for nuc[0]..nuc[num_nuclides-1]  [seconds]
  // useexpansionfactor: if true, return a modified 'abundance' at the end of the chain, with a weighting factor
  //                          accounting for photon energy loss from expansion since the decays occured
  //                          (This is needed to get the initial temperature)

  assert_always(num_nuclides >= 1);

  // if the meanlife is zero or negative, that indicates a stable nuclide

  // last nuclide might be stable (meanlife <= 0.)
  auto lambdas = std::make_unique<double[]>(num_nuclides);
  for (int i = 0; i < num_nuclides; i++) {
    assert_always(meanlifetimes[i] > 0. || (i == num_nuclides - 1));  // only the last nuclide can be stable
    lambdas[i] = (meanlifetimes[i] > 0.) ? 1. / meanlifetimes[i] : 0.;
  }

  double lambdaproduct = 1.;
  for (int j = 0; j < num_nuclides - 1; j++) {
    lambdaproduct *= lambdas[j];
  }

  double sum = 0;
  for (int j = 0; j < num_nuclides; j++) {
    double denominator = 1.;
    for (int p = 0; p < num_nuclides; p++) {
      if (p != j) {
        denominator *= (lambdas[p] - lambdas[j]);
      }
    }

    if (!useexpansionfactor)  // abundance output
    {
      sum += exp(-lambdas[j] * timediff) / denominator;
    } else {
      if (lambdas[j] > 0.) {
        const double sumtermtop =
            (1 + meanlifetimes[j] / timediff) * exp(-timediff / meanlifetimes[j]) - meanlifetimes[j] / timediff;
        sum += sumtermtop / denominator;
      }
    }
  }

  const double lastabund = firstinitabund * lambdaproduct * sum;

  return lastabund;
}

__host__ __device__ static double get_nuc_massfrac(const int modelgridindex, const int z, const int a,
                                                   const double time)
// Get the mass fraction of a nuclide accounting for all decays including those of its parent and grandparent.
// e.g., Co56 abundance may first increase with time due to Ni56 decays, then decease due to Co56 decay
// Can be called for stable nuclides that are one daughters of the radioactive nuclide list e.g., Fe56
// For stable nuclides, abundance returned only comes from other decays (some could be included in init model elem frac)
{
  if (z < 1)  // skip FAKE_GAM_LINE_ID
  {
    return 0.;
  }
  assert_always(time >= 0.);

  const double t_afterinit = time - grid::get_t_model();
  const bool nuc_exists_z_a = nuc_exists(z, a);

  // decay chains include all paths from radionuclides to other radionuclides (including trivial size-one chains)

  double nuctotal = 0.;  // abundance or decay rate, depending on mode parameter
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
    const int z_top = decaypaths[decaypathindex].z[0];
    const int a_top = decaypaths[decaypathindex].a[0];
    const int z_end = decaypaths[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
    const int a_end = decaypaths[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];

    // match 4He contribution from alpha decay of any nucleus
    if (z != 2 || a != 4 ||
        decaypaths[decaypathindex].decaytypes[get_decaypathlength(decaypathindex) - 1] != DECAYTYPE_ALPHA) {
      if (nuc_exists_z_a && !(z_end == z && a_end == a))  // requested nuclide is in network, so match last nuc in chain
      {
        continue;
      }

      if (!nuc_exists_z_a &&
          !nuc_is_parent(z_end, a_end, z,
                         a))  // requested nuclide not in network, so match daughter of last nucleus in chain
      {
        continue;
      }
    }

    const double top_initabund = grid::get_modelinitradioabund(modelgridindex, z_top, a_top) / nucmass(z_top, a_top);
    assert_always(top_initabund >= 0.) if (top_initabund <= 0.) { continue; }

    int decaypathlength = get_decaypathlength(decaypathindex);
    auto meanlifetimes = std::make_unique<double[]>(decaypathlength + 1);
    for (int i = 0; i < decaypathlength; i++) {
      meanlifetimes[i] = get_meanlife(decaypaths[decaypathindex].z[i], decaypaths[decaypathindex].a[i]);
    }

    int fulldecaypathlength = decaypathlength;
    if (!nuc_exists_z_a ||
        (z == 2 && a == 4 &&
         decaypaths[decaypathindex].decaytypes[get_decaypathlength(decaypathindex) - 1] == DECAYTYPE_ALPHA)) {
      // the nuclide is past the end of the chain, in case requested (Z, A) is stable and not in the radionuclides
      meanlifetimes[decaypathlength] = -1.;
      fulldecaypathlength = decaypathlength + 1;
    }

    const double massfraccontrib =
        (get_decaypath_branchproduct(decaypathindex) *
         calculate_decaychain(top_initabund, meanlifetimes, fulldecaypathlength, t_afterinit, false) * nucmass(z, a));
    // assert_always(massfraccontrib >= 0.);
    nuctotal += massfraccontrib;
  }

  // stable nuclei in the network will not have a size-one decay path associated with them,
  // so we need to contribute the initial abundance as-is (no decay)
  if (nuc_exists_z_a && get_meanlife(z, a) <= 0.) {
    nuctotal += grid::get_modelinitradioabund(modelgridindex, z, a);
  }

  return nuctotal;
}

__host__ __device__ static double get_endecay_to_tinf_per_ejectamass_at_time(const int modelgridindex,
                                                                             const int decaypathindex,
                                                                             const double time)
// returns decay energy [erg/g] that would be released from time tstart [s] to time infinity by a given decaypath
{
  // e.g. Ni56 -> Co56, represents the decay of Co56 nuclei
  // that were produced from Ni56 in the initial abundance.
  // Decays from Co56 due to the initial abundance of Co56 are not counted here,
  // nor is the energy from Ni56 decays

  assert_always(decaypathindex >= 0);
  assert_always(decaypathindex < get_num_decaypaths());

  const int z_top = decaypaths[decaypathindex].z[0];
  const int a_top = decaypaths[decaypathindex].a[0];
  // if it's a single-nuclide decay chain, then contribute the initial abundance, otherwise contribute
  // all ancestors

  const int decaypathlength = get_decaypathlength(decaypathindex);
  auto meanlifetimes = std::make_unique<double[]>(decaypathlength + 1);

  for (int i = 0; i < decaypathlength; i++) {
    meanlifetimes[i] = get_meanlife(decaypaths[decaypathindex].z[i], decaypaths[decaypathindex].a[i]);
  }
  // the nuclide past the end of the chain radionuclide
  meanlifetimes[decaypathlength] = -1.;  // nuclide at the end is a sink, so treat it as stable (even if it's not)

  const double top_initabund = grid::get_modelinitradioabund(modelgridindex, z_top, a_top) / nucmass(z_top, a_top);
  assert_always(top_initabund >= 0.) if (top_initabund <= 0.) { return 0.; }
  const double t_afterinit = time - grid::get_t_model();

  // count the number of chain-top nuclei that haven't decayed past the end of the chain

  const double abund_endplusone =
      calculate_decaychain(top_initabund, meanlifetimes, decaypathlength + 1, t_afterinit, false);
  const double ndecays_remaining = get_decaypath_branchproduct(decaypathindex) * (top_initabund - abund_endplusone);

  // // alternative: add up the ancestor abundances that will eventually cause decays at the end of chain
  // double ndecays_remaining = 0.;
  // for (int c = 1; c <= decaypathlength; c++)
  // {
  //   ndecays_remaining += calculate_decaychain(top_initabund, meanlifetimes, c, t_afterinit);
  // }

  const double endecay = ndecays_remaining * get_decaypath_lastnucdecayenergy(decaypathindex);

  return endecay;
}

__host__ __device__ double get_endecay_per_ejectamass_t0_to_time_withexpansion_chain_numerical(const int modelgridindex,
                                                                                               const int decaypathindex,
                                                                                               const double tstart)
// just here as as check on the analytic result from get_endecay_per_ejectamass_t0_to_time_withexpansion()
// this version does an Euler integration
{
  double min_meanlife = -1;
  for (int i = 0; i < get_decaypathlength(decaypathindex); i++) {
    const double meanlife = get_meanlife(decaypaths[decaypathindex].z[i], decaypaths[decaypathindex].a[i]);
    if (min_meanlife < 0. or meanlife < min_meanlife) {
      min_meanlife = meanlife;
    }
  }

  const int nsteps = ceil((tstart - grid::get_t_model()) / min_meanlife) * 100000;  // min steps across the meanlifetime
  double chain_endecay = 0.;
  double last_chain_endecay = -1.;
  double last_t = -1.;
  for (int i = 0; i < nsteps; i++) {
    const double t = grid::get_t_model() + (tstart - grid::get_t_model()) * i / nsteps;
    const double chain_endecay_t = get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, decaypathindex, t);
    if (last_chain_endecay >= 0) {
      const double chain_step_endecay_diff = last_chain_endecay - chain_endecay_t;
      const double expansionfactor =
          0.5 * (t + last_t) / tstart;  // photons lose energy as 1/t for homologous expansion
      chain_endecay += chain_step_endecay_diff * expansionfactor;
    }
    last_chain_endecay = chain_endecay_t;
    last_t = t;
  }

  const double chain_endecay_noexpansion =
      (get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, decaypathindex, grid::get_t_model()) -
       get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, decaypathindex, tstart));

  printout("  chain_endecay:              %g\n", chain_endecay);
  printout("  chain_endecay_noexpansion:  %g\n", chain_endecay_noexpansion);
  printout("  expansion energy factor:    %g\n", chain_endecay / chain_endecay_noexpansion);

  return chain_endecay;
}

__host__ __device__ double get_endecay_per_ejectamass_t0_to_time_withexpansion(const int modelgridindex,
                                                                               const double tstart)
// calculate the decay energy per unit mass [erg/g] released from time t_model to tstart, accounting for
// the photon energy loss due to expansion between time of decays and tstart (equation 18 of Lucy 2005)
{
  double tot_endecay = 0.;
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
    if (get_endecay_to_tinf_per_ejectamass_at_time(modelgridindex, decaypathindex, grid::get_t_model()) <= 0.) {
      // skip unused chains
      continue;
    }
    // printout_decaypath(decaypathindex);
    // get_endecay_per_ejectamass_t0_to_time_withexpansion_chain_numerical(modelgridindex, decaypathindex, tstart);

    const int decaypathlength = get_decaypathlength(decaypathindex);
    auto meanlifetimes = std::make_unique<double[]>(decaypathlength + 1);
    for (int i = 0; i < decaypathlength; i++) {
      meanlifetimes[i] = get_meanlife(decaypaths[decaypathindex].z[i], decaypaths[decaypathindex].a[i]);
    }
    // the nuclide past the end of the chain radionuclide
    meanlifetimes[decaypathlength] = -1.;  // nuclide at the end is a sink, so treat it as stable (even if it's not)

    // const double numerator = calculate_decaychain(1., meanlifetimes, decaypathlength + 1, tdiff, true);
    // const double factor = numerator / calculate_decaychain(1., meanlifetimes, decaypathlength + 1, tdiff, false);
    // printout("  Analytical expansion factor: %g\n", factor);

    const int z_top = decaypaths[decaypathindex].z[0];
    const int a_top = decaypaths[decaypathindex].a[0];

    const double top_initabund = grid::get_modelinitradioabund(modelgridindex, z_top, a_top) / nucmass(z_top, a_top);

    const double chain_endecay =
        (get_decaypath_branchproduct(decaypathindex) *
         calculate_decaychain(top_initabund, meanlifetimes, decaypathlength + 1, tstart - grid::get_t_model(), true) *
         get_decaypath_lastnucdecayenergy(decaypathindex));

    // printout("  Analytical chain_endecay: %g\n", chain_endecay);
    tot_endecay += chain_endecay;
  }

  return tot_endecay;
}

__host__ __device__ static double get_endecay_per_ejectamass_between_times(const int mgi, const int decaypathindex,
                                                                           double tlow, double thigh)
// get decay energy per mass [erg/g] released by a decaypath between times tlow [s] and thigh [s]
{
  assert_always(tlow <= thigh);
  const double energy_tlow = get_endecay_to_tinf_per_ejectamass_at_time(mgi, decaypathindex, tlow);
  const double energy_thigh = get_endecay_to_tinf_per_ejectamass_at_time(mgi, decaypathindex, thigh);
  assert_always(energy_tlow >= energy_thigh);
  const double endiff = energy_tlow - energy_thigh;
  assert_always(std::isfinite(endiff));
  return endiff;
}

static double calculate_simtime_endecay_per_ejectamass(const int mgi, const int decaypathindex)
// calculate the decay energy released during the simulation time per unit mass [erg/g]
{
  assert_testmodeonly(mgi < grid::get_npts_model());
#ifdef NO_INITIAL_PACKETS
  // get decay energy released from t=tmin to tmax
  const double simtime_endecay =
      get_endecay_per_ejectamass_between_times(mgi, decaypathindex, globals::tmin, globals::tmax);
#else
  // get decay energy released from t=0 to tmax
  const double simtime_endecay =
      get_endecay_per_ejectamass_between_times(mgi, decaypathindex, grid::get_t_model(), globals::tmax);
#endif
  return simtime_endecay;
}

__host__ __device__ static double get_simtime_endecay_per_ejectamass(const int mgi, const int decaypathindex)
// get the decay energy released during the simulation time per unit mass [erg/g]
{
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(mgi);
  const double chainendecay = decaypath_energy_per_mass[nonemptymgi * get_num_decaypaths() + decaypathindex];
  assert_testmodeonly(chainendecay >= 0.);
  assert_testmodeonly(std::isfinite(chainendecay));
  return chainendecay;
}

__host__ __device__ static double get_chain_decay_power_per_ejectamass(const int decaypathindex,
                                                                       const int modelgridindex, const double time)
// total decay power per mass [erg/s/g] for a given decaypath
{
  // only decays at the end of the chain contributed from the initial abundance of the top of the chain are counted
  // (these can be can be same for a chain of length one)

  const int z_top = decaypaths[decaypathindex].z[0];
  const int a_top = decaypaths[decaypathindex].a[0];
  const int z_end = decaypaths[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
  const int a_end = decaypaths[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];
  ;

  const double top_initabund = grid::get_modelinitradioabund(modelgridindex, z_top, a_top);
  assert_always(top_initabund >= 0.) if (top_initabund <= 0.) { return 0.; }

  const double t_afterinit = time - grid::get_t_model();

  int decaypathlength = get_decaypathlength(decaypathindex);
  auto meanlifetimes = std::make_unique<double[]>(decaypathlength);
  for (int i = 0; i < decaypathlength; i++) {
    meanlifetimes[i] = get_meanlife(decaypaths[decaypathindex].z[i], decaypaths[decaypathindex].a[i]);
  }

  // contribution to the end nuclide abundance from the top of chain (could be a length-one chain Z,A_top = Z,A_end
  // so contribution would be from init abundance only)
  const double endnucabund = get_decaypath_branchproduct(decaypathindex) *
                             calculate_decaychain(top_initabund, meanlifetimes, decaypathlength, t_afterinit, false);

  const double endecay = get_decaypath_lastnucdecayenergy(decaypathindex);

  const double decaypower = endecay * endnucabund / get_meanlife(z_end, a_end) / nucmass(z_top, a_top);

  assert_always(decaypower >= 0.);
  assert_always(std::isfinite(decaypower));

  return decaypower;
}

__host__ __device__ double get_modelcell_simtime_endecay_per_mass(const int mgi)
// get the density at time tmin of decay energy that will
// be released during the simulation time range [erg/cm3]
{
  double endecay_per_mass = 0.;
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
    endecay_per_mass += get_simtime_endecay_per_ejectamass(mgi, decaypathindex);
  }
  return endecay_per_mass;
}

void setup_decaypath_energy_per_mass(void) {
  const int nonempty_npts_model = grid::get_nonempty_npts_model();
  printout("[info] mem_usage: allocating %.1f MB for decaypath_energy_per_mass...",
           nonempty_npts_model * get_num_decaypaths() * sizeof(double) / 1024. / 1024.);
#ifdef MPI_ON
  int my_rank_cells = nonempty_npts_model / globals::node_nprocs;
  // rank_in_node 0 gets any remainder
  if (globals::rank_in_node == 0) {
    my_rank_cells += nonempty_npts_model - (my_rank_cells * globals::node_nprocs);
  }
  MPI_Aint size = my_rank_cells * get_num_decaypaths() * sizeof(double);

  int disp_unit = sizeof(double);
  assert_always(MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node,
                                        &decaypath_energy_per_mass, &win_decaypath_energy_per_mass) == MPI_SUCCESS);
  assert_always(MPI_Win_shared_query(win_decaypath_energy_per_mass, 0, &size, &disp_unit, &decaypath_energy_per_mass) ==
                MPI_SUCCESS);
#else
  decaypath_energy_per_mass =
      static_cast<double *>(malloc(nonempty_npts_model * get_num_decaypaths() * sizeof(double)));
#endif
  printout("done.\n");

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  printout("Calculating for decaypath_energy_per_mass for all cells...");
  const int num_decaypaths = get_num_decaypaths();
  for (int nonemptymgi = 0; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
    if (nonemptymgi % globals::node_nprocs == globals::rank_in_node) {
      const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
      for (int decaypathindex = 0; decaypathindex < num_decaypaths; decaypathindex++) {
        decaypath_energy_per_mass[nonemptymgi * num_decaypaths + decaypathindex] =
            calculate_simtime_endecay_per_ejectamass(mgi, decaypathindex);
      }
    }
  }
  printout("done.\n");

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void free_decaypath_energy_per_mass(void) {
#ifdef MPI_ON
  if (win_decaypath_energy_per_mass != MPI_WIN_NULL) {
    MPI_Win_free(&win_decaypath_energy_per_mass);
    win_decaypath_energy_per_mass = MPI_WIN_NULL;
  }
#else
  if (decaypath_energy_per_mass != NULL) {
    free(decaypath_energy_per_mass);
    decaypath_energy_per_mass = NULL;
  }
#endif
}

__host__ __device__ double get_particle_injection_rate(const int modelgridindex, const double t, const int decaytype)
// energy release rate in form of kinetic energy of positrons, electrons, and alpha particles in [erg/s/g]
{
  double dep_sum = 0.;
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
    const int z = get_nuc_z(nucindex);
    if (z < 1) {
      continue;
    }
    const int a = get_nuc_a(nucindex);
    const double meanlife = get_meanlife(z, a);
    if (meanlife < 0.) {
      continue;
    }
    const double en_particles = nucdecayenergyparticle(z, a, decaytype);
    if (en_particles > 0.) {
      const double nucdecayrate =
          get_nuc_massfrac(modelgridindex, z, a, t) / meanlife * get_nuc_decaybranchprob(z, a, decaytype);
      assert_always(nucdecayrate >= 0);
      dep_sum += nucdecayrate * en_particles / nucmass(z, a);
    }
  }

  assert_always(std::isfinite(dep_sum));

  return dep_sum;
}

double get_qdot_modelcell(const int modelgridindex, const double t, const int decaytype)
// energy release rate [erg/s/g] including everything (even neutrinos that are ignored elsewhere)
{
  double qdot = 0.;
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
    const int z = get_nuc_z(nucindex);
    if (z < 1) {
      continue;
    }
    const int a = get_nuc_a(nucindex);
    const double meanlife = get_meanlife(z, a);
    if (meanlife <= 0) {
      continue;
    }
    const double q_decay = nucdecayenergyqval(z, a, decaytype) * get_nuc_decaybranchprob(z, a, decaytype);
    if (q_decay <= 0.) {
      continue;
    }
    const double nucdecayrate = get_nuc_massfrac(modelgridindex, z, a, t) / meanlife;
    assert_always(nucdecayrate >= 0);
    qdot += nucdecayrate * q_decay / nucmass(z, a);
  }

  return qdot;
}

double get_global_etot_t0_tinf(void) {
  double etot_tinf = 0.;
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
    const int z_top = decaypaths[decaypathindex].z[0];
    const int a_top = decaypaths[decaypathindex].a[0];

    etot_tinf += (get_decaypath_branchproduct(decaypathindex) * grid::get_totmassradionuclide(z_top, a_top) /
                  nucmass(z_top, a_top) * get_decaypath_lastnucdecayenergy(decaypathindex));
  }
  return etot_tinf;
}

__host__ __device__ void update_abundances(const int modelgridindex, const int timestep, const double t_current)
/// Updates the mass fractions of elements associated with the decay sequence
/// Parameters: - modelgridindex: the grid cell for which to update the abundances
///             - t_current: current time (here mid of current timestep)
{
  assert_always(!globals::homogeneous_abundances);  // no longer supported

  printout("update_abundances for cell %d timestep %d\n", modelgridindex, timestep);

  for (int element = get_nelements() - 1; element >= 0; element--) {
    const int atomic_number = get_element(element);
    std::set<int> a_isotopes;
    // for the current element,
    // the mass fraction sum of radioactive isotopes, and stable nuclei coming from other decays
    double isomassfracsum = 0.;
    double isomassfrac_on_nucmass_sum = 0.;
    for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
      const int nuc_z = get_nuc_z(nucindex);
      const int a = get_nuc_a(nucindex);
      if (nuc_z == atomic_number) {
        // this nucleus is an isotope of the element
        if (!a_isotopes.count(a)) {
          a_isotopes.insert(a);
          const double nuc_massfrac = get_nuc_massfrac(modelgridindex, atomic_number, a, t_current);
          isomassfracsum += nuc_massfrac;
          isomassfrac_on_nucmass_sum += nuc_massfrac / nucmass(atomic_number, a);
        }
      } else {
        // check if the nucleus decays off the network but into the selected element
        for (int dectypeindex = 0; dectypeindex < DECAYTYPE_COUNT; dectypeindex++) {
          const int daughter_z = decay_daughter_z(nuc_z, a, dectypeindex);
          const int daughter_a = decay_daughter_a(nuc_z, a, dectypeindex);
          if (daughter_z == atomic_number && !nuc_exists(daughter_z, daughter_a) &&
              get_nuc_decaybranchprob(nuc_z, a, dectypeindex) > 0.) {
            if (!a_isotopes.count(daughter_a)) {
              a_isotopes.insert(daughter_a);
              // nuclide decays into correct atomic number but outside of the radionuclide list
              // note: there could also be stable isotopes of this element included in stable_initabund(z), but
              // here we only count the contribution from decays
              const double nuc_massfrac = get_nuc_massfrac(modelgridindex, atomic_number, daughter_a, t_current);
              isomassfracsum += nuc_massfrac;
              isomassfrac_on_nucmass_sum += nuc_massfrac / nucmass(atomic_number, daughter_a);
            }
          }
        }
      }
    }

    if (atomic_number == 2 && !nuc_exists(2, 4) && !a_isotopes.count(4)) {
      // 4He will not be identified as a daughter nucleus of above decays, so add it in
      const double nuc_massfrac = get_nuc_massfrac(modelgridindex, 2, 4, t_current);
      isomassfracsum += nuc_massfrac;
      isomassfrac_on_nucmass_sum += nuc_massfrac / nucmass(2, 4);
    }

    const double stable_init_massfrac = grid::get_stable_initabund(modelgridindex, element);
    isomassfracsum += stable_init_massfrac;
    isomassfrac_on_nucmass_sum += stable_init_massfrac / globals::elements[element].initstablemeannucmass;

    grid::set_elem_abundance(modelgridindex, element, isomassfracsum);
    grid::set_element_meanweight(modelgridindex, element, isomassfracsum / isomassfrac_on_nucmass_sum);
  }

  // consider calling calculate_electron_densities() here

  // double initnucfracsum = 0.;
  // double nucfracsum = 0.;
  // for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++)
  // {
  //   const int z = get_nuc_z(nucindex);
  //   const int a = get_nuc_a(nucindex);
  //   if (z < 1) // FAKE_GAM_LINE_ID
  //     continue;
  //   initnucfracsum += grid::get_modelinitradioabund(modelgridindex, z, a);
  //   nucfracsum += get_nuc_massfrac(modelgridindex, z, a, t_current);
  //
  //   // printout_nuclidename(z, a);
  //   // printout(" init: %g now: %g\n", grid::get_modelinitradioabund(modelgridindex, z, a),
  //   get_nuc_massfrac(modelgridindex, z, a, t_current));
  //
  //   for (int dectypeindex = 0; dectypeindex < DECAYTYPE_COUNT; dectypeindex++)
  //   {
  //     if (!nuc_exists(decay_daughter_z(z, a, dectypeindex), decay_daughter_a(z, a, dectypeindex)) &&
  //         get_nuc_decaybranchprob(z, a, dectypeindex) > 0.)
  //     {
  //       // printout_nuclidename(decay_daughter_z(z, a), decay_daughter_a(z, a));
  //       // printout("(stable) init: 0 now: %g\n", get_nuc_massfrac(modelgridindex, decay_daughter_z(z, a),
  //       decay_daughter_a(z, a), t_current));
  //       // this decay steps off the nuclide list, so add its daughter abundance to the total
  //       nucfracsum += get_nuc_massfrac(modelgridindex, decay_daughter_z(z, a, dectypeindex), decay_daughter_a(z, a,
  //       dectypeindex), t_current);
  //     }
  //   }
  // }

  // printout("initnucfracsum %g\n", initnucfracsum);
  // printout("nucfracsum %g\n", nucfracsum);

  // assert_always(fabs(nucfracsum - initnucfracsum) < 0.001); // decays shouldn't change nuclear mass fraction sum

  nonthermal::calculate_deposition_rate_density(modelgridindex, timestep);
  // printout("model cell %d at t_current %g has frac: Ni %g Co %g Fe %g, stable: Ni %g Co %g Fe %g\n",
  //          modelgridindex, t_current,
  //          grid::get_elem_abundance(modelgridinex, get_elementindex(28)),
  //          grid::get_elem_abundance(modelgridinex, get_elementindex(27)),
  //          grid::get_elem_abundance(modelgridinex, get_elementindex(26)),
  //          get_fnistabel(modelgridindex), get_fcostable(modelgridindex), get_ffestable(modelgridindex));
}

void fprint_nuc_abundances(FILE *estimators_file, const int modelgridindex, const double t_current, const int element) {
  const double rho = grid::get_rho(modelgridindex);

  const int atomic_number = get_element(element);
  std::set<int> a_isotopes;
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
    const int nuc_z = get_nuc_z(nucindex);
    const int a = get_nuc_a(nucindex);
    if (nuc_z == atomic_number) {
      if (!a_isotopes.count(a)) {
        a_isotopes.insert(a);
        // radioactive isotope of the element
        const double massfrac = get_nuc_massfrac(modelgridindex, atomic_number, a, t_current);
        if (massfrac > 0) {
          const double numberdens = massfrac / nucmass(atomic_number, a) * rho;

          fprintf(estimators_file, "  %s%d: %9.3e", get_elname(atomic_number), a, numberdens);
        }
      }
    } else {
      for (int dectypeindex = 0; dectypeindex < DECAYTYPE_COUNT; dectypeindex++) {
        const int daughter_z = decay_daughter_z(nuc_z, a, dectypeindex);
        const int daughter_a = decay_daughter_a(nuc_z, a, dectypeindex);
        if (!nuc_exists(daughter_z, daughter_a) && daughter_z == atomic_number &&
            get_nuc_decaybranchprob(nuc_z, a, dectypeindex) > 0.) {
          if (!a_isotopes.count(a)) {
            a_isotopes.insert(a);
            // nuclide decays into correct atomic number but outside of the radionuclide list. Daughter is assumed
            // stable
            const double massfrac = get_nuc_massfrac(modelgridindex, atomic_number, a, t_current);
            const double numberdens = massfrac / nucmass(nuc_z, a) * rho;
            fprintf(estimators_file, "  %s%d: %9.3e", get_elname(atomic_number), a, numberdens);
          }
        }
      }
    }
  }

  // factor to convert convert mass fraction to number density
  const double otherstablemassfrac = grid::get_stable_initabund(modelgridindex, element);
  if (otherstablemassfrac > 0) {
    const double meannucmass = globals::elements[element].initstablemeannucmass;
    const double otherstable_numberdens = otherstablemassfrac / meannucmass * grid::get_rho(modelgridindex);
    fprintf(estimators_file, "  %s_otherstable: %9.3e", get_elname(atomic_number), otherstable_numberdens);
  }
  fprintf(estimators_file, "\n");
}

void setup_radioactive_pellet(const double e0, const int mgi, struct packet *pkt_ptr) {
  const int num_decaypaths = get_num_decaypaths();
  auto cumulative_endecay = std::make_unique<double[]>(num_decaypaths);
  double endecaysum = 0.;
  for (int decaypathindex = 0; decaypathindex < num_decaypaths; decaypathindex++) {
    endecaysum += get_simtime_endecay_per_ejectamass(mgi, decaypathindex);
    cumulative_endecay[decaypathindex] = endecaysum;
  }
  assert_testmodeonly(cumulative_endecay[num_decaypaths - 1] > 0.);

  double total_endecay_per_ejectamass = cumulative_endecay[num_decaypaths - 1];
#ifndef NO_INITIAL_PACKETS
  if (USE_MODEL_INITIAL_ENERGY) {
    total_endecay_per_ejectamass += grid::get_initenergyq(mgi);
  }
#endif

  const double zrand_chain = gsl_rng_uniform(rng) * total_endecay_per_ejectamass;

  if (zrand_chain >= cumulative_endecay[num_decaypaths - 1]) {
    assert_always(USE_MODEL_INITIAL_ENERGY);

    pkt_ptr->prop_time = globals::tmin;
    pkt_ptr->tdecay = globals::tmin;
    // pkt_ptr->type = TYPE_PRE_KPKT;
    pkt_ptr->type = TYPE_RADIOACTIVE_PELLET;
    pkt_ptr->e_cmf = e0;
    pkt_ptr->nu_cmf = e0 / H;
    pkt_ptr->pellet_nucindex = -1;
    pkt_ptr->pellet_decaytype = -1;
    return;
  }

  int decaypathindex = -1;

  for (int i = 0; i < num_decaypaths; i++) {
    if (cumulative_endecay[i] > zrand_chain) {
      decaypathindex = i;
      break;
    }
  }

  assert_always(decaypathindex >= 0);  // Failed to select chain

#ifdef NO_INITIAL_PACKETS
  const double tdecaymin = globals::tmin;
#else
  const double tdecaymin = grid::get_t_model();  // allow decays before the first timestep
#endif

  if (UNIFORM_PELLET_ENERGIES) {
    pkt_ptr->tdecay = sample_decaytime(decaypathindex, tdecaymin, globals::tmax);
    pkt_ptr->e_cmf = e0;
  } else {
#ifndef NO_INITIAL_PACKETS
    assert_always(!USE_MODEL_INITIAL_ENERGY);  // not implemented
#endif
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
  const int z = decaypaths[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
  const int a = decaypaths[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];
  const int decaytype = decaypaths[decaypathindex].decaytypes[get_decaypathlength(decaypathindex) - 1];

  pkt_ptr->type = TYPE_RADIOACTIVE_PELLET;
  pkt_ptr->pellet_nucindex = get_nuc_index(z, a);
  pkt_ptr->pellet_decaytype = decaytype;

  const double zrand = gsl_rng_uniform(rng);
  pkt_ptr->originated_from_particlenotgamma =
      (zrand >= nucdecayenergygamma(z, a) / (nucdecayenergygamma(z, a) + nucdecayenergyparticle(z, a, decaytype)));
  pkt_ptr->nu_cmf = nucdecayenergyparticle(z, a, decaytype) / H;
}

void cleanup(void) { free_decaypath_energy_per_mass(); }

}  // namespace decay