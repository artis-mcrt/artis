#include "decay.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "atomic.h"
#include "gammapkt.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "nonthermal.h"
#include "sn3d.h"

namespace decay {

constexpr std::string_view elsymbols[] = {
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",  "Mg", "Al",  "Si", "P",   "S",
    "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",  "Cu", "Zn",  "Ga", "Ge",  "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",  "Pd", "Ag",  "Cd", "In",  "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu", "Gd",  "Tb", "Dy",  "Ho",
    "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",  "Hg", "Tl",  "Pb", "Bi",  "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",  "Bk", "Cf",  "Es", "Fm",  "Md",
    "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo"};
constexpr int Z_MAX = 119;

struct nuclide {
  int z = -1;                    // atomic number
  int a = -1;                    // mass number
  double meanlife = -1;          // mean lifetime before decay [s]
  double endecay_electron = 0.;  // average energy per beta- decay in kinetic energy of emitted electons [erg]
  double endecay_positron = 0.;  // average energy per beta+ decay in kinetic energy of emitted positrons [erg]
  double endecay_gamma = 0.;     // average energy per decay in gamma rays [erg]
  double endecay_alpha = 0.;     // average energy per alpha decay in kinetic energy of alpha particles [erg]
  std::array<double, decaytypes::DECAYTYPE_COUNT> endecay_q = {
      0.};  // Q-value (reactant minus product energy) for each decay type
  std::array<double, decaytypes::DECAYTYPE_COUNT> branchprobs = {0.};  // branch probability of each decay type
};

std::vector<struct nuclide> nuclides;

constexpr auto decay_daughter_z(const int z_parent, const int /*a_parent*/, const int decaytype) -> int
// check if (z_parent, a_parent) is a parent of (z, a)
{
  assert_always(decaytype >= 0);
  assert_always(decaytype < decaytypes::DECAYTYPE_COUNT);

  switch (static_cast<enum decaytypes>(decaytype)) {
    case decaytypes::DECAYTYPE_ALPHA: {
      return z_parent - 2;  // lose two protons and two neutrons
    }
    case decaytypes::DECAYTYPE_BETAPLUS:
    case decaytypes::DECAYTYPE_ELECTRONCAPTURE: {
      return z_parent - 1;  // lose a proton, gain a neutron
    }
    case decaytypes::DECAYTYPE_BETAMINUS: {
      return z_parent + 1;  // lose a neutron, gain a proton
    }
    case decaytypes::DECAYTYPE_NONE: {
      return -1;  // no daughter
    }
    case decaytypes::DECAYTYPE_COUNT: {
      assert_always(false);
    }
  }
  return -1;  // no daughter
}

constexpr auto decay_daughter_a(const int /*z_parent*/, const int a_parent, const int decaytype) -> int
// check if (z_parent, a_parent) is a parent of (z, a)
{
  switch (static_cast<enum decaytypes>(decaytype)) {
    case decaytypes::DECAYTYPE_ALPHA: {
      return a_parent - 4;  // lose two protons and two neutrons
    }
    case decaytypes::DECAYTYPE_BETAPLUS:
    case decaytypes::DECAYTYPE_ELECTRONCAPTURE:
    case decaytypes::DECAYTYPE_BETAMINUS: {
      return a_parent;  // swap a neutron to proton or vice-versa
    }
    case decaytypes::DECAYTYPE_NONE: {
      return -1;  // no daughter
    }
    case decaytypes::DECAYTYPE_COUNT: {
      assert_always(false);
    }
  }
  return -1;  // no daughter
}

// a decay path follows the contribution from an initial nuclear abundance
// to another (daughter of last nuclide in decaypath) via decays
// every different path within the network is considered, e.g. 56Ni -> 56Co -> 56Fe is separate to 56Ni -> 56Co
struct decaypath {
  std::vector<int> z{};         // atomic number
  std::vector<int> a{};         // mass number
  std::vector<int> nucindex{};  // index into nuclides list
  std::vector<int> decaytypes{};
  std::vector<double> lambdas{};
  double branchproduct{};  // product of all branching factors along the path set by calculate_decaypath_branchproduct()

  [[nodiscard]] auto final_daughter_a() const -> int { return decay_daughter_a(z.back(), a.back(), decaytypes.back()); }
  [[nodiscard]] auto final_daughter_z() const -> int { return decay_daughter_z(z.back(), a.back(), decaytypes.back()); }
};

std::vector<struct decaypath> decaypaths;

// decaypath_energy_per_mass points to an array of length npts_model * num_decaypaths
// the index [mgi * num_decaypaths + i] will hold the decay energy per mass [erg/g] released by chain i in cell mgi
// during the simulation time range
double *decaypath_energy_per_mass = nullptr;
#ifdef MPI_ON
MPI_Win win_decaypath_energy_per_mass = MPI_WIN_NULL;
#endif

auto get_num_nuclides() -> int { return nuclides.size(); }

auto get_elname(const int z) -> const char * {
  assert_testmodeonly(z <= Z_MAX);
  return &elsymbols[z].front();
}

auto get_nuc_z(int nucindex) -> int {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < get_num_nuclides());
  return nuclides[nucindex].z;
}

auto get_nuc_a(int nucindex) -> int {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < get_num_nuclides());
  return nuclides[nucindex].a;
}

static auto get_nucindex_or_neg_one(int z, int a) -> int
// get the nuclide array index from the atomic number and mass number
{
  assert_testmodeonly(get_num_nuclides() > 0);
  const int num_nuclides = get_num_nuclides();

  for (int nucindex = 0; nucindex < num_nuclides; nucindex++) {
    if (nuclides[nucindex].z == z && nuclides[nucindex].a == a) {
      return nucindex;
    }
  }
  return -1;  // nuclide not found
}

auto get_nucindex(int z, int a) -> int
// get the nuclide array index from the atomic number and mass number
{
  const int nucindex = get_nucindex_or_neg_one(z, a);
  if (nucindex >= 0) {
    return nucindex;
  }
  printout("Could not find nuclide Z=%d A=%d\n", z, a);
  assert_always(false);  // nuclide not found
  return -1;
}

auto nuc_exists(int z, int a) -> bool
// check if nuclide exists in the simulation
{
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
    if (nuclides[nucindex].z == z && nuclides[nucindex].a == a) {
      return true;
    }
  }
  return false;
}

static auto get_meanlife(int nucindex) -> double {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < get_num_nuclides());
  return nuclides[nucindex].meanlife;
}

static void printout_nuclidename(const int z, const int a) { printout("(Z=%d)%s%d", z, get_elname(z), a); }

static void printout_nuclidemeanlife(const int z, const int a) {
  const int nucindex = get_nucindex_or_neg_one(z, a);
  const bool exists = (nucindex >= 0);
  if (exists && get_meanlife(nucindex) > 0.) {
    printout("[tau %.1es]", get_meanlife(nucindex));
  } else if (exists) {
    printout("[stable,in_net]");
  } else {
    printout("[stable,offnet]");
  }
}

static auto get_nuc_decaybranchprob(const int nucindex, const int decaytype) -> double {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < get_num_nuclides());
  assert_testmodeonly(decaytype >= 0);
  assert_testmodeonly(decaytype < decaytypes::DECAYTYPE_COUNT);
  return nuclides[nucindex].branchprobs[decaytype];
}

static auto get_nuc_decaybranchprob(const int z_parent, const int a_parent, const int decaytype) -> double {
  return get_nuc_decaybranchprob(get_nucindex(z_parent, a_parent), decaytype);
}

static auto nuc_is_parent(const int z_parent, const int a_parent, const int z, const int a) -> bool
// check if (z_parent, a_parent) is a parent of (z, a)
{
  assert_testmodeonly(nuc_exists(z_parent, a_parent));
  // each radioactive nuclide is limited to one daughter nuclide
  return std::any_of(all_decaytypes.begin(), all_decaytypes.end(), [=](const auto decaytype) {
    return decay_daughter_z(z_parent, a_parent, decaytype) == z && decay_daughter_a(z_parent, a_parent, decaytype) == a;
  });
}

auto nucdecayenergygamma(int nucindex) -> double
// average energy per decay in the form of gamma rays [erg]
{
  return nuclides[nucindex].endecay_gamma;
}
auto nucdecayenergygamma(int z, int a) -> double { return nucdecayenergygamma(get_nucindex(z, a)); }

void set_nucdecayenergygamma(int nucindex, double value)
// set average energy per decay in the form of gamma rays [erg]
{
  nuclides[nucindex].endecay_gamma = value;
}

static auto nucdecayenergyparticle(const int nucindex, const int decaytype) -> double
// decay energy in the form of kinetic energy of electrons, positrons, or alpha particles,
// depending on the relevant decay type (but not including neutrinos)
{
  assert_testmodeonly(decaytype >= 0);
  assert_testmodeonly(decaytype < decaytypes::DECAYTYPE_COUNT);

  switch (decaytype) {
    case decaytypes::DECAYTYPE_ALPHA: {
      return nuclides[nucindex].endecay_alpha;
    }
    case decaytypes::DECAYTYPE_BETAPLUS: {
      return nuclides[nucindex].endecay_positron;
    }
    case decaytypes::DECAYTYPE_ELECTRONCAPTURE: {
      return 0.;
    }
    case decaytypes::DECAYTYPE_BETAMINUS: {
      return nuclides[nucindex].endecay_electron;
    }
    default: {
      return 0.;
    }
  }
}

static auto nucdecayenergytotal(const int z, const int a) -> double
// average energy (erg) per decay in the form of gammas and particles [erg]
{
  const int nucindex = get_nucindex(z, a);
  double endecay = 0.;
  endecay += nuclides[nucindex].endecay_gamma;
  for (const auto decaytype : all_decaytypes) {
    endecay += nucdecayenergyparticle(nucindex, decaytype) * get_nuc_decaybranchprob(nucindex, decaytype);
  }

  return endecay;
}

static auto nucdecayenergy(int nucindex, int decaytype) -> double
// contributed energy release per decay [erg] for decaytype (e.g. decaytypes::DECAYTYPE_BETAPLUS)
// (excludes neutrinos!)
{
  const double endecay = nuclides[nucindex].endecay_gamma + nucdecayenergyparticle(nucindex, decaytype);

  return endecay;
}

static auto nucdecayenergyqval(int nucindex, int decaytype) -> double {
  return nuclides[nucindex].endecay_q[decaytype];
}

auto nucmass(int z, int a) -> double {
  assert_testmodeonly(z > 0);
  assert_testmodeonly(a >= z);

  return a * MH;

  // return nuclides[nucindex].amass;
}

static auto get_num_decaypaths() -> int { return static_cast<int>(decaypaths.size()); }

static auto get_decaypathlength(const decaypath &dpath) -> int { return static_cast<int>(dpath.z.size()); }
static auto get_decaypathlength(int decaypathindex) -> int { return get_decaypathlength(decaypaths[decaypathindex]); }

static auto calculate_decaypath_branchproduct(int decaypathindex) -> double
// return the product of all branching factors in the decay path
{
  double branchprod = 1.;
  for (int i = 0; i < get_decaypathlength(decaypathindex); i++) {
    branchprod = branchprod * get_nuc_decaybranchprob(decaypaths[decaypathindex].nucindex[i],
                                                      decaypaths[decaypathindex].decaytypes[i]);
  }
  return branchprod;
}

static auto get_decaypath_lastnucdecayenergy(const int decaypathindex) -> double
// a decaypath's energy is the decay energy of the last nuclide and decaytype in the chain
{
  const int nucindex_end = decaypaths[decaypathindex].nucindex[get_decaypathlength(decaypathindex) - 1];
  const int decaytype_end = decaypaths[decaypathindex].decaytypes[get_decaypathlength(decaypathindex) - 1];
  return nucdecayenergy(nucindex_end, decaytype_end);
}

static void printout_decaytype(const int decaytype) {
  switch (decaytype) {
    case decaytypes::DECAYTYPE_ALPHA: {
      printout("alpha");
      break;
    }
    case decaytypes::DECAYTYPE_BETAPLUS: {
      printout("beta+");
      break;
    }
    case decaytypes::DECAYTYPE_ELECTRONCAPTURE: {
      printout("ec");
      break;
    }
    case decaytypes::DECAYTYPE_BETAMINUS: {
      printout("beta-");
      break;
    }
    case decaytypes::DECAYTYPE_NONE: {
      printout("none");
      break;
    }
    default:
      break;
  }
}

static void printout_decaypath(const int decaypathindex) {
  assert_always(!decaypaths.empty());
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
  const int last_decaytype = decaypaths[decaypathindex].decaytypes[get_decaypathlength(decaypathindex) - 1];
  const int end_z = decaypaths[decaypathindex].final_daughter_z();
  const int end_a = decaypaths[decaypathindex].final_daughter_a();
  printout(" -> ");
  printout_decaytype(last_decaytype);
  printout(" -> ");
  printout_nuclidename(end_z, end_a);
  printout_nuclidemeanlife(end_z, end_a);

  printout("\n");
}

static void extend_lastdecaypath()
// follow decays at the ends of the current list of decaypaths,
// to get decaypaths from all descendants
{
  const int startdecaypathindex = decaypaths.size() - 1;

  const int daughter_z = decaypaths[startdecaypathindex].final_daughter_z();
  const int daughter_a = decaypaths[startdecaypathindex].final_daughter_a();
  if (nuc_exists(daughter_z, daughter_a)) {
    const int daughter_nucindex = get_nucindex(daughter_z, daughter_a);
    for (enum decaytypes dectypeindex2 : all_decaytypes) {
      if (get_nuc_decaybranchprob(daughter_nucindex, dectypeindex2) == 0.) {
        continue;
      }

      // check for nuclide in existing path, which would indicate a loop
      for (int i = 0; i < get_decaypathlength(startdecaypathindex); i++) {
        if (decaypaths[startdecaypathindex].z[i] == daughter_z && decaypaths[startdecaypathindex].a[i] == daughter_a) {
          printout("\nERROR: Loop found in nuclear decay chain.\n");
          abort();
        }
      }

      decaypaths.push_back(decaypaths[startdecaypathindex]);
      const int lastindex = decaypaths.size() - 1;

      decaypaths[lastindex].z.push_back(daughter_z);
      decaypaths[lastindex].a.push_back(daughter_a);
      decaypaths[lastindex].nucindex.push_back(daughter_nucindex);
      decaypaths[lastindex].decaytypes.push_back(dectypeindex2);

      extend_lastdecaypath();
    }
  }
}

static auto operator<(const struct decaypath &d1, const struct decaypath &d2) -> bool
// true if d1 < d2
// chains are sorted by mass number, then atomic number, then length
{
  const int d1_length = get_decaypathlength(d1);
  const int d2_length = get_decaypathlength(d2);
  const int smallestpathlength = std::min(d1_length, d2_length);
  for (int i = 0; i < smallestpathlength; i++) {
    if (d1.a[i] < d2.a[i]) {
      return true;
    }
    if (d1.a[i] > d2.a[i]) {
      return false;
    }
    if (d1.z[i] < d2.z[i]) {
      return true;
    }
    if (d1.z[i] > d2.z[i]) {
      return false;
    }
  }
  // one is an extension of the other
  return d1_length < d2_length;
}

static void find_decaypaths(const std::vector<int> &custom_zlist, const std::vector<int> &custom_alist,
                            std::vector<struct nuclide> &standard_nuclides) {
  decaypaths.clear();
  for (int startnucindex = 0; startnucindex < get_num_nuclides(); startnucindex++) {
    const int z = get_nuc_z(startnucindex);
    const int a = get_nuc_a(startnucindex);

    for (const auto decaytype : all_decaytypes) {
      if (get_nuc_decaybranchprob(startnucindex, decaytype) == 0. || get_meanlife(startnucindex) <= 0.) {
        continue;
      }
      bool is_custom_nuclide = false;
      for (size_t i = 0; i < custom_zlist.size(); i++) {
        if ((z == custom_zlist[i]) && (a == custom_alist[i])) {
          is_custom_nuclide = true;
          break;
        }
      }
      // skip path if it doesn't start from a nuclide in the custom or standard input lists
      if (!is_custom_nuclide &&
          !std::any_of(standard_nuclides.begin(), standard_nuclides.end(),
                       [z, a](const auto &stdnuc) { return (z == stdnuc.z) && (a == stdnuc.a); })) {
        continue;
      }

      decaypaths.push_back({.z = std::vector<int>(1, z),
                            .a = std::vector<int>(1, a),
                            .nucindex = std::vector<int>(1, startnucindex),
                            .decaytypes = std::vector<int>(1, decaytype)});

      extend_lastdecaypath();  // take this single step chain and find all descendants
    }
  }

  std::sort(decaypaths.begin(), decaypaths.end());

  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
    int const decaypathlength = get_decaypathlength(decaypathindex);

    for (int i = 0; i < decaypathlength; i++) {
      const double meanlife = get_meanlife(decaypaths[decaypathindex].nucindex[i]);

      // last nuclide might be stable (meanlife <= 0.)
      assert_always(meanlife > 0. || (i == decaypathlength - 1));  // only the last nuclide can be stable
      const double lambda = (meanlife > 0.) ? 1. / meanlife : 0.;

      decaypaths[decaypathindex].lambdas.push_back(lambda);
    }

    // the nuclide past the end of the path is a used as a sink, so treat it as stable (even if it's not)
    decaypaths[decaypathindex].lambdas.push_back(0.);

    decaypaths[decaypathindex].branchproduct = calculate_decaypath_branchproduct(decaypathindex);
  }
  decaypaths.shrink_to_fit();
}

static void filter_unused_nuclides(const std::vector<int> &custom_zlist, const std::vector<int> &custom_alist,
                                   std::vector<struct nuclide> &standard_nuclides) {
  // remove nuclides that are not a standard or custom input-specified nuclide, or connected to these by decays
  nuclides.erase(
      std::remove_if(nuclides.begin(), nuclides.end(),
                     [&](const auto &nuc) {
                       // keep nucleus if it is in the standard list
                       if (std::any_of(standard_nuclides.begin(), standard_nuclides.end(), [&](const auto &stdnuc) {
                             return (stdnuc.z == nuc.z) && (stdnuc.a == nuc.a);
                           })) {
                         return false;
                       }
                       // keep nucleus if it is in the custom list
                       for (size_t i = 0; i < custom_zlist.size(); i++) {
                         if ((nuc.z == custom_zlist[i]) && (nuc.a == custom_alist[i])) {
                           return false;
                         }
                       }

                       // keep if it is connected by decays to one of the standard or custom input-specified nuclides
                       for (const auto &decaypath : decaypaths) {
                         if (decaypath.final_daughter_z() == nuc.z && decaypath.final_daughter_a() == nuc.a) {
                           return false;
                         }

                         for (size_t i = 0; i < decaypath.z.size(); i++) {
                           if (decaypath.z[i] == nuc.z && decaypath.a[i] == nuc.a) {
                             // decay path starts with input nuc and nuc is in the decay path, so keep it
                             return false;
                           };
                         }
                       }
                       printout("removing unused nuclide (Z=%d)%s-%d\n", nuc.z, get_elname(nuc.z), nuc.a);
                       return true;
                     }),
      nuclides.end());
  nuclides.shrink_to_fit();

  // update the nuclide indicies in the decay paths after we possibly removed some nuclides
  for (auto &decaypath : decaypaths) {
    for (int i = 0; i < get_decaypathlength(decaypath); i++) {
      decaypath.nucindex[i] = get_nucindex(decaypath.z[i], decaypath.a[i]);
    }
  }
}

auto get_nucstring_z(const std::string &strnuc) -> int
// convert something like Ni56 to integer 28
{
  std::string elcode = strnuc;
  elcode.erase(std::remove_if(elcode.begin(), elcode.end(), &isdigit), elcode.end());

  for (int z = 0; z <= Z_MAX; z++) {
    if (strcmp(elcode.c_str(), get_elname(z)) == 0)  // first to letters match el symbol
    {
      return z;
    }
  }
  printout("Could not get atomic number of '%s' '%s'\n", strnuc.c_str(), elcode.c_str());
  assert_always(false);  // could not match to an element
  return -1;
}

auto get_nucstring_a(const std::string &strnuc) -> int
// convert something like Ni56 to integer 56
{
  // find first digit character
  size_t i = 0;
  for (; i < strnuc.length(); i++) {
    if (isdigit(strnuc[i]) != 0) {
      break;
    }
  }

  // remove the non-digit charts
  const std::string strmassnum = strnuc.substr(i);

  const int a = std::atoi(strmassnum.c_str());
  assert_always(a > 0);
  return a;
}

void init_nuclides(const std::vector<int> &custom_zlist, const std::vector<int> &custom_alist) {
  // add all nuclides and decays, and later trim any irrelevant ones (not connected to input-specified nuclei) later

  assert_always(custom_zlist.size() == custom_alist.size());

  // Ni57
  nuclides.push_back({.z = 28, .a = 57, .meanlife = 51.36 * 60});
  nuclides.back().endecay_positron = 0.354 * MEV;
  nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 0.436;
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1. - 0.436;

  // Ni56
  nuclides.push_back({.z = 28, .a = 56, .meanlife = 8.80 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;

  // Co56
  nuclides.push_back({.z = 27, .a = 56, .meanlife = 113.7 * DAY});
  nuclides.back().endecay_positron = 0.63 * MEV;
  nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 0.19;
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1 - 0.19;

  // Cr48
  nuclides.push_back({.z = 24, .a = 48, .meanlife = 1.29602 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;

  // V48
  nuclides.push_back({.z = 23, .a = 48, .meanlife = 23.0442 * DAY});
  nuclides.back().endecay_positron = 0.290 * MEV * 0.499;
  nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 1.;

  // Co57
  nuclides.push_back({.z = 27, .a = 57, .meanlife = 392.03 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;

  // Fe52
  nuclides.push_back({.z = 26, .a = 52, .meanlife = 0.497429 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;

  // Mn52
  nuclides.push_back({.z = 25, .a = 52, .meanlife = 0.0211395 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;

  auto standard_nuclides = nuclides;

  if (!custom_alist.empty()) {
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

      assert_always(!nuc_exists(z, a));
      nuclides.push_back({.z = z, .a = a, .meanlife = tau_sec});
      nuclides.back().branchprobs[DECAYTYPE_BETAMINUS] = 1.;
      nuclides.back().endecay_q[DECAYTYPE_BETAMINUS] = q_mev * MEV;
      nuclides.back().endecay_electron = e_elec_mev * MEV;
      nuclides.back().endecay_gamma = e_gamma_mev * MEV;
      // printout("betaminus file: Adding (Z=%d)%s-%d endecay_electron %g endecay_gamma %g tau_s %g\n",
      //          z, get_elname(z), a, e_elec_mev, e_gamma_mev, tau_sec);
      assert_always(e_elec_mev >= 0.);
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

      bool const keeprow = ((branch_alpha > 0. || branch_beta > 0.) && halflife > 0.);
      if (keeprow) {
        const double tau_sec = halflife / log(2);
        int alphanucindex = -1;
        if (nuc_exists(z, a)) {
          alphanucindex = get_nucindex(z, a);
          // printout("compare z %d a %d e_gamma_mev1 %g e_gamma_mev2 %g\n", z, a, nucdecayenergygamma(z, a) / MEV,
          // e_gamma_mev); printout("compare z %d a %d tau1 %g tau2 %g\n", z, a, get_meanlife(z, a), tau_sec);
          // printout("compare z %d a %d e_beta_mev1 %g e_beta_mev2 %g\n", z, a, nuclides[get_nucindex(z,
          // a)].endecay_positron / MEV, e_beta_mev);
        } else {
          nuclides.push_back({.z = z, .a = a, .meanlife = tau_sec, .endecay_gamma = e_gamma_mev * MEV});
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

  // add any extra nuclides that were specified but not in the decay data files
  for (int i = 0; i < static_cast<int>(custom_alist.size()); i++) {
    const int z = custom_zlist[i];
    const int a = custom_alist[i];
    if (!nuc_exists(z, a)) {
      // printout("Adding Z %d A %d with no decay data (assuming stable)\n", z, a);
      nuclides.push_back({.z = z, .a = a, .meanlife = -1});
    }
  }

  printout("Number of nuclides before filtering: %d\n", get_num_nuclides());
  find_decaypaths(custom_zlist, custom_alist, standard_nuclides);
  filter_unused_nuclides(custom_zlist, custom_alist, standard_nuclides);

  printout("Number of nuclides:  %d\n", get_num_nuclides());

  int maxdecaypathlength = 0;
  for (auto &decaypath : decaypaths) {
    maxdecaypathlength = std::max(maxdecaypathlength, get_decaypathlength(decaypath));
  }
  printout("Number of decay paths: %d (max length %d)\n", get_num_decaypaths(), maxdecaypathlength);

  /// Read in data for gamma ray lines and make a list of them in energy order.
  gammapkt::init_gamma_linelist();

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

static auto sample_decaytime(const int decaypathindex, const double tdecaymin, const double tdecaymax) -> double {
  double tdecay = -1;
  const double t_model = grid::get_t_model();
  // rejection method. draw random times until they are within the correct range.
  while (tdecay <= tdecaymin || tdecay >= tdecaymax) {
    tdecay = t_model;  // can't decay before initial model snapshot time

    for (int i = 0; i < get_decaypathlength(decaypathindex); i++) {
      const int nucindex = decaypaths[decaypathindex].nucindex[i];
      const double zrand = rng_uniform_pos();
      tdecay += -get_meanlife(nucindex) * log(zrand);
    }
  }
  return tdecay;
}

static auto calculate_decaychain(const double firstinitabund, const std::vector<double> &lambdas,
                                 const int num_nuclides, const double timediff, const bool useexpansionfactor)
    -> double {
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
            (1 + 1 / lambdas[j] / timediff) * exp(-timediff * lambdas[j]) - 1. / lambdas[j] / timediff;
        sum += sumtermtop / denominator;
      }
    }
  }

  const double lastabund = firstinitabund * lambdaproduct * sum;

  return lastabund;
}

static auto get_nuc_massfrac(const int modelgridindex, const int z, const int a, const double time) -> double
// Get the mass fraction of a nuclide accounting for all decays including those of its parent and grandparent.
// e.g., Co56 abundance may first increase with time due to Ni56 decays, then decease due to Co56 decay
// Can be called for stable nuclides that are one daughters of the radioactive nuclide list e.g., Fe56
// For stable nuclides, abundance returned only comes from other decays (some could be included in init model elem frac)
{
  assert_always(time >= 0.);

  const double t_afterinit = time - grid::get_t_model();
  const int nucindex = get_nucindex_or_neg_one(z, a);
  const bool nuc_exists_z_a = (nucindex >= 0);

  // decay chains include all paths from radionuclides to other radionuclides (including trivial size-one chains)

  double nuctotal = 0.;  // abundance or decay rate, depending on mode parameter
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
    const int z_end = decaypaths[decaypathindex].z[get_decaypathlength(decaypathindex) - 1];
    const int a_end = decaypaths[decaypathindex].a[get_decaypathlength(decaypathindex) - 1];

    // match 4He abundance to alpha decay of any nucleus (no continue), otherwise check daughter nuclide matches
    if (z != 2 || a != 4 ||
        decaypaths[decaypathindex].decaytypes[get_decaypathlength(decaypathindex) - 1] != decaytypes::DECAYTYPE_ALPHA) {
      if (nuc_exists_z_a && (z_end != z || a_end != a))  // requested nuclide is in network, so match last nuc in chain
      {
        continue;
      }

      // if requested nuclide is not in network then match daughter of last nucleus in chain
      if (!nuc_exists_z_a && !nuc_is_parent(z_end, a_end, z, a)) {
        continue;
      }
    }

    const int z_top = decaypaths[decaypathindex].z[0];
    const int a_top = decaypaths[decaypathindex].a[0];
    const int nucindex_top = decaypaths[decaypathindex].nucindex[0];

    const double top_initabund = grid::get_modelinitradioabund(modelgridindex, nucindex_top) / nucmass(z_top, a_top);
    assert_always(top_initabund >= 0.);
    if (top_initabund <= 0.) {
      continue;
    }

    int const decaypathlength = get_decaypathlength(decaypathindex);

    int fulldecaypathlength = decaypathlength;
    // if the nuclide is out of network, it's one past the end of the chain
    if (!nuc_exists_z_a || (z == 2 && a == 4 &&
                            decaypaths[decaypathindex].decaytypes[get_decaypathlength(decaypathindex) - 1] ==
                                decaytypes::DECAYTYPE_ALPHA)) {
      fulldecaypathlength = decaypathlength + 1;
    }

    const double massfraccontrib = (decaypaths[decaypathindex].branchproduct *
                                    calculate_decaychain(top_initabund, decaypaths[decaypathindex].lambdas,
                                                         fulldecaypathlength, t_afterinit, false) *
                                    nucmass(z, a));
    // assert_always(massfraccontrib >= 0.);
    nuctotal += massfraccontrib;
  }

  // for stable nuclei in the network, we need to contribute the initial abundance
  if (nuc_exists_z_a && get_meanlife(nucindex) <= 0.) {
    nuctotal += grid::get_modelinitradioabund(modelgridindex, nucindex);
  }

  return nuctotal;
}

static auto get_endecay_to_tinf_per_ejectamass_at_time(const int modelgridindex, const int decaypathindex,
                                                       const double time) -> double
// returns decay energy [erg/g] that would be released from time tstart [s] to time infinity by a given decaypath
{
  // e.g. Ni56 -> Co56, represents the decay of Co56 nuclei
  // that were produced from Ni56 in the initial abundance.
  // Decays from Co56 due to the initial abundance of Co56 are not counted here,
  // nor is the energy from Ni56 decays

  assert_testmodeonly(decaypathindex >= 0);
  assert_testmodeonly(decaypathindex < get_num_decaypaths());

  const int z_top = decaypaths[decaypathindex].z[0];
  const int a_top = decaypaths[decaypathindex].a[0];
  const int nucindex_top = decaypaths[decaypathindex].nucindex[0];
  // if it's a single-nuclide decay chain, then contribute the initial abundance, otherwise contribute
  // all ancestors

  const double top_initabund = grid::get_modelinitradioabund(modelgridindex, nucindex_top) / nucmass(z_top, a_top);
  if (top_initabund <= 0.) {
    return 0.;
  }
  assert_testmodeonly(top_initabund >= 0.);

  const int decaypathlength = get_decaypathlength(decaypathindex);

  const double t_afterinit = time - grid::get_t_model();

  // count the number of chain-top nuclei that haven't decayed past the end of the chain

  const double abund_endplusone =
      calculate_decaychain(top_initabund, decaypaths[decaypathindex].lambdas, decaypathlength + 1, t_afterinit, false);
  const double ndecays_remaining = decaypaths[decaypathindex].branchproduct * (top_initabund - abund_endplusone);

  // // alternative: add up the ancestor abundances that will eventually cause decays at the end of chain
  // double ndecays_remaining = 0.;
  // for (int c = 1; c <= decaypathlength; c++)
  // {
  //   ndecays_remaining += calculate_decaychain(top_initabund, decaypaths[decaypathindex].lambdas, c, t_afterinit);
  // }

  const double endecay = ndecays_remaining * get_decaypath_lastnucdecayenergy(decaypathindex);

  return endecay;
}

auto get_endecay_per_ejectamass_t0_to_time_withexpansion_chain_numerical(const int modelgridindex,
                                                                         const int decaypathindex, const double tstart)
    -> double
// just here as as check on the analytic result from get_endecay_per_ejectamass_t0_to_time_withexpansion()
// this version does an Euler integration
{
  double min_meanlife = -1;
  for (int i = 0; i < get_decaypathlength(decaypathindex); i++) {
    const double meanlife = get_meanlife(decaypaths[decaypathindex].nucindex[i]);
    if (min_meanlife < 0. or meanlife < min_meanlife) {
      min_meanlife = meanlife;
    }
  }

  const int nsteps = static_cast<int>(ceil((tstart - grid::get_t_model()) / min_meanlife) *
                                      100000);  // min steps across the meanlifetime
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

auto get_endecay_per_ejectamass_t0_to_time_withexpansion(const int modelgridindex, const double tstart) -> double
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

    // const double numerator = calculate_decaychain(1., meanlifetimes, decaypathlength + 1, tdiff, true);
    // const double factor = numerator / calculate_decaychain(1., meanlifetimes, decaypathlength + 1, tdiff, false);
    // printout("  Analytical expansion factor: %g\n", factor);

    const int z_top = decaypaths[decaypathindex].z[0];
    const int a_top = decaypaths[decaypathindex].a[0];
    const int nucindex_top = decaypaths[decaypathindex].nucindex[0];

    const double top_initabund = grid::get_modelinitradioabund(modelgridindex, nucindex_top) / nucmass(z_top, a_top);

    const double chain_endecay = (decaypaths[decaypathindex].branchproduct *
                                  calculate_decaychain(top_initabund, decaypaths[decaypathindex].lambdas,
                                                       decaypathlength + 1, tstart - grid::get_t_model(), true) *
                                  get_decaypath_lastnucdecayenergy(decaypathindex));

    // printout("  Analytical chain_endecay: %g\n", chain_endecay);
    tot_endecay += chain_endecay;
  }

  return tot_endecay;
}

static auto get_endecay_per_ejectamass_between_times(const int mgi, const int decaypathindex, double tlow, double thigh)
    -> double
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

static auto calculate_simtime_endecay_per_ejectamass(const int mgi, const int decaypathindex) -> double
// calculate the decay energy released during the simulation time per unit mass [erg/g]
{
  assert_testmodeonly(mgi < grid::get_npts_model());

  if constexpr (!INITIAL_PACKETS_ON) {
    // get decay energy released from t=tmin to tmax
    return get_endecay_per_ejectamass_between_times(mgi, decaypathindex, globals::tmin, globals::tmax);
  } else {
    // get decay energy released from t=0 to tmax
    return get_endecay_per_ejectamass_between_times(mgi, decaypathindex, grid::get_t_model(), globals::tmax);
  }
}

static auto get_simtime_endecay_per_ejectamass(const int mgi, const int decaypathindex) -> double
// get the decay energy released during the simulation time per unit mass [erg/g]
{
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(mgi);
  const double chainendecay = decaypath_energy_per_mass[nonemptymgi * get_num_decaypaths() + decaypathindex];
  assert_testmodeonly(chainendecay >= 0.);
  assert_testmodeonly(std::isfinite(chainendecay));
  return chainendecay;
}

static auto get_decaypath_power_per_ejectamass(const int decaypathindex, const int modelgridindex, const double time)
    -> double
// total decay power per mass [erg/s/g] for a given decaypath
{
  // only decays at the end of the chain contributed from the initial abundance of the top of the chain are counted
  // (these can be can be same for a chain of length one)

  const int z_top = decaypaths[decaypathindex].z[0];
  const int a_top = decaypaths[decaypathindex].a[0];
  const int nucindex_top = decaypaths[decaypathindex].nucindex[0];

  const double top_initabund = grid::get_modelinitradioabund(modelgridindex, nucindex_top);
  assert_always(top_initabund >= 0.) if (top_initabund <= 0.) { return 0.; }

  const int nucindex_end = decaypaths[decaypathindex].nucindex[get_decaypathlength(decaypathindex) - 1];

  const double t_afterinit = time - grid::get_t_model();

  int const decaypathlength = get_decaypathlength(decaypathindex);

  // contribution to the end nuclide abundance from the top of chain (could be a length-one chain Z,A_top = Z,A_end
  // so contribution would be from init abundance only)
  const double endnucabund =
      decaypaths[decaypathindex].branchproduct *
      calculate_decaychain(top_initabund, decaypaths[decaypathindex].lambdas, decaypathlength, t_afterinit, false);

  const double endecay = get_decaypath_lastnucdecayenergy(decaypathindex);

  const double decaypower = endecay * endnucabund / get_meanlife(nucindex_end) / nucmass(z_top, a_top);

  assert_always(decaypower >= 0.);
  assert_always(std::isfinite(decaypower));

  return decaypower;
}

auto get_modelcell_simtime_endecay_per_mass(const int mgi) -> double
// get the density at time tmin of decay energy that will
// be released during the simulation time range [erg/cm3]
{
  double endecay_per_mass = 0.;
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
    endecay_per_mass += get_simtime_endecay_per_ejectamass(mgi, decaypathindex);
  }
  return endecay_per_mass;
}

void setup_decaypath_energy_per_mass() {
  const int nonempty_npts_model = grid::get_nonempty_npts_model();
  printout(
      "[info] mem_usage: decaypath_energy_per_mass[nonempty_npts_model*num_decaypaths] occupies %.1f MB (node "
      "shared)...",
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

  printout("Calculating decaypath_energy_per_mass for all cells...");
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

void free_decaypath_energy_per_mass() {
#ifdef MPI_ON
  if (win_decaypath_energy_per_mass != MPI_WIN_NULL) {
    printout("[info] mem_usage: decaypath_energy_per_mass was freed\n");
    MPI_Win_free(&win_decaypath_energy_per_mass);
    win_decaypath_energy_per_mass = MPI_WIN_NULL;
  }
#else
  if (decaypath_energy_per_mass != nullptr) {
    printout("[info] mem_usage: decaypath_energy_per_mass was freed\n");
    free(decaypath_energy_per_mass);
    decaypath_energy_per_mass = nullptr;
  }
#endif
}

auto get_particle_injection_rate(const int modelgridindex, const double t, const int decaytype) -> double
// energy release rate in form of kinetic energy of positrons, electrons, and alpha particles in [erg/s/g]
{
  double dep_sum = 0.;
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
    const int z = get_nuc_z(nucindex);
    if (z < 1) {
      continue;
    }
    const int a = get_nuc_a(nucindex);
    const double meanlife = get_meanlife(nucindex);
    if (meanlife < 0.) {
      continue;
    }
    const double en_particles = nucdecayenergyparticle(nucindex, decaytype);
    if (en_particles > 0.) {
      const double nucdecayrate =
          get_nuc_massfrac(modelgridindex, z, a, t) / meanlife * get_nuc_decaybranchprob(nucindex, decaytype);
      assert_always(nucdecayrate >= 0);
      dep_sum += nucdecayrate * en_particles / nucmass(z, a);
    }
  }

  assert_always(std::isfinite(dep_sum));

  return dep_sum;
}

auto get_qdot_modelcell(const int modelgridindex, const double t, const int decaytype) -> double
// energy release rate [erg/s/g] including everything (even neutrinos that are ignored elsewhere)
{
  double qdot = 0.;
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
    const int z = get_nuc_z(nucindex);
    if (z < 1) {
      continue;
    }
    const int a = get_nuc_a(nucindex);
    const double meanlife = get_meanlife(nucindex);
    if (meanlife <= 0) {
      continue;
    }
    const double q_decay = nucdecayenergyqval(nucindex, decaytype) * get_nuc_decaybranchprob(nucindex, decaytype);
    if (q_decay <= 0.) {
      continue;
    }
    const double nucdecayrate = get_nuc_massfrac(modelgridindex, z, a, t) / meanlife;
    assert_always(nucdecayrate >= 0);
    qdot += nucdecayrate * q_decay / nucmass(z, a);
  }

  return qdot;
}

auto get_global_etot_t0_tinf() -> double {
  double etot_tinf = 0.;
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
    const int z_top = decaypaths[decaypathindex].z[0];
    const int a_top = decaypaths[decaypathindex].a[0];

    etot_tinf += (decaypaths[decaypathindex].branchproduct * grid::get_totmassradionuclide(z_top, a_top) /
                  nucmass(z_top, a_top) * get_decaypath_lastnucdecayenergy(decaypathindex));
  }
  return etot_tinf;
}

void update_abundances(const int modelgridindex, const int timestep, const double t_current)
/// Updates the mass fractions of elements associated with the decay sequence
/// Parameters: - modelgridindex: the grid cell for which to update the abundances
///             - t_current: current time (here mid of current timestep)
{
  assert_always(!globals::homogeneous_abundances);  // no longer supported

  printout("update_abundances for cell %d timestep %d\n", modelgridindex, timestep);

  for (int element = get_nelements() - 1; element >= 0; element--) {
    const int atomic_number = get_atomicnumber(element);
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
        if (!a_isotopes.contains(a)) {
          a_isotopes.insert(a);
          const double nuc_massfrac = get_nuc_massfrac(modelgridindex, atomic_number, a, t_current);
          isomassfracsum += nuc_massfrac;
          isomassfrac_on_nucmass_sum += nuc_massfrac / nucmass(atomic_number, a);
        }
      } else {
        // check if the nucleus decays off the network but into the selected element
        for (const auto decaytype : all_decaytypes) {
          const int daughter_z = decay_daughter_z(nuc_z, a, decaytype);
          const int daughter_a = decay_daughter_a(nuc_z, a, decaytype);
          if (daughter_z == atomic_number && !nuc_exists(daughter_z, daughter_a) &&
              get_nuc_decaybranchprob(nuc_z, a, decaytype) > 0.) {
            if (!a_isotopes.contains(daughter_a)) {
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

    if (atomic_number == 2 && !nuc_exists(2, 4) && (!a_isotopes.contains(4))) {
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
  //   initnucfracsum += grid::get_modelinitradioabund(modelgridindex, z, a);
  //   nucfracsum += get_nuc_massfrac(modelgridindex, z, a, t_current);
  //
  //   // printout_nuclidename(z, a);
  //   // printout(" init: %g now: %g\n", grid::get_modelinitradioabund(modelgridindex, z, a),
  //   get_nuc_massfrac(modelgridindex, z, a, t_current));
  //
  //   for (int dectypeindex = 0; dectypeindex < decaytypes::DECAYTYPE_COUNT; dectypeindex++)
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
}

void fprint_nuc_abundances(FILE *estimators_file, const int modelgridindex, const double t_current, const int element) {
  const double rho = grid::get_rho(modelgridindex);

  const int atomic_number = get_atomicnumber(element);
  std::set<int> a_isotopes;  // ensure we don't repeat isotopes
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
    const int nuc_z = get_nuc_z(nucindex);
    const int nuc_a = get_nuc_a(nucindex);
    if (nuc_z == atomic_number) {  // isotope of this element is on the network
      if (!a_isotopes.contains(nuc_a)) {
        a_isotopes.insert(nuc_a);
        // radioactive isotope of the element
        const double massfrac = get_nuc_massfrac(modelgridindex, atomic_number, nuc_a, t_current);
        if (massfrac > 0) {
          const double numberdens = massfrac / nucmass(atomic_number, nuc_a) * rho;

          fprintf(estimators_file, "  %s%d: %9.3e", get_elname(atomic_number), nuc_a, numberdens);
        }
      }
    } else {  // not the element that we want, but check if a decay produces it
      for (const auto decaytype : all_decaytypes) {
        const int daughter_z = decay_daughter_z(nuc_z, nuc_a, decaytype);
        const int daughter_a = decay_daughter_a(nuc_z, nuc_a, decaytype);
        // if the nucleus exists, it will be picked up by the upper condition
        if (daughter_z == atomic_number && !nuc_exists(daughter_z, daughter_a) &&
            get_nuc_decaybranchprob(nucindex, decaytype) > 0.) {
          if (!a_isotopes.contains(nuc_a)) {
            a_isotopes.insert(nuc_a);
            // nuclide decays into correct atomic number but outside of the radionuclide list. Daughter is assumed
            // stable
            const double massfrac = get_nuc_massfrac(modelgridindex, atomic_number, nuc_a, t_current);
            const double numberdens = massfrac / nucmass(nuc_z, nuc_a) * rho;
            fprintf(estimators_file, "  %s%d: %9.3e", get_elname(atomic_number), nuc_a, numberdens);
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

  // decay channels include all radioactive decay paths, and possibly also an initial cell energy channel
  const int num_decaychannels = num_decaypaths + ((INITIAL_PACKETS_ON && USE_MODEL_INITIAL_ENERGY) ? 1 : 0);

  auto cumulative_en_sum = std::make_unique<double[]>(num_decaychannels);
  double energysum = 0.;

  // add the radioactive decay paths
  for (int decaypathindex = 0; decaypathindex < num_decaypaths; decaypathindex++) {
    energysum += get_simtime_endecay_per_ejectamass(mgi, decaypathindex);
    cumulative_en_sum[decaypathindex] = energysum;
  }

  if (num_decaychannels > num_decaypaths) {
    // the t_model / tmin expansion factor was already applied at model read in
    // so "init" here means at tmin
    energysum += grid::get_initenergyq(mgi);
    cumulative_en_sum[num_decaychannels - 1] = energysum;
  }

  assert_testmodeonly(cumulative_en_sum[num_decaychannels - 1] > 0.);

  const double zrand_en = rng_uniform() * cumulative_en_sum[num_decaychannels - 1];

  // first cumulative_en_sum[i] such that cumulative_en_sum[i] > zrand_en
  const double *const upperval =
      std::upper_bound(&cumulative_en_sum[0], &cumulative_en_sum[num_decaychannels], zrand_en);
  const ptrdiff_t decaychannelindex = upperval - &cumulative_en_sum[0];

  assert_always(decaychannelindex >= 0);
  assert_always(decaychannelindex < num_decaychannels);

  // initial cell energy selected
  if (decaychannelindex >= num_decaypaths) {
    assert_always(decaychannelindex == num_decaypaths);  // only one non-radioactive channel for now
    assert_always(USE_MODEL_INITIAL_ENERGY);
    assert_always(INITIAL_PACKETS_ON);

    pkt_ptr->prop_time = globals::tmin;
    pkt_ptr->tdecay = globals::tmin;
    pkt_ptr->type = TYPE_RADIOACTIVE_PELLET;
    pkt_ptr->e_cmf = e0;
    pkt_ptr->nu_cmf = e0 / H;
    pkt_ptr->pellet_nucindex = -1;
    pkt_ptr->pellet_decaytype = -1;
    return;
  }

  const int decaypathindex = decaychannelindex;

  // possibly allow decays before the first timestep
  const double tdecaymin = !INITIAL_PACKETS_ON ? globals::tmin : grid::get_t_model();

  if constexpr (UNIFORM_PELLET_ENERGIES) {
    pkt_ptr->tdecay = sample_decaytime(decaypathindex, tdecaymin, globals::tmax);
    pkt_ptr->e_cmf = e0;
  } else {
    // use uniform decay time distribution and scale the packet energies instead.
    // keeping the pellet decay rate constant will give better statistics at late times
    // when very little energy and few packets are released
    const double zrand = rng_uniform();
    pkt_ptr->tdecay = zrand * tdecaymin + (1. - zrand) * globals::tmax;

    // we need to scale the packet energy up or down according to decay rate at the randomly selected time.
    // e0 is the average energy per packet for this cell and decaypath, so we scale this up or down
    // according to: decay power at this time relative to the average decay power
    const double avgpower = get_simtime_endecay_per_ejectamass(mgi, decaypathindex) / (globals::tmax - tdecaymin);
    assert_always(avgpower > 0.);
    assert_always(std::isfinite(avgpower));
    pkt_ptr->e_cmf = e0 * get_decaypath_power_per_ejectamass(decaypathindex, mgi, pkt_ptr->tdecay) / avgpower;
    assert_always(pkt_ptr->e_cmf >= 0);
    assert_always(std::isfinite(pkt_ptr->e_cmf));
  }

  // final decaying nuclide at the end of the chain
  const int pathlength = get_decaypathlength(decaypathindex);
  const int nucindex = decaypaths[decaypathindex].nucindex[pathlength - 1];
  const int decaytype = decaypaths[decaypathindex].decaytypes[pathlength - 1];

  pkt_ptr->type = TYPE_RADIOACTIVE_PELLET;
  pkt_ptr->pellet_nucindex = nucindex;
  pkt_ptr->pellet_decaytype = decaytype;

  const auto engamma = nucdecayenergygamma(nucindex);
  const auto enparticle = nucdecayenergyparticle(nucindex, decaytype);

  const double zrand = rng_uniform();
  pkt_ptr->originated_from_particlenotgamma = (zrand >= engamma / (engamma + enparticle));
  pkt_ptr->nu_cmf = enparticle / H;  // will be overwritten for gamma rays, but affects the thermalisation of particles
}

void cleanup() { free_decaypath_energy_per_mass(); }

}  // namespace decay