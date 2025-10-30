#include "decay.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <numeric>
#include <ranges>
#include <set>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "gammapkt.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "packet.h"
#include "random.h"
#include "sn3d.h"

namespace decay {

namespace {

constexpr auto elsymbols = std::array<const std::string, 119>{
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",  "Mg", "Al",  "Si", "P",   "S",
    "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",  "Cu", "Zn",  "Ga", "Ge",  "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",  "Pd", "Ag",  "Cd", "In",  "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu", "Gd",  "Tb", "Dy",  "Ho",
    "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",  "Hg", "Tl",  "Pb", "Bi",  "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",  "Bk", "Cf",  "Es", "Fm",  "Md",
    "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo"};
constexpr int Z_MAX = elsymbols.size() - 1;

struct Nuclide {
  int z{-1};  // atomic number
  int a{-1};  // mass number
  double meanlife{-1};  // mean lifetime before decay [s]
  double endecay_electron{0.};  // average energy per beta- decay in kinetic energy of emitted electons [erg]
  double endecay_positron{0.};  // average energy per beta+ decay in kinetic energy of emitted positrons [erg]
  double endecay_gamma{0.};  // average energy per decay in gamma rays [erg]
  double endecay_alpha{0.};  // average energy per alpha decay in kinetic energy of alpha particles [erg]
  std::array<double, decaytypes::DECAYTYPE_COUNT> endecay_q = {
      0.};  // Q-value (reactant minus product energy) for each decay type
  std::array<double, decaytypes::DECAYTYPE_COUNT> branchprobs = {0.};  // branch probability of each decay type
};

[[nodiscard]] constexpr auto decay_daughter_z_a(const int z_parent, const int a_parent, const int decaytype)
    -> std::tuple<int, int> {
  assert_always(decaytype >= 0);
  assert_always(decaytype < decaytypes::DECAYTYPE_COUNT);

  switch (static_cast<enum decaytypes>(decaytype)) {
    case decaytypes::DECAYTYPE_ALPHA: {
      return {z_parent - 2, a_parent - 4};  // lose two protons and two neutrons
    }
    case decaytypes::DECAYTYPE_BETAPLUS:
    case decaytypes::DECAYTYPE_ELECTRONCAPTURE: {
      return {z_parent - 1, a_parent};  // lose a proton, gain a neutron
    }
    case decaytypes::DECAYTYPE_BETAMINUS: {
      return {z_parent + 1, a_parent};  // lose a neutron, gain a proton
    }
    case decaytypes::DECAYTYPE_NONE: {
      return {-1, -1};  // no daughter
    }
    case decaytypes::DECAYTYPE_COUNT: {
      assert_always(false);
    }
  }
  return {-1, -1};  // no daughter
}

// a decay path follows the contribution from an initial nuclear abundance
// to another (daughter of last nuclide in decaypath) via decays
// every different path within the network is considered, e.g. 56Ni -> 56Co -> 56Fe is separate to 56Ni -> 56Co
struct DecayPath {
  std::vector<int> z;  // atomic number
  std::vector<int> a;  // mass number
  std::vector<int> nucindex;  // index into nuclides list
  std::vector<int> decaytypes;
  std::vector<double> lambdas;
  double branchproduct{
      0.};  // product of all branching factors along the path set by calculate_decaypath_branchproduct()
};

std::vector<Nuclide> nuclides;
std::vector<DecayPath> decaypaths;

// decaypath_energy_per_mass points to an array of length npts_model * num_decaypaths
// the index [mgi * num_decaypaths + i] will hold the decay energy per mass [erg/g] released by chain i in cell mgi
// during the simulation time range tmin to tmax
std::span<double> decaypath_energy_per_mass{};
MPI_Win win_decaypath_energy_per_mass{MPI_WIN_NULL};

// Get the probability that a decay of decaytype occurs
[[nodiscard]] auto get_nuc_decaybranchprob(const int nucindex, const int decaytype) -> double {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < get_num_nuclides());
  assert_testmodeonly(decaytype >= 0);
  assert_testmodeonly(decaytype < decaytypes::DECAYTYPE_COUNT);
  return nuclides[nucindex].branchprobs[decaytype];
}

[[nodiscard]] auto nucmass(int nucindex) -> double { return get_nuc_a(nucindex) * MH; }

[[nodiscard]] auto get_nuc_z_a(const int nucindex) -> std::tuple<int, int> {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < get_num_nuclides());
  return {nuclides[nucindex].z, nuclides[nucindex].a};
}

// get the nuclide array index from the atomic number and mass number
[[nodiscard]] auto get_nucindex_or_neg_one(const int z, const int a) -> int {
  assert_testmodeonly(get_num_nuclides() > 0);
  const auto num_nuclides = get_num_nuclides();

  for (int nucindex = 0; nucindex < num_nuclides; nucindex++) {
    if (nuclides[nucindex].z == z && nuclides[nucindex].a == a) {
      return nucindex;
    }
  }
  return -1;  // nuclide not found
}

[[nodiscard]] auto get_meanlife(const int nucindex) -> double {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < get_num_nuclides());
  return nuclides[nucindex].meanlife;
}

void printout_nuclidename(const int z, const int a) { printlog("(Z={}){}{}", z, get_elname(z), a); }

void printout_nuclidemeanlife(const int z, const int a) {
  const int nucindex = get_nucindex(z, a);
  if (get_meanlife(nucindex) > 0.) {
    printlog("[tau {:.1e}s]", get_meanlife(nucindex));
  } else {
    printlog("[stable]");
  }
}

// decay energy in the form of kinetic energy of electrons, positrons, or alpha particles,
// depending on the relevant decay type (but not including neutrinos)
[[nodiscard]] auto nucdecayenergyparticle(const int nucindex, const int decaytype) -> double {
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

// average energy (erg) per decay in the form of gammas and particles [erg]
[[nodiscard]] auto nucdecayenergytotal(const int z, const int a) -> double {
  const int nucindex = get_nucindex(z, a);
  const auto endecay_particles = std::accumulate(
      all_decaytypes.cbegin(), all_decaytypes.cend(), 0., [nucindex](const double ensum, const auto& decaytype) {
        return ensum + (nucdecayenergyparticle(nucindex, decaytype) * get_nuc_decaybranchprob(nucindex, decaytype));
      });

  return nuclides[nucindex].endecay_gamma + endecay_particles;
}

// contributed energy release per decay [erg] for decaytype (e.g. decaytypes::DECAYTYPE_BETAPLUS) (excludes neutrinos!)
[[nodiscard]] auto nucdecayenergy(const int nucindex, const int decaytype) -> double {
  const double endecay = nuclides[nucindex].endecay_gamma + nucdecayenergyparticle(nucindex, decaytype);

  return endecay;
}

[[nodiscard]] auto nucdecayenergyqval(const int nucindex, const int decaytype) -> double {
  return nuclides[nucindex].endecay_q[decaytype];
}

[[nodiscard]] auto get_num_decaypaths() -> int { return static_cast<int>(decaypaths.size()); }

[[nodiscard]] auto get_decaypathlength(const DecayPath& dpath) -> int {
  return static_cast<int>(dpath.decaytypes.size());
}
[[nodiscard]] auto get_decaypathlength(const int decaypathindex) -> int {
  return get_decaypathlength(decaypaths[decaypathindex]);
}

// return the product of all branching factors in the decay path
[[nodiscard]] auto calculate_decaypath_branchproduct(const DecayPath& decaypath) -> double {
  double branchprod = 1.;
  const auto decaypathlength = get_decaypathlength(decaypath);
  for (int i = 0; i < decaypathlength - 1; i++) {
    branchprod = branchprod * get_nuc_decaybranchprob(decaypath.nucindex[i], decaypath.decaytypes[i]);
  }
  return branchprod;
}

// a decaypath's energy is the decay energy of the last nuclide and decaytype in the chain
[[nodiscard]] auto get_decaypath_lastnucdecayenergy(const DecayPath& dpath) -> double {
  const auto secondlastindex = dpath.nucindex.size() - 2;
  assert_testmodeonly(dpath.decaytypes[secondlastindex] != DECAYTYPE_NONE);
  assert_testmodeonly(dpath.decaytypes.back() == DECAYTYPE_NONE);
  return nucdecayenergy(dpath.nucindex[secondlastindex], dpath.decaytypes[secondlastindex]);
}

[[nodiscard]] auto get_decaypath_lastnucdecayenergy(const int decaypathindex) -> double {
  return get_decaypath_lastnucdecayenergy(decaypaths[decaypathindex]);
}

[[nodiscard]] auto get_str_decaytype(const int decaytype) -> std::string {
  switch (decaytype) {
    case decaytypes::DECAYTYPE_ALPHA: {
      return "alpha";
    }
    case decaytypes::DECAYTYPE_BETAPLUS: {
      return "beta+";
    }
    case decaytypes::DECAYTYPE_ELECTRONCAPTURE: {
      return "ec";
    }
    case decaytypes::DECAYTYPE_BETAMINUS: {
      return "beta-";
    }
    case decaytypes::DECAYTYPE_NONE: {
      return "none";
    }
    default:
      return "unknown";
  }
}

void printout_decaypath(const int decaypathindex) {
  assert_always(!decaypaths.empty());
  const auto& decaypath = decaypaths[decaypathindex];
  printlog(" decaypath {}: ", decaypathindex);

  for (const auto [decaytype, z, a] : std::views::zip(decaypath.decaytypes, decaypath.z, decaypath.a)) {
    printout_nuclidename(z, a);
    printout_nuclidemeanlife(z, a);

    if (decaytype != DECAYTYPE_NONE) {
      printlog(" -> {} -> ", get_str_decaytype(decaytype));
    }
  }

  printlnlog("");
}

// follow decays at the ends of the current list of decaypaths
// to get decaypaths from all descendants
void extend_lastdecaypath() {
  const int startdecaypathindex = static_cast<int>(decaypaths.size() - 1);

  const int end_nucindex = decaypaths[startdecaypathindex].nucindex.back();
  if ((get_meanlife(end_nucindex) <= 0.)) {
    // daughter is stable: no extension possible
    return;
  }
  const int prev_end_z = decaypaths[startdecaypathindex].z.back();
  const int prev_end_a = decaypaths[startdecaypathindex].a.back();
  for (const auto decaytypeindex : all_decaytypes) {
    if (get_nuc_decaybranchprob(end_nucindex, decaytypeindex) == 0.) {
      continue;
    }

    const auto [daughter_z, daughter_a] = decay_daughter_z_a(prev_end_z, prev_end_a, decaytypeindex);
    // check for nuclide in existing path, which would indicate a loop
    for (const auto [z, a] : std::views::zip(decaypaths[startdecaypathindex].z, decaypaths[startdecaypathindex].a)) {
      if (z == daughter_z && a == daughter_a) {
        printlnlog("\nERROR: Loop found in nuclear decay chain.");
        std::abort();
      }
    }
    const auto daughter_nucindex = get_nucindex(daughter_z, daughter_a);
    decaypaths.push_back(decaypaths[startdecaypathindex]);
    decaypaths.back().z.push_back(daughter_z);
    decaypaths.back().a.push_back(daughter_a);
    decaypaths.back().nucindex.push_back(daughter_nucindex);
    decaypaths.back().decaytypes.back() = decaytypeindex;  // replace the DECAYTYPE_NONE at end with this decay type
    decaypaths.back().decaytypes.push_back(DECAYTYPE_NONE);  // add new DECAYTYPE_NONE at end

    extend_lastdecaypath();
  }
}

void find_decaypaths(const std::vector<int>& custom_zlist, const std::vector<int>& custom_alist,
                     const std::vector<Nuclide>& standard_nuclides) {
  decaypaths.clear();
  for (int startnucindex = 0; startnucindex < get_num_nuclides(); startnucindex++) {
    if (get_meanlife(startnucindex) <= 0.) {
      continue;  // skip stable nuclides as start points
    }
    const int z = get_nuc_z(startnucindex);
    const int a = get_nuc_a(startnucindex);

    for (const auto decaytype : all_decaytypes) {
      if (get_nuc_decaybranchprob(startnucindex, decaytype) == 0.) {
        continue;
      }
      bool is_custom_nuclide = false;
      for (auto i = 0z; i < std::ssize(custom_zlist); i++) {
        if ((z == custom_zlist[i]) && (a == custom_alist[i])) {
          is_custom_nuclide = true;
          break;
        }
      }
      // skip path if it doesn't start from a nuclide in the custom or standard input lists
      // i.e. the first nuclide will have zero initial abundance anyway
      if (!is_custom_nuclide && !std::ranges::any_of(standard_nuclides, [z, a](const auto& stdnuc) {
            return (z == stdnuc.z) && (a == stdnuc.a);
          })) {
        continue;
      }

      const auto [daughter_z, daughter_a] = decay_daughter_z_a(z, a, decaytype);

      decaypaths.push_back({.z = {z, daughter_z},
                            .a = {a, daughter_a},
                            .nucindex = {startnucindex, get_nucindex(daughter_z, daughter_a)},
                            .decaytypes = {decaytype, DECAYTYPE_NONE},
                            .lambdas = {},
                            .branchproduct = 0.});

      extend_lastdecaypath();  // take this single step chain and find all descendants
    }
  }

  std::ranges::SORT_OR_STABLE_SORT(decaypaths, [](const DecayPath& d1, const DecayPath& d2) {
    // true if d1 < d2
    // chains are sorted by mass number, then atomic number, then length
    const int d1_length = get_decaypathlength(d1);
    const int d2_length = get_decaypathlength(d2);
    // -1 to ignore last item, which keeps bit-identical results as before when final daughter nuclide was not included
    // TODO: it would probably be better to sort by all items in reverse order
    const int smallestpathlength = std::min(d1_length, d2_length) - 1;
    for (int i = 0; i < smallestpathlength; i++) {
      if (d1.a[i] != d2.a[i]) {
        return d1.a[i] < d2.a[i];
      }
      if (d1.z[i] != d2.z[i]) {
        return d1.z[i] < d2.z[i];
      }
    }

    // one is an extension of the other, so place the shorter one first
    return d1_length < d2_length;
  });

  for (auto& decaypath : decaypaths) {
    // all nuclei in the path (except for the last one, which is allowed to be stable) must have a mean life >0
    assert_always(std::all_of(decaypath.nucindex.cbegin(), decaypath.nucindex.cend() - 1,
                              [](const auto nucindex) { return get_meanlife(nucindex) > 0.; }));

    assert_always(decaypath.decaytypes.back() == DECAYTYPE_NONE);

    // convert mean lifetimes to decay constants
    decaypath.lambdas.resize(decaypath.nucindex.size(), -1.);
    std::ranges::transform(decaypath.nucindex, decaypath.lambdas.begin(), [](const auto nucindex) {
      const double meanlife = get_meanlife(nucindex);
      // last nuclide might be stable (meanlife <= 0.)
      const double lambda = (meanlife > 0.) ? 1. / meanlife : 0.;
      return lambda;
    });

    decaypath.branchproduct = calculate_decaypath_branchproduct(decaypath);
  }
  decaypaths.shrink_to_fit();
}

// remove nuclides that are not a standard or custom input-specified nuclide, or connected to these by decays
void filter_unused_nuclides(const std::vector<int>& custom_zlist, const std::vector<int>& custom_alist,
                            const std::vector<Nuclide>& standard_nuclides) {
  std::erase_if(nuclides, [&](const auto& nuc) {
    // keep nucleus if it is in the standard list
    if (std::ranges::any_of(standard_nuclides,
                            [&](const auto& stdnuc) { return (stdnuc.z == nuc.z) && (stdnuc.a == nuc.a); })) {
      return false;
    }
    // keep nucleus if it is in the custom list
    for (const auto [z, a] : std::views::zip(custom_zlist, custom_alist)) {
      if ((z == nuc.z) && (a == nuc.a)) {
        return false;
      }
    }

    const bool in_any_decaypath = std::ranges::any_of(decaypaths, [&nuc](const auto& decaypath) {
      for (const auto [z, a] : std::views::zip(decaypath.z, decaypath.a)) {
        if (z == nuc.z && a == nuc.a) {
          // nuc is in the decay path
          return true;
        }
      }

      // return true if nuc is the final daughter of a decay path
      return (decaypath.z.back() == nuc.z && decaypath.a.back() == nuc.a);
    });

    if (in_any_decaypath) {
      return false;
    }

    printout("removing unused nuclide (Z=%d)%s%d\n", nuc.z, get_elname(nuc.z).c_str(), nuc.a);
    return true;
  });
  nuclides.shrink_to_fit();

  // update the nuclide indices in the decay paths after we possibly removed some nuclides
  for (auto& decaypath : decaypaths) {
    std::ranges::transform(decaypath.z, decaypath.a, decaypath.nucindex.begin(),
                           [](const auto z, const auto a) { return get_nucindex(z, a); });
  }
}

auto sample_decaytime(const int decaypathindex, const double tdecaymin, const double tdecaymax) -> double {
  double tdecay = -1;
  const double t_model = grid::get_t_model();
  // rejection method. draw random times with the right distribution until they are within the correct range.
  while (tdecay <= tdecaymin || tdecay >= tdecaymax) {
    tdecay = t_model;  // can't decay before initial model snapshot time
    const int maxlength = get_decaypathlength(decaypaths[decaypathindex]) - 1;
    for (int i = 0; i < maxlength; i++) {
      tdecay +=
          -get_meanlife(decaypaths[decaypathindex].nucindex[i]) * std::log(static_cast<double>(rng_uniform_pos()));
    }
  }
  return tdecay;
}

// calculate final number abundance from multiple decays, e.g., Ni56 -> Co56 -> Fe56 (nuc[0] -> nuc[1] -> nuc[2])
// the top nuclide initial abundance is set and the chain-end abundance is returned (all intermediates nuclides
// are assumed to start with zero abundance)
// note: first and last can be nuclide can be the same if num_nuclides==1, reducing to simple decay formula
//
// timediff:           time elapsed for decays [seconds]
// numnuclides:        number of items in lambdas to use
// lambdas:            array of 1/(mean lifetime) for nuc[0]..nuc[num_nuclides-1]  [seconds^-1]
// useexpansionfactor: if true, return a modified 'abundance' at the end of the chain, with a weighting factor
//                          accounting for adiabatic loss from expansion since the decays occurred
//                          (This is needed to get the initial temperature)
constexpr auto calculate_decaychain(const double firstinitabund, const std::vector<double>& lambdas,
                                    const int num_nuclides, const double timediff, const bool useexpansionfactor)
    -> double {
  assert_testmodeonly(num_nuclides >= 1);
  assert_testmodeonly(std::ssize(lambdas) >= num_nuclides);

  double lambdaproduct = 1.;
  for (int j = 0; j < (num_nuclides - 1); j++) {
    lambdaproduct *= lambdas[j];
  }

  double sum = 0;
  for (int j = 0; j < num_nuclides; j++) {
    const auto lambda_j = lambdas[j];
    double denominator = 1.;
    for (int p = 0; p < num_nuclides; p++) {
      if (p != j) {
        denominator *= (lambdas[p] - lambda_j);
      }
    }

    if (!useexpansionfactor) {
      // get abundance output
      sum += exp(-lambda_j * timediff) / denominator;
    } else {
      if (lambda_j > 0.) {
        const double sumtermtop =
            ((1 + 1 / lambda_j / timediff) * exp(-timediff * lambda_j)) - (1. / lambda_j / timediff);
        sum += sumtermtop / denominator;
      }
    }
  }

  const double lastabund = firstinitabund * lambdaproduct * sum;

  return lastabund;
}

// Get the mass fraction of a nuclide accounting for all decays and initial abundances.
// e.g., Co56 abundance may first increase with time due to Ni56 decays, then decease due to Co56 decay
auto get_nuc_massfrac(const int nonemptymgi, const int nucindex, const double time) -> double {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  assert_always(time >= 0.);

  const double t_afterinit = time - grid::get_t_model();
  const auto [z, a] = get_nuc_z_a(nucindex);
  const bool nuc_is_stable = (get_meanlife(nucindex) <= 0.);
  const auto lambda = nuc_is_stable ? 0. : 1. / get_meanlife(nucindex);

  double nuctotal = 0.;  // abundance or decay rate, depending on mode parameter
  for (const auto& decaypath : decaypaths) {
    const auto last_decaytype = decaypath.decaytypes[decaypath.nucindex.size() - 2];
    // match 4He abundance to alpha decay of any nucleus (no continue), otherwise check daughter nuclide matches
    if (z != 2 || a != 4 || last_decaytype != decaytypes::DECAYTYPE_ALPHA) {
      if (decaypath.z.back() != z || decaypath.a.back() != a) {
        continue;
      }
    }

    const int nucindex_top = decaypath.nucindex[0];

    const double top_initabund = grid::get_modelinitnucmassfrac(modelgridindex, nucindex_top) / nucmass(nucindex_top);
    assert_always(top_initabund >= 0.);
    if (top_initabund <= 0.) {
      continue;
    }

    const int decaypathlength = get_decaypathlength(decaypath);
    auto lambdas = std::vector<double>(decaypath.lambdas);
    if (z == 2 && a == 4) {
      // treat the end nuclide as stable He4
      lambdas.back() = 0.;
    }

    const double massfraccontrib =
        (decaypath.branchproduct * calculate_decaychain(top_initabund, lambdas, decaypathlength, t_afterinit, false) *
         nucmass(nucindex));
    nuctotal += massfraccontrib;
  }

  // add the initial abundance
  nuctotal += grid::get_modelinitnucmassfrac(modelgridindex, nucindex) * exp(-t_afterinit * lambda);

  return nuctotal;
}

// Get the decay energy [erg/g] that would be released from time tstart [s] to time infinity by a given decaypath
// e.g. Ni56 -> Co56, represents the decays of Co56 nuclei that were produced from Ni56 in the initial abundance.
// Decays from Co56 due to the initial abundance of Co56 are not counted here, nor is the energy from Ni56 decays
auto get_endecay_to_tinf_per_ejectamass_at_time(const int modelgridindex, const int decaypathindex, const double time)
    -> double {
  assert_testmodeonly(decaypathindex >= 0);
  assert_testmodeonly(decaypathindex < get_num_decaypaths());
  const auto& decaypath = decaypaths[decaypathindex];

  const auto nucmass_top = nucmass(decaypath.nucindex[0]);
  const int nucindex_top = decaypath.nucindex[0];

  const double top_initabund = grid::get_modelinitnucmassfrac(modelgridindex, nucindex_top) / nucmass_top;
  if (top_initabund <= 0.) {
    return 0.;
  }
  assert_testmodeonly(top_initabund >= 0.);

  const int decaypathlength = get_decaypathlength(decaypathindex);

  const double t_afterinit = time - grid::get_t_model();

  // count the number of chain-top nuclei that haven't decayed past the end of the chain
  auto lambdas = std::vector<double>(decaypath.lambdas);
  // treat the end nuclide as stable to count how many got produced
  lambdas[decaypathlength - 1] = 0.;

  const double abund_endsink = calculate_decaychain(top_initabund, lambdas, decaypathlength, t_afterinit, false);
  const double ndecays_remaining = decaypath.branchproduct * (top_initabund - abund_endsink);
  // TODO ensure non-negative due to numerical precision?

  const double endecay = ndecays_remaining * get_decaypath_lastnucdecayenergy(decaypathindex);

  return endecay;
}

// Simple Euler integration as a check on the analytic result from
// get_endecay_per_ejectamass_tmodel_to_time_withexpansion()
auto get_endecay_per_ejectamass_tmodel_to_time_withexpansion_chain_numerical(const int nonemptymgi,
                                                                             const int decaypathindex,
                                                                             const double tstart) -> double {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double min_meanlife =
      std::ranges::min(decaypaths[decaypathindex].nucindex |
                       std::views::transform([](const auto nucindex) { return get_meanlife(nucindex); }) |
                       std::views::filter([](const auto meanlife) { return meanlife > 0.; }));
  // min steps across the meanlifetime
  const int nsteps = static_cast<int>(ceil((tstart - grid::get_t_model()) / min_meanlife) * 100000);

  double chain_endecay = 0.;
  double last_chain_endecay = -1.;
  double last_t = -1.;
  for (int i = 0; i < nsteps; i++) {
    const double t = grid::get_t_model() + ((tstart - grid::get_t_model()) * i / nsteps);
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

  printlnlog("  chain_endecay:              {:g}", chain_endecay);
  printlnlog("  chain_endecay_noexpansion:  {:g}", chain_endecay_noexpansion);
  printlnlog("  expansion energy factor:    {:g}", chain_endecay / chain_endecay_noexpansion);

  return chain_endecay;
}

// get decay energy per mass [erg/g] released by a decaypath between times tlow [s] and thigh [s]
auto get_endecay_per_ejectamass_between_times(const int mgi, const int decaypathindex, const double tlow,
                                              const double thigh) -> double {
  assert_testmodeonly(mgi < grid::get_npts_model());
  assert_testmodeonly(tlow <= thigh);
  const double energy_tlow = get_endecay_to_tinf_per_ejectamass_at_time(mgi, decaypathindex, tlow);
  const double energy_thigh = get_endecay_to_tinf_per_ejectamass_at_time(mgi, decaypathindex, thigh);
  assert_always(energy_tlow >= energy_thigh);
  const double endiff = energy_tlow - energy_thigh;
  assert_always(std::isfinite(endiff));
  return endiff;
}

// get the decay energy released during the simulation [(tmodel if initial packets else tmin) to tmax] per unit mass
// [erg/g]
auto get_simtime_endecay_per_ejectamass(const int nonemptymgi, const int decaypathindex) -> double {
  const double chainendecay = decaypath_energy_per_mass[(nonemptymgi * get_num_decaypaths()) + decaypathindex];
  assert_testmodeonly(chainendecay >= 0.);
  assert_testmodeonly(std::isfinite(chainendecay));
  return chainendecay;
}

// Get the total decay power per mass [erg/s/g] for a given decaypath
// We only count the power from the last decay in the chain to avoid double counting of decay energy (all sub paths are
// handled separately)
[[nodiscard]] auto get_decaypath_power_per_ejectamass(const int decaypathindex, const int nonemptymgi,
                                                      const double time) -> double {
  // only decays at the end of the chain contributed from the initial abundance of the top of the chain are counted
  // (these can be can be same for a chain of length one)

  const auto& decaypath = decaypaths[decaypathindex];
  const int nucindex_top = decaypath.nucindex[0];
  const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);

  const double top_initabund = grid::get_modelinitnucmassfrac(modelgridindex, nucindex_top);
  assert_always(top_initabund >= 0.);
  if (top_initabund <= 0.) {
    return 0.;
  }

  const int lastdecay_nucindex = decaypath.nucindex[get_decaypathlength(decaypathindex) - 2];

  const double t_afterinit = time - grid::get_t_model();

  const int decaypathlength = get_decaypathlength(decaypathindex);

  // contribution to the end nuclide abundance from the top of chain (could be a length-one chain Z,A_top = Z,A_end
  // so contribution would be from init abundance only)
  const double decayingnucabund =
      decaypath.branchproduct *
      calculate_decaychain(top_initabund, decaypath.lambdas, decaypathlength - 1, t_afterinit, false);

  const double endecay = get_decaypath_lastnucdecayenergy(decaypathindex);

  const double decaypower = endecay * decayingnucabund / get_meanlife(lastdecay_nucindex) / nucmass(nucindex_top);

  assert_always(decaypower >= 0.);
  assert_always(std::isfinite(decaypower));

  return decaypower;
}

auto write_nuclides_list() {
  auto nuclides_file = std::fstream("nuclides.out", std::ofstream::out | std::ofstream::trunc);
  assert_always(nuclides_file.is_open());
  nuclides_file << "#nucindex Z A\n";
  for (int nucindex = 0; nucindex < get_num_nuclides(); nucindex++) {
    nuclides_file << nucindex << ' ' << get_nuc_z(nucindex) << ' ' << get_nuc_a(nucindex) << '\n';
  }
}

}  // anonymous namespace

[[nodiscard]] auto get_decay_neutrino_frac(const int nucindex, const int decaytype) -> double {
  // subtract fraction of other decay products, nucdecayenergy() excludes neutrinos!
  const double nu_frac = 1. - (nucdecayenergy(nucindex, decaytype) / nuclides[nucindex].endecay_q[decaytype]);
  return nu_frac;
}

[[gnu::pure]] [[nodiscard]] auto get_num_nuclides() -> ptrdiff_t { return std::ssize(nuclides); }

[[nodiscard]] auto get_elname(const int z) -> std::string {
  assert_testmodeonly(z <= Z_MAX);
  return elsymbols[z];
}

[[nodiscard]] auto get_nuc_z(const int nucindex) -> int {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < get_num_nuclides());
  return nuclides[nucindex].z;
}

[[nodiscard]] auto get_nuc_a(const int nucindex) -> int {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < get_num_nuclides());
  return nuclides[nucindex].a;
}

// get the nuclide array index from the atomic number and mass number
[[nodiscard]] auto get_nucindex(const int z, const int a) -> int {
  const int nucindex = get_nucindex_or_neg_one(z, a);
  if (nucindex >= 0) {
    return nucindex;
  }
  assert_always(false);  // nuclide not found
  return -1;
}

// check if nuclide exists in the simulation
[[nodiscard]] auto nuc_exists(const int z, const int a) -> bool { return get_nucindex_or_neg_one(z, a) >= 0; }

// average energy per decay in the form of gamma rays [erg]
[[nodiscard]] __host__ __device__ auto nucdecayenergygamma(const int nucindex) -> double {
  return nuclides[nucindex].endecay_gamma;
}

[[nodiscard]] __host__ __device__ auto nucdecayenergygamma(const int z, const int a) -> double {
  return nucdecayenergygamma(get_nucindex(z, a));
}

// set average energy per decay in the form of gamma rays [erg]
void set_nucdecayenergygamma(const int nucindex, const double value) { nuclides[nucindex].endecay_gamma = value; }

// convert something like Ni56 to integer 28
auto get_nucstring_z(const std::string& strnuc) -> int {
  std::string elcode = strnuc;
  std::erase_if(elcode, &isdigit);

  for (int z = 0; z <= Z_MAX; z++) {
    if (elcode == get_elname(z)) {
      return z;
    }
  }
  assert_always(false);  // could not match to an element
  return -1;
}

// convert something like Ni56 to integer 56
auto get_nucstring_a(const std::string& strnuc) -> int {
  // find first digit character
  auto i = 0zU;
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

// add all nuclides and decays, and later trim any irrelevant ones (not connected to input-specified nuclei)
void init_nuclides(const std::vector<int>& custom_zlist, const std::vector<int>& custom_alist) {
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

  const auto standard_nuclides = nuclides;

  // any nuclides in the custom list that are not in the standard list need beta and alpha decay data

  bool use_custom_nuclides = false;
  for (auto i = 0z; i < std::ssize(custom_zlist); i++) {
    if (custom_zlist[i] < 0 || custom_alist[i] < 0) {
      continue;
    }
    const bool in_std_list = std::ranges::any_of(standard_nuclides, [=](const auto& stdnuc) {
      return (custom_zlist[i] == stdnuc.z) && (custom_alist[i] == stdnuc.a);
    });
    if (!in_std_list) {
      use_custom_nuclides = true;
      break;
    }
  }

  if (use_custom_nuclides) {
    auto fbetaminus = fstream_required("betaminusdecays.txt", std::ios::in);
    std::string line;
    while (get_noncommentline(fbetaminus, line)) {
      // energies are average per beta decay
      // columns: # A, Z, Q[MeV], E_gamma[MeV], E_elec[MeV], E_neutrino[MeV], meanlife[s]
      int a = -1;
      int z = -1;
      double q_beta_mev = 0.;
      double e_gamma_mev = 0.;
      double e_elec_mev = 0.;
      double e_neutrino = 0.;
      double tau_sec = 0.;
      std::stringstream(line) >> a >> z >> q_beta_mev >> e_gamma_mev >> e_elec_mev >> e_neutrino >> tau_sec;
      if (q_beta_mev > 0.) {
        assert_always(!nuc_exists(z, a));
        nuclides.push_back({.z = z, .a = a, .meanlife = tau_sec});
        nuclides.back().branchprobs[DECAYTYPE_BETAMINUS] = 1.;
        nuclides.back().endecay_q[DECAYTYPE_BETAMINUS] = q_beta_mev * MEV;
        nuclides.back().endecay_electron = e_elec_mev * MEV;
        nuclides.back().endecay_gamma = e_gamma_mev * MEV;
        assert_always(e_elec_mev >= 0.);
      }
    }

    auto falpha = fstream_required("alphadecays.txt", std::ios::in);
    assert_always(falpha.is_open());
    if (!nuc_exists(2, 4)) {
      nuclides.push_back({.z = 2, .a = 4, .meanlife = -1});
    }
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

      const bool keeprow = ((branch_alpha > 0. || branch_beta > 0.) && halflife > 0.);
      if (keeprow) {
        const double tau_sec = halflife / std::numbers::ln2;
        int alphanucindex = -1;
        if (nuc_exists(z, a)) {
          alphanucindex = get_nucindex(z, a);
        } else {
          nuclides.push_back({.z = z, .a = a, .meanlife = tau_sec, .endecay_gamma = e_gamma_mev * MEV});
          alphanucindex = static_cast<int>(nuclides.size() - 1);
        }
        nuclides[alphanucindex].endecay_alpha = e_alpha_mev * MEV;
        nuclides[alphanucindex].branchprobs[DECAYTYPE_BETAMINUS] = branch_beta;
        nuclides[alphanucindex].endecay_q[DECAYTYPE_BETAMINUS] = Q_total_betadec * MEV;
        nuclides[alphanucindex].branchprobs[DECAYTYPE_ALPHA] = branch_alpha;
        nuclides[alphanucindex].endecay_q[DECAYTYPE_ALPHA] = Q_total_alphadec * MEV;
      }
    }
  }

  std::set<std::tuple<int, int>> nuclides_ensure_list{};

  // add any extra nuclides that were specified but not in the decay data files
  for (const auto [z, a] : std::views::zip(custom_zlist, custom_alist)) {
    nuclides_ensure_list.insert({z, a});
  }

  // ensure any daughters nuclides are included
  for (const auto& nuc : nuclides) {
    for (const auto& decaytype : all_decaytypes) {
      if (nuc.branchprobs[decaytype] > 0.) {
        const auto [z_daughter, a_daughter] = decay_daughter_z_a(nuc.z, nuc.a, decaytype);
        nuclides_ensure_list.insert({z_daughter, a_daughter});
      }
    }
  }
  // add any required nuclides (assume they are stable)
  for (const auto& [z, a] : nuclides_ensure_list) {
    if (!nuc_exists(z, a)) {
      nuclides.push_back({.z = z, .a = a, .meanlife = -1});
    }
  }

  std::set<std::tuple<int, int>> seen_z_a{};
  for (const auto& nuc : nuclides) {
    assert_always(!seen_z_a.contains({nuc.z, nuc.a}));  // duplicate nuclide detected
    seen_z_a.insert({nuc.z, nuc.a});
  }

  printlnlog("Number of nuclides before filtering: {}", get_num_nuclides());
  find_decaypaths(custom_zlist, custom_alist, standard_nuclides);
  filter_unused_nuclides(custom_zlist, custom_alist, standard_nuclides);

  printlnlog("Number of nuclides: {}", get_num_nuclides());

  const int maxdecaypathlength = std::ranges::fold_left(decaypaths, 0, [](const int maxlen, const auto& decaypath) {
    return std::max(maxlen, get_decaypathlength(decaypath));
  });

  printlnlog("Number of decay paths: {} (max length {})", get_num_decaypaths(), maxdecaypathlength);
  constexpr bool print_decaypaths = false;
  if (print_decaypaths) {
    for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
      printout_decaypath(decaypathindex);
    }
  }

  // Read in data for gamma ray lines and make a list of them in energy order.
  gammapkt::init_gamma_data();

  // TODO: generalise this to all included nuclides
  printlnlog("decayenergy(NI56), decayenergy(CO56), decayenergy_gamma(CO56): {:g}, {:g}, {:g}",
             nucdecayenergytotal(28, 56) / MEV, nucdecayenergytotal(27, 56) / MEV, nucdecayenergygamma(27, 56) / MEV);
  printlnlog("decayenergy(NI57), decayenergy_gamma(NI57), nucdecayenergy(CO57): {:g}, {:g}, {:g}",
             nucdecayenergytotal(28, 57) / MEV, nucdecayenergygamma(28, 57) / MEV, nucdecayenergytotal(27, 57) / MEV);
  printlnlog("decayenergy(CR48), decayenergy(V48): {:g} {:g}", nucdecayenergytotal(24, 48) / MEV,
             nucdecayenergytotal(23, 48) / MEV);
  printlnlog("decayenergy(FE52), decayenergy(MN52): {:g} {:g}", nucdecayenergytotal(26, 52) / MEV,
             nucdecayenergytotal(25, 52) / MEV);

  if (globals::my_rank == 0 && !globals::simulation_continued_from_saved) {
    write_nuclides_list();
  }
}

// calculate the decay energy per unit mass [erg/g] released from time t_model (can be before tmin) to tstart,
// accounting for the photon energy loss due to expansion between time of decays and tstart (equation 18 of Lucy 2005)
auto get_endecay_per_ejectamass_tmodel_to_time_withexpansion(const int nonemptymgi, const double tstart) -> double {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  double tot_endecay = 0.;
  for (const auto& decaypath : decaypaths) {
    const int decaypathlength = get_decaypathlength(decaypath);

    const int nucindex_top = decaypath.nucindex[0];

    const double top_initabund = grid::get_modelinitnucmassfrac(modelgridindex, nucindex_top) / nucmass(nucindex_top);
    auto lambdas = std::vector<double>(decaypath.lambdas);
    // treat the end nuclide as stable to count how many got produced
    lambdas[decaypathlength - 1] = 0.;

    const double chain_endecay =
        (decaypath.branchproduct *
         calculate_decaychain(top_initabund, lambdas, decaypathlength, tstart - grid::get_t_model(), true) *
         get_decaypath_lastnucdecayenergy(decaypath));

    tot_endecay += chain_endecay;
  }

  return tot_endecay;
}

// get the decay energy that will be released during the simulation time range [erg/g]
auto get_modelcell_simtime_endecay_per_mass(const int nonemptymgi) -> double {
  double endecay_per_mass = 0.;
  for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
    endecay_per_mass += get_simtime_endecay_per_ejectamass(nonemptymgi, decaypathindex);
  }
  return endecay_per_mass;
}

void setup_decaypath_energy_per_mass() {
  const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();
  printlog(
      "[info] mem_usage: decaypath_energy_per_mass[nonempty_npts_model*num_decaypaths] occupies {:.1f} MB (node shared "
      "memory)...",
      nonempty_npts_model * get_num_decaypaths() * sizeof(double) / 1024. / 1024.);
  std::tie(decaypath_energy_per_mass, win_decaypath_energy_per_mass) =
      MPI_shared_malloc_span_keepwin<double>(nonempty_npts_model * get_num_decaypaths());
  printlnlog("done.");

  MPI_Barrier(MPI_COMM_WORLD);
  const auto time_min_decay = INITIAL_PACKETS_ON ? grid::get_t_model() : globals::tmin;
  const ptrdiff_t num_decaypaths = get_num_decaypaths();
  for (int nonemptymgi = 0; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
    if (nonemptymgi % globals::node_nprocs == globals::rank_in_node) {
      const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
      for (int decaypathindex = 0; decaypathindex < num_decaypaths; decaypathindex++) {
        decaypath_energy_per_mass[(nonemptymgi * num_decaypaths) + decaypathindex] =
            get_endecay_per_ejectamass_between_times(mgi, decaypathindex, time_min_decay, globals::tmax);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

void free_decaypath_energy_per_mass() {
  if (win_decaypath_energy_per_mass != MPI_WIN_NULL) {
    printlnlog("[info] mem_usage: decaypath_energy_per_mass was freed");
    MPI_Win_free(&win_decaypath_energy_per_mass);
    win_decaypath_energy_per_mass = MPI_WIN_NULL;
  }
  decaypath_energy_per_mass = {};
}

// energy release rate in form of kinetic energy of positrons, electrons, and alpha particles in [erg/s/g]
[[nodiscard]] auto get_particle_injection_rate(const int nonemptymgi, const double t, const int decaytype) -> double {
  double dep_sum = 0.;
  const auto num_nuclides = get_num_nuclides();
  for (int nucindex = 0; nucindex < num_nuclides; nucindex++) {
    const double meanlife = get_meanlife(nucindex);
    if (meanlife < 0.) {
      continue;
    }
    const double en_particles = nucdecayenergyparticle(nucindex, decaytype);
    if (en_particles > 0.) {
      const double nucdecayrate =
          get_nuc_massfrac(nonemptymgi, nucindex, t) / meanlife * get_nuc_decaybranchprob(nucindex, decaytype);
      assert_testmodeonly(nucdecayrate >= 0);
      dep_sum += nucdecayrate * en_particles / nucmass(nucindex);
    }
  }

  assert_always(std::isfinite(dep_sum));

  return dep_sum;
}

// energy release rate in form of gamma-rays in [erg/s/g]
[[nodiscard]] auto get_gamma_emission_rate(const int nonemptymgi, const double t) -> double {
  double eps_gamma_sum = 0.;
  const auto num_nuclides = get_num_nuclides();
  for (int nucindex = 0; nucindex < num_nuclides; nucindex++) {
    const double meanlife = get_meanlife(nucindex);
    if (meanlife < 0.) {
      continue;
    }
    const double en_gamma = nucdecayenergygamma(nucindex);
    if (en_gamma > 0.) {
      const double nucdecayrate = get_nuc_massfrac(nonemptymgi, nucindex, t) / meanlife;
      assert_testmodeonly(nucdecayrate >= 0);
      eps_gamma_sum += nucdecayrate * en_gamma / nucmass(nucindex);
    }
  }

  assert_always(std::isfinite(eps_gamma_sum));

  return eps_gamma_sum;
}

// energy release rate [erg/s/g] including everything (even neutrinos that are ignored elsewhere)
[[nodiscard]] auto get_qdot_modelcell(const int nonemptymgi, const double t, const int decaytype) -> double {
  double qdot = 0.;
  const auto num_nuclides = get_num_nuclides();
  for (int nucindex = 0; nucindex < num_nuclides; nucindex++) {
    const double meanlife = get_meanlife(nucindex);
    if (meanlife < 0.) {
      continue;
    }
    const double q_decay = nucdecayenergyqval(nucindex, decaytype) * get_nuc_decaybranchprob(nucindex, decaytype);
    if (q_decay <= 0.) {
      continue;
    }
    const double nucdecayrate = get_nuc_massfrac(nonemptymgi, nucindex, t) / meanlife;
    assert_always(nucdecayrate >= 0);
    qdot += nucdecayrate * q_decay / nucmass(nucindex);
  }

  return qdot;
}

// total decay energy [erg] that will be released from all decay paths in the model from snapshot time until time
// infinity
auto get_global_etot_tmodel_tinf() -> double {
  double etot_tinf = 0.;
  for (const auto& decaypath : decaypaths) {
    const int nucindex_top = decaypath.nucindex[0];
    const int z_top = decaypath.z[0];
    const int a_top = decaypath.a[0];
    etot_tinf += (decaypath.branchproduct * grid::get_totmassnuclide_tmodel(z_top, a_top) / nucmass(nucindex_top) *
                  get_decaypath_lastnucdecayenergy(decaypath));
  }
  assert_always(std::isfinite(etot_tinf));
  assert_always(etot_tinf > 0.);
  return etot_tinf;
}

// Update the mass fractions of elements using the current abundances of nuclides
void update_abundances(const int nonemptymgi, const double t_current) {
  for (int element = get_nelements() - 1; element >= 0; element--) {
    const int atomic_number = get_atomicnumber(element);

    // the mass fraction sum of radioactive isotopes, and stable nuclei coming from other decays for the current element
    double isomassfracsum = 0.;
    double isomassfrac_on_nucmass_sum = 0.;
    const auto num_nuclides = get_num_nuclides();
    for (int nucindex = 0; nucindex < num_nuclides; nucindex++) {
      if (get_nuc_z(nucindex) == atomic_number) {
        const double nuc_massfrac = get_nuc_massfrac(nonemptymgi, nucindex, t_current);
        isomassfracsum += nuc_massfrac;
        isomassfrac_on_nucmass_sum += nuc_massfrac / nucmass(nucindex);
      }
    }

    const double otherstablemassfrac = grid::get_otherstable_initabund(nonemptymgi, element);
    isomassfracsum += otherstablemassfrac;
    isomassfrac_on_nucmass_sum += otherstablemassfrac / globals::elements[element].initstablemeannucmass;

    grid::set_elem_abundance(nonemptymgi, element, static_cast<float>(isomassfracsum));
    if (isomassfrac_on_nucmass_sum > 0.) {
      grid::set_element_meanweight(nonemptymgi, element,
                                   static_cast<float>(isomassfracsum / isomassfrac_on_nucmass_sum));
    } else {
      // avoid a divide by zero
      grid::set_element_meanweight(nonemptymgi, element, globals::elements[element].initstablemeannucmass);
    }
  }

  grid::set_nnetot(nonemptymgi);
}

void output_nuc_abundances(std::ostream& estimators_file, const int nonemptymgi, const double t_current,
                           const int element) {
  const double rho = grid::get_rho(nonemptymgi);

  const int atomic_number = get_atomicnumber(element);
  std::set<std::tuple<int, int>> a_isotopes;  // so that we output sorted by mass number
  const auto num_nuclides = get_num_nuclides();
  for (int nucindex = 0; nucindex < num_nuclides; nucindex++) {
    const auto [nuc_z, nuc_a] = get_nuc_z_a(nucindex);
    if (nuc_z == atomic_number) {
      a_isotopes.insert({nuc_a, nucindex});
    }
  }

  for (const auto& [nuc_a, nucindex] : a_isotopes) {
    const double massfrac = get_nuc_massfrac(nonemptymgi, nucindex, t_current);
    if (massfrac > 0) {
      const double numberdens = massfrac / nucmass(nucindex) * rho;
      estimators_file << std::format("  {}{}: {:9.3e}", get_elname(atomic_number), nuc_a, numberdens);
    }
  }

  const double otherstablemassfrac = grid::get_otherstable_initabund(nonemptymgi, element);
  if (otherstablemassfrac > 0) {
    // factor to convert convert mass fraction to number density
    const double meannucmass = globals::elements[element].initstablemeannucmass;
    const double otherstable_numberdens = otherstablemassfrac / meannucmass * grid::get_rho(nonemptymgi);
    estimators_file << std::format("  {}_otherstable: {:9.3e}", get_elname(atomic_number), otherstable_numberdens);
  }
  estimators_file << '\n';
}

void setup_radioactive_pellet(const double e_cmf_per_packet, const int nonemptymgi, Packet& pkt) {
  const int num_decaypaths = get_num_decaypaths();

  // decay channels include all radioactive decay paths, and possibly also an initial cell energy channel
  const int num_decaychannels = num_decaypaths + ((INITIAL_PACKETS_ON && USE_MODEL_INITIAL_ENERGY) ? 1 : 0);

  auto cumulative_en_sum = std::vector<double>(num_decaychannels, 0.);
  double energysum = 0.;

  // add the radioactive decay paths
  for (int decaypathindex = 0; decaypathindex < num_decaypaths; decaypathindex++) {
    energysum += get_simtime_endecay_per_ejectamass(nonemptymgi, decaypathindex);
    cumulative_en_sum[decaypathindex] = energysum;
  }

  if (num_decaychannels > num_decaypaths) {
    // the t_model / tmin expansion factor was already applied at model read in
    // so "init" here means at tmin
    const auto mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    energysum += grid::get_initenergyq(mgi);
    cumulative_en_sum[num_decaychannels - 1] = energysum;
  }

  assert_testmodeonly(cumulative_en_sum[num_decaychannels - 1] > 0.);

  const double zrand_en = rng_uniform() * cumulative_en_sum[num_decaychannels - 1];

  // first decaychannelindex such that cumulative_en_sum[decaychannelindex] > zrand_en
  const int decaychannelindex =
      static_cast<int>(std::ranges::upper_bound(cumulative_en_sum, zrand_en) - cumulative_en_sum.cbegin());

  assert_always(decaychannelindex >= 0);
  assert_always(decaychannelindex < num_decaychannels);

  // initial cell energy selected
  if (decaychannelindex >= num_decaypaths) {
    assert_always(decaychannelindex == num_decaypaths);  // only one non-radioactive channel for now
    assert_always(USE_MODEL_INITIAL_ENERGY);
    assert_always(INITIAL_PACKETS_ON);

    pkt.prop_time = globals::tmin;
    pkt.tdecay = globals::tmin;
    pkt.type = TYPE_RADIOACTIVE_PELLET;
    pkt.e_cmf = e_cmf_per_packet;
    pkt.nu_cmf = e_cmf_per_packet / H;
    pkt.pellet_nucindex = -1;
    pkt.pellet_decaytype = -1;
    return;
  }

  const int decaypathindex = decaychannelindex;

  // possibly allow decays before the first timestep
  const double tdecaymin = !INITIAL_PACKETS_ON ? globals::tmin : grid::get_t_model();

  if constexpr (UNIFORM_PELLET_ENERGIES) {
    pkt.tdecay = sample_decaytime(decaypathindex, tdecaymin, globals::tmax);
    pkt.e_cmf = e_cmf_per_packet;
  } else {
    // use uniform decay time distribution and scale the packet energies instead.
    // keeping the pellet decay rate constant will give better statistics at late times
    // when very little energy and few packets are released
    const double zrand = rng_uniform();
    pkt.tdecay = (zrand * tdecaymin) + ((1. - zrand) * globals::tmax);

    // we need to scale the packet energy up or down according to decay rate at the randomly selected time.
    // e_cmf_average is the average energy per packet for this cell and decaypath, so we scale this up or down
    // according to: decay power at this time relative to the average decay power
    const double avgpower =
        get_simtime_endecay_per_ejectamass(nonemptymgi, decaypathindex) / (globals::tmax - tdecaymin);
    assert_always(avgpower > 0.);
    pkt.e_cmf =
        e_cmf_per_packet * get_decaypath_power_per_ejectamass(decaypathindex, nonemptymgi, pkt.tdecay) / avgpower;
    assert_always(pkt.e_cmf >= 0);
  }

  pkt.type = TYPE_RADIOACTIVE_PELLET;

  // final decaying nuclide at the end of the chain (one before the end, which is the daughter of the last decay)
  const int pathlength = get_decaypathlength(decaypathindex);
  pkt.pellet_nucindex = decaypaths[decaypathindex].nucindex[pathlength - 2];
  pkt.pellet_decaytype = decaypaths[decaypathindex].decaytypes[pathlength - 2];

  const auto engamma = nucdecayenergygamma(pkt.pellet_nucindex);
  const auto enparticle = nucdecayenergyparticle(pkt.pellet_nucindex, pkt.pellet_decaytype);

  pkt.originated_from_particlenotgamma = (rng_uniform() >= engamma / (engamma + enparticle));
  if (pkt.originated_from_particlenotgamma) {
    // particle (positron, electron, or alpha) emitted
    pkt.nu_cmf = enparticle / H;
  } else {
    // gamma ray emitted
    pkt.nu_cmf = gammapkt::choose_gamma_ray(pkt.pellet_nucindex);
  }
}

void cleanup() { free_decaypath_energy_per_mass(); }

}  // namespace decay
