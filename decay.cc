// Radioactive decays: nuclide properties, decay paths, time-dependent nuclide abundances
// (Bateman equation solutions), and the energy release rates and decay-time sampling used to
// set up the radioactive energy pellets.

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
#include <functional>
#include <iostream>
#include <iterator>
#include <numbers>
#include <numeric>
#include <print>
#include <ranges>
#include <set>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "gammapkt.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "mpi_logging.h"
#include "packet.h"
#include "random.h"
#include "sn3d.h"

namespace decay {

namespace {

constexpr std::array elsymbols{
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",  "Mg", "Al",  "Si", "P",   "S",
    "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",  "Cu", "Zn",  "Ga", "Ge",  "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",  "Pd", "Ag",  "Cd", "In",  "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu", "Gd",  "Tb", "Dy",  "Ho",
    "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",  "Hg", "Tl",  "Pb", "Bi",  "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",  "Bk", "Cf",  "Es", "Fm",  "Md",
    "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo",
};
constexpr int Z_MAX = elsymbols.size() - 1;

struct DecayDaughter {
  int z{-1};
  int a{-1};
  double probability{0.};
};

struct Nuclide {
  int z{-1};  // atomic number
  int a{-1};  // mass number
  double meanlife{-1};  // mean lifetime before decay [s]
  double endecay_electron{0.};  // average energy per beta- decay in kinetic energy of emitted electrons [erg]
  double endecay_positron{0.};  // average energy per beta+ decay in kinetic energy of emitted positrons [erg]
  double endecay_gamma{0.};  // average energy per decay in gamma rays [erg]
  double endecay_alpha{0.};  // average energy per alpha decay in kinetic energy of alpha particles [erg]
  double endecay_fission{0.};  // average energy per fission decay in kinetic energy of fission fragments [erg]
  std::array<double, DecayType::DECAYTYPE_COUNT> endecay_q = {
      0., 0., 0., 0., 0., 0.,
  };  // Q-value (reactant minus product energy) for each decay type
  std::array<double, DecayType::DECAYTYPE_COUNT> branchprobs = {
      0., 0., 0., 0., 0., 0.,
  };  // branch probability of each decay type

  // (Z, A, probability) of fission daughters
  std::vector<DecayDaughter> fission_daughters_z_a_prob{};  // NOLINT(readability-redundant-member-init)

  // sum of daughter probabilities for all decay types
  // default to 1.0 for the single-daughter decays and replace for fission
  double decay_daughters_probsum{1.};
};

// a decay path follows the contribution from an initial nuclear abundance
// to another (daughter of last nuclide in decaypath) via decays
// every different path within the network is considered, e.g. 56Ni -> 56Co -> 56Fe is separate to 56Ni -> 56Co
struct DecayPath {
  std::vector<int> z;  // atomic number
  std::vector<int> a;  // mass number
  std::vector<int> nucindex;  // index into nuclides list
  std::vector<DecayType> decaytypes;
  std::vector<double> lambdas;
  double branchproduct{
      0.,
  };  // product of all branching factors along the path set by calculate_decaypath_branchproduct()
};

std::vector<Nuclide> nuclides;
std::vector<DecayPath> decaypaths;
std::vector<bool> alldecaytypes_is_used;

[[nodiscard]] constexpr auto decay_daughters_z_a_prob(const int z_parent, const int a_parent, const DecayType decaytype)
    -> std::vector<DecayDaughter> {
  assert_always(decaytype >= 0);
  assert_always(decaytype < DecayType::DECAYTYPE_COUNT);

  switch (static_cast<enum DecayType>(decaytype)) {
    case DecayType::DECAYTYPE_ALPHA: {
      return {
          DecayDaughter{
              .z = z_parent - 2,
              .a = a_parent - 4,
              .probability = 1.,
          },
      };  // lose two protons and two neutrons (He4 handled separately)
    }
    case DecayType::DECAYTYPE_BETAPLUS:
    case DecayType::DECAYTYPE_ELECTRONCAPTURE: {
      return {DecayDaughter{.z = z_parent - 1, .a = a_parent, .probability = 1.}};  // lose a proton, gain a neutron
    }
    case DecayType::DECAYTYPE_BETAMINUS: {
      return {DecayDaughter{.z = z_parent + 1, .a = a_parent, .probability = 1.}};  // lose a neutron, gain a proton
    }
    case DecayType::DECAYTYPE_SPONTFISSION: {
      const auto nucindex = get_nucindex(z_parent, a_parent);
      assert_always(!nuclides[nucindex].fission_daughters_z_a_prob.empty());
      return nuclides[nucindex].fission_daughters_z_a_prob;
    }
    case DecayType::DECAYTYPE_NONE: {
      return {DecayDaughter{.z = -1, .a = -1, .probability = 0.}};  // no daughter
    }
    case DecayType::DECAYTYPE_COUNT: {
      assert_always(false);
    }
  }
  return {DecayDaughter{.z = -1, .a = -1, .probability = 0.}};  // no daughter
}

// Get the probability that a decay of decaytype occurs
[[nodiscard]] auto get_nuc_decaybranchprob(const int nucindex, const DecayType decaytype) -> double {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < std::ssize(nuclides));
  assert_testmodeonly(decaytype >= 0);
  assert_testmodeonly(decaytype < DecayType::DECAYTYPE_COUNT);
  return nuclides[nucindex].branchprobs[decaytype];
}

[[nodiscard]] auto nucmass(int nucindex) -> double { return get_nuc_a(nucindex) * MH; }

[[nodiscard]] auto get_nuc_z_a(const int nucindex) -> std::tuple<int, int> {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < std::ssize(nuclides));
  return {nuclides[nucindex].z, nuclides[nucindex].a};
}

// get the nuclide array index from the atomic number and mass number
[[nodiscard]] auto get_nucindex_or_neg_one(const int z, const int a) -> int {
  assert_testmodeonly(std::ssize(nuclides) > 0);
  const auto it = std::ranges::find_if(nuclides, [z, a](const auto& nuc) { return nuc.z == z && nuc.a == a; });
  return it != nuclides.end() ? static_cast<int>(std::ranges::distance(nuclides.begin(), it)) : -1;
}

[[nodiscard]] auto get_meanlife(const int nucindex) -> double {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < std::ssize(nuclides));
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
[[nodiscard]] auto nucdecayenergyparticle(const int nucindex, const DecayType decaytype) -> double {
  assert_testmodeonly(decaytype >= 0);
  assert_testmodeonly(decaytype < DecayType::DECAYTYPE_COUNT);

  switch (decaytype) {
    case DecayType::DECAYTYPE_ALPHA: {
      return nuclides[nucindex].endecay_alpha;
    }
    case DecayType::DECAYTYPE_BETAPLUS: {
      return nuclides[nucindex].endecay_positron;
    }
    case DecayType::DECAYTYPE_ELECTRONCAPTURE: {
      return 0.;
    }
    case DecayType::DECAYTYPE_BETAMINUS: {
      return nuclides[nucindex].endecay_electron;
    }
    case DecayType::DECAYTYPE_SPONTFISSION: {
      return nuclides[nucindex].endecay_fission;
    }
    case DecayType::DECAYTYPE_NONE: {
      return 0.;
    }
    default: {
      assert_always(false);
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
[[nodiscard]] auto nucdecayenergy(const int nucindex, const DecayType decaytype) -> double {
  const double endecay = nuclides[nucindex].endecay_gamma + nucdecayenergyparticle(nucindex, decaytype);

  return endecay;
}

[[nodiscard]] auto nucdecayenergyqval(const int nucindex, const DecayType decaytype) -> double {
  return nuclides[nucindex].endecay_q[decaytype];
}

[[nodiscard]] auto get_num_decaypaths() -> int { return static_cast<int>(decaypaths.size()); }

// a decaypath's energy is the decay energy of the last nuclide and decaytype in the chain
[[nodiscard]] auto get_decaypath_lastdecayenergy(const DecayPath& decaypath) -> double {
  const auto secondlastindex = decaypath.nucindex.size() - 2;
  assert_testmodeonly(decaypath.decaytypes[secondlastindex] != DecayType::DECAYTYPE_NONE);
  assert_testmodeonly(decaypath.lambdas[secondlastindex] > 0.);
  assert_testmodeonly(decaypath.decaytypes.back() == DecayType::DECAYTYPE_NONE);
  // Normalize decay energy by decay_daughters_probsum to properly distribute energy among probabilistic fission product
  // outcomes. This avoids multiply counting the decay energy, since each daughter product will have it's own
  // separate decaypath that differs only in the last nuclide.
  const auto nucindex = decaypath.nucindex[secondlastindex];
  return nucdecayenergy(nucindex, decaypath.decaytypes[secondlastindex]) / nuclides[nucindex].decay_daughters_probsum;
}

[[nodiscard]] auto get_str_decaytype(const DecayType decaytype) -> std::string_view {
  switch (decaytype) {
    case DecayType::DECAYTYPE_ALPHA: {
      return "alpha";
    }
    case DecayType::DECAYTYPE_BETAPLUS: {
      return "beta+";
    }
    case DecayType::DECAYTYPE_ELECTRONCAPTURE: {
      return "ec";
    }
    case DecayType::DECAYTYPE_BETAMINUS: {
      return "beta-";
    }
    case DecayType::DECAYTYPE_SPONTFISSION: {
      return "sf";
    }
    case DecayType::DECAYTYPE_NONE: {
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

    if (decaytype != DecayType::DECAYTYPE_NONE) {
      printlog(" -> {} -> ", get_str_decaytype(decaytype));
    }
  }

  printlnlog("");
}

// follow decays at the ends of the current list of decaypaths
// to get decaypaths from all descendants
void extend_lastdecaypath(std::vector<DecayPath>& localdecaypaths) {
  const auto initial_last_decaypath = localdecaypaths.back();

  const int end_nucindex = initial_last_decaypath.nucindex.back();
  if ((get_meanlife(end_nucindex) <= 0.)) {
    // daughter is stable: no extension possible
    return;
  }
  const int prev_end_z = initial_last_decaypath.z.back();
  const int prev_end_a = initial_last_decaypath.a.back();
  for (const auto decaytypeindex : all_decaytypes) {
    if (get_nuc_decaybranchprob(end_nucindex, decaytypeindex) == 0.) {
      continue;
    }

    for (const auto& daughter : decay_daughters_z_a_prob(prev_end_z, prev_end_a, decaytypeindex)) {
      // check for nuclide in existing path, which would indicate a loop
      for (const auto [z, a] : std::views::zip(initial_last_decaypath.z, initial_last_decaypath.a)) {
        if (z == daughter.z && a == daughter.a) {
          std::string chainstr;
          for (const auto [chain_z, chain_a] : std::views::zip(initial_last_decaypath.z, initial_last_decaypath.a)) {
            chainstr += std::format("(Z={},A={}) -> ", chain_z, chain_a);
          }
          printlnlog("[error] loop in nuclear decay chain: {}(Z={},A={}) already occurred. aborting", chainstr,
                     daughter.z, daughter.a);
          std::abort();
        }
      }
      const auto daughter_nucindex = get_nucindex(daughter.z, daughter.a);
      auto newdecaypath = initial_last_decaypath;
      newdecaypath.z.push_back(daughter.z);
      newdecaypath.a.push_back(daughter.a);
      newdecaypath.nucindex.push_back(daughter_nucindex);
      newdecaypath.decaytypes.back() = decaytypeindex;  // replace the DECAYTYPE_NONE at end with this decay type
      newdecaypath.decaytypes.push_back(DecayType::DECAYTYPE_NONE);  // add new DECAYTYPE_NONE at end
      newdecaypath.branchproduct *= get_nuc_decaybranchprob(end_nucindex, decaytypeindex) * daughter.probability;
      localdecaypaths.push_back(newdecaypath);

      extend_lastdecaypath(localdecaypaths);
    }
  }
}

auto find_decaypaths(const std::span<const int> custom_zlist, const std::span<const int> custom_alist,
                     const std::vector<Nuclide>& standard_nuclides) -> std::vector<DecayPath> {
  std::vector<DecayPath> localdecaypaths;
  for (int startnucindex = 0; startnucindex < std::ssize(nuclides); startnucindex++) {
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
      for (auto i = 0Z; i < std::ssize(custom_zlist); i++) {
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

      for (const auto& daughter : decay_daughters_z_a_prob(z, a, decaytype)) {
        localdecaypaths.push_back({
            .z = {z, daughter.z},
            .a = {a, daughter.a},
            .nucindex = {startnucindex, get_nucindex(daughter.z, daughter.a)},
            .decaytypes = {decaytype, DecayType::DECAYTYPE_NONE},
            .lambdas = {},
            .branchproduct = get_nuc_decaybranchprob(startnucindex, decaytype) * daughter.probability,
        });

        extend_lastdecaypath(localdecaypaths);  // take this single step chain and find all descendants
      }
    }
  }

  std::ranges::SORT_OR_STABLE_SORT(localdecaypaths, [](const DecayPath& d1, const DecayPath& d2) {
    // true if d1 < d2
    // chains are sorted by mass number, then atomic number, then length
    const auto d1_length = std::ssize(d1.z);
    const auto d2_length = std::ssize(d2.z);
    // -1 to ignore last item, which keeps bit-identical results as before when final daughter nuclide was not included
    // TODO: it would probably be better to sort by all items in reverse order
    const auto smallestpathlength = std::min(d1_length, d2_length) - 1;
    for (auto i = 0Z; i < smallestpathlength; i++) {
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

  for (auto& decaypath : localdecaypaths) {
    // all nuclei in the path (except for the last one, which is allowed to be stable) must have a mean life >0
    assert_always(std::all_of(decaypath.nucindex.cbegin(), decaypath.nucindex.cend() - 1,
                              [](const auto nucindex) { return get_meanlife(nucindex) > 0.; }));

    assert_always(decaypath.decaytypes.back() == DecayType::DECAYTYPE_NONE);

    // convert mean lifetimes to decay constants
    decaypath.lambdas.resize(decaypath.nucindex.size(), -1.);
    std::ranges::transform(decaypath.nucindex, decaypath.lambdas.begin(), [](const auto nucindex) {
      const double meanlife = get_meanlife(nucindex);
      // last nuclide might be stable (meanlife <= 0.)
      const double lambda = (meanlife > 0.) ? 1. / meanlife : 0.;
      return lambda;
    });
  }
  localdecaypaths.shrink_to_fit();

  return localdecaypaths;
}

// remove nuclides that are not a standard or custom input-specified nuclide, or connected to these by decays
void filter_unused_nuclides(const std::span<const int> custom_zlist, const std::span<const int> custom_alist,
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

    return !in_any_decaypath;
  });
  nuclides.shrink_to_fit();

  // update the nuclide indices in the decay paths after we possibly removed some nuclides
  for (auto& decaypath : decaypaths) {
    std::ranges::transform(decaypath.z, decaypath.a, decaypath.nucindex.begin(),
                           [](const auto z, const auto a) { return get_nucindex(z, a); });
  }
}

auto sample_decaytime(const int decaypathindex, const double tdecaymin, const double tdecaymax, rngstate_type& rngstate)
    -> double {
  double tdecay = -1;
  const double t_model = grid::get_t_model();
  // rejection method. draw random times with the right distribution until they are within the correct range.
  while (tdecay <= tdecaymin || tdecay >= tdecaymax) {
    tdecay = t_model;  // can't decay before initial model snapshot time
    const auto maxlength = std::ssize(decaypaths[decaypathindex].nucindex) - 1;
    for (auto i = 0Z; i < maxlength; i++) {
      tdecay -= get_meanlife(decaypaths[decaypathindex].nucindex[i]) *
                std::log(static_cast<double>(rng_uniform_pos(rngstate)));
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
// lambdas:            array of 1/(mean lifetime) for nuc[0]..nuc[num_nuclides-1]  [seconds^-1]
// useexpansionfactor: if true, return a modified 'abundance' at the end of the chain, with a weighting factor
//                          accounting for adiabatic loss from expansion since the decays occurred
//                          (This is needed to get the initial temperature)
constexpr auto calculate_decaychain(const double firstinitabund, const std::span<const double> lambdas,
                                    const double timediff, const bool useexpansionfactor) -> double {
  const int num_nuclides = static_cast<int>(lambdas.size());
  assert_testmodeonly(num_nuclides >= 1);

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

    // the Bateman solution is singular when two nuclides in a chain have equal decay constants
    assert_always(std::abs(denominator) > 0.);
    if (!useexpansionfactor) {
      // get abundance output
      sum += exp(-lambda_j * timediff) / denominator;
    } else {
      if (lambda_j > 0.) {
        const double sumtermtop =
            ((1 + (1 / lambda_j / timediff)) * exp(-timediff * lambda_j)) - (1. / lambda_j / timediff);
        sum += sumtermtop / denominator;
      }
    }
  }

  const double lastabund = firstinitabund * lambdaproduct * sum;
  assert_always(std::isfinite(lastabund));
  return lastabund;
}

// Get the mass fraction of a nuclide accounting for all decays and initial abundances.
// e.g., Co56 abundance may first increase with time due to Ni56 decays, then decrease due to Co56 decay
auto get_nuc_massfrac(const int nonemptymgi, const int nucindex, const double time) -> double {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  assert_always(time >= 0.);

  const double t_afterinit = time - grid::get_t_model();
  const auto [z, a] = get_nuc_z_a(nucindex);

  double nuctotal = 0.;  // abundance or decay rate, depending on mode parameter
  for (const auto& decaypath : decaypaths) {
    const auto last_decaytype = decaypath.decaytypes[decaypath.nucindex.size() - 2];
    // match 4He abundance to alpha decay of any nucleus (no continue), otherwise check daughter nuclide matches
    if ((z != 2 || a != 4 || last_decaytype != DecayType::DECAYTYPE_ALPHA) &&
        (decaypath.z.back() != z || decaypath.a.back() != a)) {
      continue;
    }

    const int nucindex_top = decaypath.nucindex[0];

    const double top_initabund = grid::get_modelinitnucmassfrac(modelgridindex, nucindex_top) / nucmass(nucindex_top);
    if (top_initabund <= 0.) {
      continue;
    }

    auto lambdas = decaypath.lambdas;
    if (z == 2 && a == 4) {
      // treat the end nuclide as stable He4
      lambdas[lambdas.size() - 1] = 0.;
    }

    const double massfraccontrib =
        (decaypath.branchproduct * calculate_decaychain(top_initabund, lambdas, t_afterinit, false) *
         nucmass(nucindex));
    nuctotal += massfraccontrib;
  }

  const auto meanlife = get_meanlife(nucindex);
  const auto lambda = (meanlife > 0.) ? 1. / meanlife : 0.;
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

  const int nucindex_top = decaypath.nucindex[0];

  const double top_initabund = grid::get_modelinitnucmassfrac(modelgridindex, nucindex_top) / nucmass(nucindex_top);
  if (top_initabund <= 0.) {
    return 0.;
  }

  const double t_afterinit = time - grid::get_t_model();

  // count the number of chain-top nuclei that haven't decayed past the end of the chain
  auto lambdas = decaypath.lambdas;
  // treat the end nuclide as stable to count how many got produced
  lambdas[lambdas.size() - 1] = 0.;

  const double abund_endsink = calculate_decaychain(top_initabund, lambdas, t_afterinit, false);
  const double ndecays_remaining = decaypath.branchproduct * (top_initabund - abund_endsink);
  // TODO ensure non-negative due to numerical precision?

  const double endecay = ndecays_remaining * get_decaypath_lastdecayenergy(decaypath);

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
  // min steps across the meanlifetime. Use a wide type: for short-lived nuclides the product easily
  // exceeds INT_MAX, which would overflow (undefined behaviour) if cast to int.
  const auto nsteps = static_cast<ptrdiff_t>(std::ceil((tstart - grid::get_t_model()) / min_meanlife) * 100000);

  double chain_endecay = 0.;
  double last_chain_endecay = -1.;
  double last_t = -1.;
  for (ptrdiff_t i = 0; i < nsteps; i++) {
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
  const double endiff = energy_tlow - energy_thigh;
  assert_always(std::isfinite(endiff));
  if (endiff < 0.) {
    // if the error is larger than just roundoff, this is a problem
    assert_always((energy_tlow * (1 + 2e-5)) >= energy_thigh);
    return 0.;
  }
  return endiff;
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

  const double top_initabund = grid::get_modelinitnucmassfrac(modelgridindex, nucindex_top) / nucmass(nucindex_top);
  if (top_initabund <= 0.) {
    return 0.;
  }

  const double t_afterinit = time - grid::get_t_model();

  const auto lambdas = std::span{decaypath.lambdas}.first(decaypath.lambdas.size() - 1);  // exclude the decay daughter

  const double endecay = get_decaypath_lastdecayenergy(decaypath);
  // contribution to the end nuclide abundance from the top of chain (could be a length-one chain
  // Z,A_top = Z,A_end so contribution would be from init abundance only)
  const double decayingnucabund =
      decaypath.branchproduct * calculate_decaychain(top_initabund, lambdas, t_afterinit, false);

  const int lastdecay_nucindex = decaypath.nucindex[decaypath.nucindex.size() - 2];

  const double decaypower = endecay * decayingnucabund / get_meanlife(lastdecay_nucindex);

  assert_always(decaypower >= 0.);
  assert_always(std::isfinite(decaypower));

  return decaypower;
}

auto write_nuclides_list() {
  auto nuclides_file = fstream_required("nuclides.out", std::ofstream::out | std::ofstream::trunc);
  std::println(nuclides_file, "#nucindex Z A");
  for (int nucindex = 0; nucindex < std::ssize(nuclides); nucindex++) {
    std::println(nuclides_file, "{} {} {}", nucindex, get_nuc_z(nucindex), get_nuc_a(nucindex));
  }
}

}  // anonymous namespace

[[nodiscard]] auto get_decay_neutrino_frac(const int nucindex, const DecayType decaytype) -> double {
  // subtract fraction of other decay products, nucdecayenergy() excludes neutrinos!
  const double nu_frac = 1. - (nucdecayenergy(nucindex, decaytype) / nuclides[nucindex].endecay_q[decaytype]);
  return nu_frac;
}

[[gnu::pure]] [[nodiscard]] auto get_num_nuclides() -> ptrdiff_t { return std::ssize(nuclides); }

[[nodiscard]] auto get_elname(const int z) -> std::string {
  assert_testmodeonly(z <= Z_MAX);
  return elsymbols.at(z);
}

[[nodiscard]] auto get_nuc_z(const int nucindex) -> int {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < std::ssize(nuclides));
  return nuclides[nucindex].z;
}

[[nodiscard]] auto get_nuc_a(const int nucindex) -> int {
  assert_testmodeonly(nucindex >= 0);
  assert_testmodeonly(nucindex < std::ssize(nuclides));
  return nuclides[nucindex].a;
}

// get the nuclide array index from the atomic number and mass number
[[nodiscard]] auto get_nucindex(const int z, const int a) -> int {
  const int nucindex = get_nucindex_or_neg_one(z, a);
  if (nucindex >= 0) {
    return nucindex;
  }
  printlnlog("[error] nuclide Z={} A={} not found in nuclide list", z, a);
  assert_always(false);  // nuclide not found
  return -1;
}

// check if nuclide exists in the simulation
[[nodiscard]] auto nuc_exists(const int z, const int a) -> bool { return get_nucindex_or_neg_one(z, a) >= 0; }

// average energy per decay in the form of gamma rays [erg]
[[nodiscard]] auto nucdecayenergygamma(const int nucindex) -> double { return nuclides[nucindex].endecay_gamma; }

[[nodiscard]] auto nucdecayenergygamma(const int z, const int a) -> double {
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
  auto i = 0ZU;
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

namespace {

// add the standard set of nuclides for the classical double-decay chains (56Ni, 57Ni, 48Cr,
// and 52Fe parents) with hardcoded decay data
void add_standard_nuclides() {
  // Ni57
  nuclides.push_back({.z = 28, .a = 57, .meanlife = 51.36 * HOUR});
  nuclides.back().endecay_positron = 0.354 * MEV;
  nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 0.436;
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1. - 0.436;
  const double ni57_q_ec = 3.261697 * MEV;
  nuclides.back().endecay_q[DECAYTYPE_ELECTRONCAPTURE] = ni57_q_ec;
  nuclides.back().endecay_q[DECAYTYPE_BETAPLUS] = ni57_q_ec - (2. * ME * CLIGHTSQUARED);

  // Ni56
  nuclides.push_back({.z = 28, .a = 56, .meanlife = 8.80 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;
  nuclides.back().endecay_q[DECAYTYPE_ELECTRONCAPTURE] = 2.132869 * MEV;

  // Co56
  nuclides.push_back({.z = 27, .a = 56, .meanlife = 113.7 * DAY});
  nuclides.back().endecay_positron = 0.63 * MEV;
  nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 0.19;
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1 - 0.19;
  const double co56_q_ec = 4.566645 * MEV;
  nuclides.back().endecay_q[DECAYTYPE_ELECTRONCAPTURE] = co56_q_ec;
  nuclides.back().endecay_q[DECAYTYPE_BETAPLUS] = co56_q_ec - (2. * ME * CLIGHTSQUARED);

  // Cr48
  nuclides.push_back({.z = 24, .a = 48, .meanlife = 1.29602 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;
  nuclides.back().endecay_q[DECAYTYPE_ELECTRONCAPTURE] = 1.656692 * MEV;

  // V48
  nuclides.push_back({.z = 23, .a = 48, .meanlife = 23.0442 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_BETAPLUS] = 0.499;
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1. - 0.499;
  nuclides.back().endecay_positron = 0.290 * MEV;
  const double v48_q_ec = 4.014947 * MEV;
  nuclides.back().endecay_q[DECAYTYPE_ELECTRONCAPTURE] = v48_q_ec;
  nuclides.back().endecay_q[DECAYTYPE_BETAPLUS] = v48_q_ec - (2. * ME * CLIGHTSQUARED);

  // Co57
  nuclides.push_back({.z = 27, .a = 57, .meanlife = 392.03 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;
  nuclides.back().endecay_q[DECAYTYPE_ELECTRONCAPTURE] = 0.836359 * MEV;

  // Fe52
  nuclides.push_back({.z = 26, .a = 52, .meanlife = 0.497429 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;
  nuclides.back().endecay_q[DECAYTYPE_ELECTRONCAPTURE] = 2.001543 * MEV;

  // Mn52
  nuclides.push_back({.z = 25, .a = 52, .meanlife = 0.0211395 * DAY});
  nuclides.back().branchprobs[DECAYTYPE_ELECTRONCAPTURE] = 1.;
  nuclides.back().endecay_q[DECAYTYPE_ELECTRONCAPTURE] = 5.085870 * MEV;
}

// add the nuclides with beta-minus decay data from betaminusdecays.txt
void read_betaminus_decaydata() {
  auto fbetaminus = fstream_required("betaminusdecays.txt", std::ios::in);
  std::string line;
  while (get_noncommentline(fbetaminus, line)) {
    // energies are average per beta decay
    // columns: # A, Z, Q[MeV], E_gamma[MeV], E_elec[MeV], E_neutrino[MeV], meanlife[s]
    int a = -1;
    int z = -1;
    double Q_betadecay_mev = 0.;
    double e_gamma_mev = 0.;
    double e_elec_mev = 0.;
    double e_neutrino = 0.;
    double tau_sec = 0.;
    assert_always(std::stringstream(line) >> a >> z >> Q_betadecay_mev >> e_gamma_mev >> e_elec_mev >> e_neutrino >>
                  tau_sec);
    if (Q_betadecay_mev > 0.) {
      assert_always(!nuc_exists(z, a));
      nuclides.push_back({.z = z, .a = a, .meanlife = tau_sec});
      nuclides.back().branchprobs[DECAYTYPE_BETAMINUS] = 1.;
      nuclides.back().endecay_q[DECAYTYPE_BETAMINUS] = Q_betadecay_mev * MEV;
      nuclides.back().endecay_electron = e_elec_mev * MEV;
      nuclides.back().endecay_gamma = e_gamma_mev * MEV;
      assert_always(e_elec_mev >= 0.);
    }
  }
}

// add or update the nuclides with alpha decay data from alphadecays.txt (also ensures that He4
// exists as a decay product)
void read_alpha_decaydata() {
  auto falpha = fstream_required("alphadecays.txt", std::ios::in);
  std::string line;
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
    double Q_alphadecay_mev = 0.;
    double Q_betadecay_mev = 0.;
    double e_alpha_mev = 0.;
    double e_gamma_mev = 0.;
    double e_beta_mev = 0.;
    assert_always(std::stringstream(line) >> a >> z >> branch_alpha >> branch_beta >> halflife >> Q_alphadecay_mev >>
                  Q_betadecay_mev >> e_alpha_mev >> e_gamma_mev >> e_beta_mev);

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
      nuclides[alphanucindex].endecay_q[DECAYTYPE_BETAMINUS] = Q_betadecay_mev * MEV;
      nuclides[alphanucindex].branchprobs[DECAYTYPE_ALPHA] = branch_alpha;
      nuclides[alphanucindex].endecay_q[DECAYTYPE_ALPHA] = Q_alphadecay_mev * MEV;
    }
  }
}

// add the nuclides with spontaneous fission decay data from fissiondecays.txt
void read_spontfission_decaydata() {
  auto ffission = fstream_required("fissiondecays.txt", std::ios::in);
  std::string line;
  while (get_noncommentline(ffission, line)) {
    int z_in = -1;
    int a_in = -1;
    double q_fission_mev = 0.;
    double e_gamma_mev = 0.;
    double e_1_mev = 0.;
    double e_2_mev = 0.;
    double m1 = 0.;
    double m2 = 0.;
    double z1 = 0.;
    double z2 = 0.;
    double tau_sec = 0.;
    assert_always(std::stringstream(line) >> a_in >> z_in >> q_fission_mev >> e_gamma_mev >> e_1_mev >> e_2_mev >> m1 >>
                  m2 >> z1 >> z2 >> tau_sec);
    assert_always(!nuc_exists(z_in, a_in));
    nuclides.push_back({.z = z_in, .a = a_in, .meanlife = tau_sec});
    nuclides.back().branchprobs[DECAYTYPE_SPONTFISSION] = 1.;
    nuclides.back().endecay_q[DECAYTYPE_SPONTFISSION] = q_fission_mev * MEV;
    nuclides.back().endecay_fission = q_fission_mev * MEV;  // will be overwritten if we have fission product data
    printlnlog("  added spontaneous fission nuclide: (Z={}){}{} meanlife {:.1e} days", z_in, get_elname(z_in), a_in,
               tau_sec / DAY);
  }
}

// read the fission product tables from fissionproducts_GEF_100keV.txt for the spontaneous
// fission nuclides that are in use
void read_fissionproduct_data() {
  auto ffission_products = fstream_required("fissionproducts_GEF_100keV.txt", std::ios::in);
  std::string line;
  while (get_noncommentline(ffission_products, line)) {
    int z_parent = -1;
    int a_parent = -1;
    assert_always(std::stringstream(line) >> z_parent >> a_parent);
    assert_always(get_noncommentline(ffission_products, line));
    double num_neutrons = 0;
    int tablesize = 0;
    double q_fission_mev = 0.;
    assert_always(std::stringstream(line) >> num_neutrons >> tablesize >> q_fission_mev);
    const int nucindex = get_nucindex_or_neg_one(z_parent, a_parent);
    const bool keep_table = (nucindex >= 0) && (nuclides[nucindex].branchprobs[DECAYTYPE_SPONTFISSION] > 0.);
    if (keep_table) {
      nuclides[nucindex].endecay_q[DECAYTYPE_SPONTFISSION] = q_fission_mev * MEV;
      nuclides[nucindex].endecay_fission = q_fission_mev * MEV;
      nuclides[nucindex].fission_daughters_z_a_prob.clear();
      nuclides[nucindex].fission_daughters_z_a_prob.reserve(tablesize);
    }

    double daughter_prob_sum = 0.;
    for (int i = 0; i < tablesize; i++) {
      assert_always(get_noncommentline(ffission_products, line));
      if (keep_table) {
        int daughter_a = -1;
        int daughter_z = -1;
        double probability_before_neutron_emission = 0.;
        double probability = 0.;
        assert_always(std::stringstream(line) >> daughter_a >> daughter_z >> probability_before_neutron_emission >>
                      probability);
        nuclides[nucindex].fission_daughters_z_a_prob.push_back(
            {.z = daughter_z, .a = daughter_a, .probability = probability});
        daughter_prob_sum += probability;
      }
    }
    if (keep_table) {
      nuclides[nucindex].decay_daughters_probsum = daughter_prob_sum;
    }
  }
}

// make sure that the input-specified nuclides and all decay daughters exist in the nuclide
// list, adding any missing ones as stable nuclides
void add_custom_and_daughter_nuclides(const std::span<const int> custom_zlist,
                                      const std::span<const int> custom_alist) {
  std::set<std::tuple<int, int>> nuclides_ensure_list{};

  // add any extra nuclides that were specified but not in the decay data files
  for (const auto [z, a] : std::views::zip(custom_zlist, custom_alist)) {
    nuclides_ensure_list.insert({z, a});
  }

  // ensure any daughters nuclides are included
  for (const auto& nuc : nuclides) {
    for (const auto& decaytype : all_decaytypes) {
      if (nuc.branchprobs[decaytype] > 0.) {
        for (const auto& daughter : decay_daughters_z_a_prob(nuc.z, nuc.a, decaytype)) {
          nuclides_ensure_list.insert({daughter.z, daughter.a});
        }
      }
    }
  }
  // add any required nuclides (assume they are stable)
  for (const auto& [z, a] : nuclides_ensure_list) {
    if (!nuc_exists(z, a)) {
      nuclides.push_back({.z = z, .a = a, .meanlife = -1});
    }
  }
}

// consistency checks on the nuclide list: no duplicate nuclides, branching probabilities that
// sum to one for unstable nuclides, and decay energies present for every active decay type
void check_nuclide_data() {
  std::set<std::tuple<int, int>> seen_z_a{};
  for (const auto& nuc : nuclides) {
    assert_always(!seen_z_a.contains({nuc.z, nuc.a}));  // duplicate nuclide detected
    seen_z_a.insert({nuc.z, nuc.a});

    const auto branchprob_sum = std::ranges::fold_left(nuc.branchprobs, 0.0, std::plus{});
    if (branchprob_sum > 0.) {
      assert_always(nuc.meanlife > 0.);  // unstable nuclide has decay modes
      assert_always(std::abs(branchprob_sum - 1.) < 1e-3);  // branching ratios sum to 100%
    } else {
      assert_always(nuc.meanlife < 0.);  // stable nuclide has no decay modes
    }

    if (nuc.branchprobs[DECAYTYPE_BETAPLUS] > 0.) {
      assert_always(nuc.endecay_positron >= 0.);
      assert_always(nuc.endecay_q[DECAYTYPE_BETAPLUS] >= 0.);
    }
    if (nuc.branchprobs[DECAYTYPE_BETAMINUS] > 0.) {
      assert_always(nuc.endecay_electron >= 0.);
      assert_always(nuc.endecay_q[DECAYTYPE_BETAMINUS] >= 0.);
    }
    if (nuc.branchprobs[DECAYTYPE_ALPHA] > 0.) {
      assert_always(nuc.endecay_alpha >= 0.);
      assert_always(nuc.endecay_q[DECAYTYPE_ALPHA] >= 0.);
    }
    if (nuc.branchprobs[DECAYTYPE_SPONTFISSION] > 0.) {
      assert_always(nuc.endecay_fission >= 0.);
      assert_always(nuc.endecay_q[DECAYTYPE_SPONTFISSION] >= 0.);
    }
    assert_always(nuc.endecay_gamma >= 0.);
  }
}

}  // anonymous namespace

// add all nuclides and decays, and later trim any irrelevant ones (not connected to input-specified nuclei)
void init_nuclides(const std::span<const int> custom_zlist, const std::span<const int> custom_alist) {
  assert_always(custom_zlist.size() == custom_alist.size());

  add_standard_nuclides();

  // deliberately a copy, not a reference: the decay-data readers below add further nuclides
  // cppcheck-suppress redundantCopyLocalConst
  const auto standard_nuclides = nuclides;

  read_betaminus_decaydata();
  read_alpha_decaydata();
  read_spontfission_decaydata();
  read_fissionproduct_data();

  add_custom_and_daughter_nuclides(custom_zlist, custom_alist);

  check_nuclide_data();

  printlnlog("Number of nuclides before filtering: {}", std::ssize(nuclides));
  decaypaths = find_decaypaths(custom_zlist, custom_alist, standard_nuclides);
  filter_unused_nuclides(custom_zlist, custom_alist, standard_nuclides);

  printlnlog("Number of nuclides: {}", std::ssize(nuclides));

  const int maxdecaypathlength = std::ranges::fold_left(decaypaths, 0ZU, [](const auto maxlen, const auto& decaypath) {
    return std::max(maxlen, decaypath.nucindex.size());
  });

  printlnlog("Number of decay paths: {} (max length {})", get_num_decaypaths(), maxdecaypathlength);
  constexpr bool print_decaypaths = false;
  if (print_decaypaths) {
    for (int decaypathindex = 0; decaypathindex < get_num_decaypaths(); decaypathindex++) {
      printout_decaypath(decaypathindex);
    }
  }

  reserve_resize(alldecaytypes_is_used, DecayType::DECAYTYPE_COUNT);
  for (DecayType decaytype : all_decaytypes) {
    alldecaytypes_is_used[decaytype] = std::ranges::any_of(
        std::views::iota(0UZ, nuclides.size()),
        [decaytype](const auto nucindex) { return get_nuc_decaybranchprob(nucindex, decaytype) > 0.; });
  }

  // Read in data for gamma ray lines and make a list of them in energy order.
  gammapkt::init_gamma_data();

  // TODO: generalise this to all included nuclides
  printlnlog("decayenergy(NI56) {:g} MeV, decayenergy(CO56) {:g} MeV, decayenergy_gamma(CO56) {:g} MeV",
             nucdecayenergytotal(28, 56) / MEV, nucdecayenergytotal(27, 56) / MEV, nucdecayenergygamma(27, 56) / MEV);
  printlnlog("decayenergy(NI57) {:g} MeV, decayenergy_gamma(NI57) {:g} MeV, decayenergy(CO57) {:g} MeV",
             nucdecayenergytotal(28, 57) / MEV, nucdecayenergygamma(28, 57) / MEV, nucdecayenergytotal(27, 57) / MEV);
  printlnlog("decayenergy(CR48) {:g} MeV, decayenergy(V48) {:g} MeV", nucdecayenergytotal(24, 48) / MEV,
             nucdecayenergytotal(23, 48) / MEV);
  printlnlog("decayenergy(FE52) {:g} MeV, decayenergy(MN52) {:g} MeV", nucdecayenergytotal(26, 52) / MEV,
             nucdecayenergytotal(25, 52) / MEV);

  if (globals::my_rank == 0 && !globals::simulation_continued_from_saved) {
    write_nuclides_list();
  }
}

[[nodiscard]] auto decaytype_is_used(const DecayType decaytype) -> bool {
  assert_testmodeonly(!all_decaytypes.empty());
  assert_testmodeonly(decaytype >= 0);
  assert_testmodeonly(decaytype < DecayType::DECAYTYPE_COUNT);
  return alldecaytypes_is_used[decaytype];
}

// calculate the decay energy per unit mass [erg/g] released from time t_model (can be before tmin) to tstart,
// accounting for the photon energy loss due to expansion between time of decays and tstart (equation 18 of Lucy 2005)
auto get_endecay_per_ejectamass_tmodel_to_time_withexpansion(const int nonemptymgi, const double tstart) -> double {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double t_afterinit = tstart - grid::get_t_model();
  double tot_endecay = 0.;
  for (const auto& decaypath : decaypaths) {
    const int nucindex_top = decaypath.nucindex[0];

    const double top_initabund = grid::get_modelinitnucmassfrac(modelgridindex, nucindex_top) / nucmass(nucindex_top);
    auto lambdas = decaypath.lambdas;
    // treat the end nuclide as stable to count how many got produced by the chain
    lambdas[lambdas.size() - 1] = 0.;

    const double chain_endecay =
        (decaypath.branchproduct * calculate_decaychain(top_initabund, lambdas, t_afterinit, true) *
         get_decaypath_lastdecayenergy(decaypath));

    tot_endecay += chain_endecay;
  }

  return tot_endecay;
}

// get the decay energy that will be released during the simulation time [(tmodel if initial packets else tmin) to tmax]
// indexed by nonemptymgi and decaypathindex [erg/g]
auto get_modelcell_simtime_endecay_per_mass(const int nonemptymgi,
                                            const std::span<const double> energy_per_mass_nonemptymgi_decaypath)
    -> double {
  const auto num_decaypaths = get_num_decaypaths();
  double endecay_per_mass = 0.;
  for (int decaypathindex = 0; decaypathindex < num_decaypaths; decaypathindex++) {
    endecay_per_mass += energy_per_mass_nonemptymgi_decaypath[(nonemptymgi * num_decaypaths) + decaypathindex];
  }
  return endecay_per_mass;
}

// energy_per_mass_nonemptymgi_decaypath is an array indexed by [nonemptymgi * num_decaypaths + i] will hold the
// decay energy per mass [erg/g] released by chain i in cell mgi during the simulation time range tmin to tmax
auto get_energy_per_mass_nonemptymgi_decaypath() -> MPI_shared_array<const double> {
  const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();
  printlog(
      "[info] mem_usage: energy_per_mass_nonemptymgi_decaypath[nonempty_npts_model*num_decaypaths] occupies {:.1f} "
      "MB (node shared memory)...",
      nonempty_npts_model * get_num_decaypaths() * sizeof(double) / 1024. / 1024.);
  auto energy_per_mass_nonemptymgi_decaypath = MPI_shared_array<double>{nonempty_npts_model * get_num_decaypaths(), 0.};
  printlnlog("done.");

  MPI_Barrier_allranks();
  const auto time_min_decay = INITIAL_PACKETS_ON ? grid::get_t_model() : globals::tmin;
  const ptrdiff_t num_decaypaths = get_num_decaypaths();
  for (int nonemptymgi = 0; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
    if (nonemptymgi % globals::node_nprocs == globals::rank_in_node) {
      const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
      for (int decaypathindex = 0; decaypathindex < num_decaypaths; decaypathindex++) {
        const double energy_per_mass =
            get_endecay_per_ejectamass_between_times(mgi, decaypathindex, time_min_decay, globals::tmax);
        assert_testmodeonly(energy_per_mass >= 0.);
        assert_testmodeonly(std::isfinite(energy_per_mass));
        energy_per_mass_nonemptymgi_decaypath[(nonemptymgi * num_decaypaths) + decaypathindex] = energy_per_mass;
      }
    }
  }

  MPI_Barrier_allranks();
  return energy_per_mass_nonemptymgi_decaypath;
}

// energy release rate in form of kinetic energy of positrons, electrons, and alpha particles in [erg/s/g]
[[nodiscard]] auto get_particle_injection_rate(const int nonemptymgi, const double t, const DecayType decaytype)
    -> double {
  double dep_sum = 0.;
  const auto num_nuclides = std::ssize(nuclides);
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
  const auto num_nuclides = std::ssize(nuclides);
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
[[nodiscard]] auto get_qdot_modelcell(const int nonemptymgi, const double t, const DecayType decaytype) -> double {
  double qdot = 0.;
  const auto num_nuclides = std::ssize(nuclides);
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
                  get_decaypath_lastdecayenergy(decaypath));
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
    const auto num_nuclides = std::ssize(nuclides);
    for (int nucindex = 0; nucindex < num_nuclides; nucindex++) {
      if (get_nuc_z(nucindex) == atomic_number) {
        const double nuc_massfrac = get_nuc_massfrac(nonemptymgi, nucindex, t_current);
        isomassfracsum += nuc_massfrac;
        isomassfrac_on_nucmass_sum += nuc_massfrac / nucmass(nucindex);
      }
    }

    const double otherstablemassfrac = grid::get_elem_untrackedstable_initmassfrac(nonemptymgi, element);
    isomassfracsum += otherstablemassfrac;
    isomassfrac_on_nucmass_sum += otherstablemassfrac / globals::elements[element].initstablemeannucmass;

    grid::set_elem_massfrac(nonemptymgi, element, static_cast<float>(isomassfracsum));
    const auto meanweight = static_cast<float>(isomassfracsum / isomassfrac_on_nucmass_sum);
    grid::set_element_meanweight(
        nonemptymgi, element,
        (std::isfinite(meanweight) && meanweight > 0.) ? meanweight : globals::elements[element].initstablemeannucmass);
  }

  grid::set_nnetot(nonemptymgi);
}

void output_isotopic_densities(std::ostream& estimators_file, const int nonemptymgi, const double t_current,
                               const int element) {
  const double rho = grid::get_rho(nonemptymgi);

  const int atomic_number = get_atomicnumber(element);
  std::set<std::tuple<int, int>> a_isotopes;  // so that we output sorted by mass number
  const auto num_nuclides = std::ssize(nuclides);
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
      std::print(estimators_file, "  {}{}: {:9.3e}", get_elname(atomic_number), nuc_a, numberdens);
    }
  }

  const double otherstablemassfrac = grid::get_elem_untrackedstable_initmassfrac(nonemptymgi, element);
  if (otherstablemassfrac > 0) {
    // factor to convert convert mass fraction to number density
    const double meannucmass = globals::elements[element].initstablemeannucmass;
    const double otherstable_numberdens = otherstablemassfrac / meannucmass * grid::get_rho(nonemptymgi);
    std::print(estimators_file, "  {}_otherstable: {:9.3e}", get_elname(atomic_number), otherstable_numberdens);
  }
  std::println(estimators_file, "");
}

void setup_radioactive_pellet(const double e_cmf_per_packet, const int nonemptymgi, Packet& pkt,
                              const std::span<const double> energy_per_mass_nonemptymgi_decaypath) {
  const int num_decaypaths = get_num_decaypaths();

  // decay channels include all radioactive decay paths, and possibly also an initial cell energy channel
  const int num_decaychannels = num_decaypaths + ((INITIAL_PACKETS_ON && USE_MODEL_INITIAL_ENERGY) ? 1 : 0);

  auto cumulative_en_sum = std::vector<double>(num_decaychannels, 0.);
  double energysum = 0.;

  // add the radioactive decay paths
  for (int decaypathindex = 0; decaypathindex < num_decaypaths; decaypathindex++) {
    energysum += energy_per_mass_nonemptymgi_decaypath[(nonemptymgi * num_decaypaths) + decaypathindex];
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

  const double zrand_en = rng_uniform(get_rngstate(pkt)) * cumulative_en_sum[num_decaychannels - 1];

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
    pkt.tdecay = sample_decaytime(decaypathindex, tdecaymin, globals::tmax, get_rngstate(pkt));
    pkt.e_cmf = e_cmf_per_packet;
  } else {
    // use uniform decay time distribution and scale the packet energies instead.
    // keeping the pellet decay rate constant will give better statistics at late times
    // when very little energy and few packets are released
    // NB: the interpolation endpoints are deliberately (tmax, tdecaymin) rather than the other way round.
    // The draw is uniform over [tdecaymin, tmax] either way; swapping them would change which decay time
    // each random number maps to, and hence every result.
    pkt.tdecay = std::lerp(globals::tmax, tdecaymin, rng_uniform(get_rngstate(pkt)));

    // we need to scale the packet energy up or down according to decay rate at the randomly selected time.
    // e_cmf_average is the average energy per packet for this cell and decaypath, so we scale this up or down
    // according to: decay power at this time relative to the average decay power
    const double avgpower = energy_per_mass_nonemptymgi_decaypath[(nonemptymgi * num_decaypaths) + decaypathindex] /
                            (globals::tmax - tdecaymin);
    assert_always(avgpower > 0.);
    pkt.e_cmf =
        e_cmf_per_packet * get_decaypath_power_per_ejectamass(decaypathindex, nonemptymgi, pkt.tdecay) / avgpower;
    assert_always(pkt.e_cmf >= 0);
  }

  pkt.type = TYPE_RADIOACTIVE_PELLET;

  // final decaying nuclide in the chain (one before the end, which is the daughter of the last decay)
  const auto pathlength = decaypaths[decaypathindex].nucindex.size();
  pkt.pellet_nucindex = decaypaths[decaypathindex].nucindex[pathlength - 2];
  pkt.pellet_decaytype = decaypaths[decaypathindex].decaytypes[pathlength - 2];

  const auto engamma = nucdecayenergygamma(pkt.pellet_nucindex);
  const auto enparticle = nucdecayenergyparticle(pkt.pellet_nucindex, static_cast<DecayType>(pkt.pellet_decaytype));

  pkt.originated_from_particlenotgamma = (rng_uniform(get_rngstate(pkt)) >= engamma / (engamma + enparticle));
  if (pkt.originated_from_particlenotgamma) {
    // particle (positron, electron, or alpha) emitted
    pkt.nu_cmf = enparticle / H;
  } else {
    // gamma ray emitted
    pkt.nu_cmf = gammapkt::choose_gamma_ray(pkt.pellet_nucindex, get_rngstate(pkt));
  }
}

}  // namespace decay
