// Declarations for the radioactive decay module (decay.cc).

#ifndef DECAY_H
#define DECAY_H

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <span>
#include <string>

#include "mpi_logging.h"
#include "packet.h"

namespace decay {

enum DecayType : int {
  DECAYTYPE_ALPHA = 0,
  DECAYTYPE_ELECTRONCAPTURE = 1,
  DECAYTYPE_BETAPLUS = 2,
  DECAYTYPE_BETAMINUS = 3,
  DECAYTYPE_NONE = 4,
  DECAYTYPE_SPONTFISSION = 5,
  DECAYTYPE_COUNT = 6,
};

constexpr std::array all_decaytypes{
    DecayType::DECAYTYPE_ALPHA,     DecayType::DECAYTYPE_ELECTRONCAPTURE, DecayType::DECAYTYPE_BETAPLUS,
    DecayType::DECAYTYPE_BETAMINUS, DecayType::DECAYTYPE_SPONTFISSION,
};

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
  // The expansion-factor result is only the energy-weighted count of decays reaching the end of the
  // chain, and so only non-negative, if the chain end is treated as a sink. Left decaying, the same
  // sum is the weighted integral of a signed derivative and may legitimately come out negative,
  // which the clamp at the end would then swallow. The one caller zeroes it; require that.
  assert_testmodeonly(!useexpansionfactor || (lambdas.back() == 0.));

  double lambdaproduct = 1.;
  for (int j = 0; j < (num_nuclides - 1); j++) {
    lambdaproduct *= lambdas[j];
  }

  double sum = 0;
  double maxabsterm = 0.;
  for (int j = 0; j < num_nuclides; j++) {
    const auto lambda_j = lambdas[j];
    double denominator = 1.;
    for (int p = 0; p < num_nuclides; p++) {
      if (p != j) {
        denominator *= (lambdas[p] - lambda_j);
      }
    }

    // the Bateman solution is singular when two nuclides in a chain have equal decay constants
    // (init_nuclides() rejects decay data that could trigger this)
    assert_always(std::abs(denominator) > 0.);
    double sumterm = 0.;
    if (!useexpansionfactor) {
      // get abundance output
      sumterm = exp(-lambda_j * timediff) / denominator;
    } else {
      // (1 + 1/x) * exp(-x) - 1/x, rearranged as exp(-x) + expm1(-x)/x. Written the first way it
      // subtracts two numbers of size 1/x to leave a result of size x/2, so the relative error
      // grows as epsilon/x^2 and it underflows to exactly zero below x of about 1e-8. 237-Np at
      // 1e5 s has x = 1e-9, so every long-lived nuclide was contributing nothing at all here.
      //
      // x is zero for a stable nuclide, and also if timediff is zero, which would need the first
      // timestep midpoint to fall exactly on the model snapshot time. Either way the term is zero:
      // expm1(-x)/x tends to -1 as x tends to zero, so the bracket tends to zero. Both forms divide
      // by x, so leaving it unguarded would give 0/0 here and inf - inf in the original.
      const double x = lambda_j * timediff;
      if (x > 0.) {
        sumterm = (exp(-x) + (std::expm1(-x) / x)) / denominator;
      }
    }
    sum += sumterm;
    const auto abssumterm = std::abs(sumterm);
    maxabsterm = (abssumterm > maxabsterm) ? abssumterm : maxabsterm;
  }

  const double lastabund = firstinitabund * lambdaproduct * sum;
  assert_always(std::isfinite(lastabund));

  // Neither quantity returned here can be negative. An abundance is a number of nuclei, and the
  // expansion-factor result is the decay energy still available after adiabatic losses, which cannot
  // remove more than the decay released. So a negative only ever means the Bateman sum has lost its
  // significance, which happens readily now that the actinide chains are connected: the terms
  // alternate in sign with magnitudes scaling as the spread of the decay constants, and a chain from
  // 238-U to 206-Pb spans twenty-one orders of magnitude. Evaluated in 60-digit precision the 237-Np
  // chain at 1e5 s gives +3.5e-54 where double precision gives -5.8e-23.
  //
  // The magnitude test additionally covers timediff = 0, where the exact sum is zero for any chain of
  // two or more nuclides and the four supernova chains land on either sign of 1e-16 purely by accident
  // of rounding. Testing only the sign there would leave the answer depending on which way it fell,
  // which a bit-reproducible build cannot have.
  if ((lastabund < 0.) || (std::abs(sum) <= (num_nuclides * std::numeric_limits<double>::epsilon() * maxabsterm))) {
    return 0.;
  }

  return lastabund;
}

void init_nuclides(std::span<const int> custom_zlist, std::span<const int> custom_alist);
[[nodiscard]] auto decaytype_is_used(DecayType decaytype) -> bool;
[[nodiscard]] auto get_nucstring_z(const std::string& strnuc) -> int;
[[nodiscard]] auto get_nucstring_a(const std::string& strnuc) -> int;
[[gnu::pure]] [[nodiscard]] auto get_num_nuclides() -> ptrdiff_t;
[[nodiscard]] auto get_elname(int z) -> std::string;
[[nodiscard]] auto get_nuc_z(int nucindex) -> int;
[[nodiscard]] auto get_nuc_a(int nucindex) -> int;
[[nodiscard]] auto get_nucindex(int z, int a) -> int;
[[nodiscard]] auto nuc_exists(int z, int a) -> bool;
[[nodiscard]] auto nucdecayenergygamma(int nucindex) -> double;
[[nodiscard]] auto nucdecayenergygamma(int z, int a) -> double;
[[nodiscard]] auto get_decay_neutrino_frac(int nucindex, DecayType decaytype) -> double;
void set_nucdecayenergygamma(int nucindex, double value);
void update_abundances(int nonemptymgi, double t_current);
[[nodiscard]] auto get_endecay_per_ejectamass_tmodel_to_time_withexpansion(int nonemptymgi, double tstart) -> double;
[[nodiscard]] auto get_modelcell_simtime_endecay_per_mass(int nonemptymgi,
                                                          std::span<const double> energy_per_mass_nonemptymgi_decaypath)
    -> double;
auto get_energy_per_mass_nonemptymgi_decaypath() -> MPI_shared_array<const double>;
[[nodiscard]] auto get_qdot_modelcell(int nonemptymgi, double t, DecayType decaytype) -> double;
[[nodiscard]] auto get_particle_injection_rate(int nonemptymgi, double t, DecayType decaytype) -> double;
[[nodiscard]] auto get_gamma_emission_rate(int nonemptymgi, double t) -> double;
[[nodiscard]] auto get_global_etot_tmodel_tinf() -> double;
void output_isotopic_densities(std::ostream& estimators_file, int nonemptymgi, double t_current, int element);
// Construct an indivisible radioactive pellet by energy-weighted sampling a decay path or the optional initial-energy
// channel, then sample its release time and record the emitting nuclide and decay type.
// Lucy (2005), doi:10.1051/0004-6361:20041656.
void setup_radioactive_pellet(double e_cmf_per_packet, int nonemptymgi, Packet& pkt,
                              std::span<const double> energy_per_mass_nonemptymgi_decaypath);

}  // namespace decay

#endif  // DECAY_H
