// Declarations for r-packet propagation (rpkt.cc) plus inline helpers for line-list searches
// and the continuum opacity structure.

#ifndef RPKT_H
#define RPKT_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <span>

#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "mpi_logging.h"
#include "packet.h"
#include "random.h"

constexpr double expopac_lambdamin = 60.;
constexpr double expopac_lambdamax = 40000.;
constexpr double expopac_deltalambda = 20.;
constexpr auto expopac_nbins = static_cast<ptrdiff_t>((expopac_lambdamax - expopac_lambdamin) / expopac_deltalambda);

// wavelength bins are ordered by ascending wavelength (descending frequency)

constexpr auto get_expopac_bin_nu_upper(const ptrdiff_t binindex) -> double {
  const auto binindex_f64 = static_cast<double>(binindex);
  const auto lambda_lower = expopac_lambdamin + (binindex_f64 * expopac_deltalambda);
  return 1e8 * CLIGHT / lambda_lower;
}

constexpr auto get_expopac_bin_nu_lower(const ptrdiff_t binindex) -> double {
  const auto binindexplusone_f64 = static_cast<double>(binindex + 1);
  const auto lambda_upper = expopac_lambdamin + (binindexplusone_f64 * expopac_deltalambda);
  return 1e8 * CLIGHT / lambda_upper;
}

static_assert(get_expopac_bin_nu_lower(0) == get_expopac_bin_nu_upper(1));  // bins are contiguous
static_assert(get_expopac_bin_nu_lower(0) < get_expopac_bin_nu_upper(0));
static_assert(get_expopac_bin_nu_upper(expopac_nbins - 1) > get_expopac_bin_nu_lower(expopac_nbins - 1));

// kappa in cm^2/g for each bin of each non-empty cell
inline MPI_shared_array<float> expansionopacities{};

struct Phixslist {
  // for either USE_LUT_PHOTOION = true or USE_ION_BFHEATING_ESTIMATORS = true. Size =
  // nbfcontinua_ground
  std::span<double> groundcont_gamma_contr;
  // needed for DETAILED_BF_ESTIMATORS_ON. Size = bfestimcount
  std::span<double> gamma_contr;
  int bfestimend{-1};
  int bfestimbegin{0};

  // unique ptrs are used instead of vectors for nvc++ compatibility in device code
  // NOLINTBEGIN(*-avoid-c-arrays)
  std::unique_ptr<double[]> _groundcont_gamma_contr;
  std::unique_ptr<double[]> _gamma_contr;
  // NOLINTEND(*-avoid-c-arrays)
};

struct ContinuumOpacity {
  double nu{-1.};  // frequency at which opacity was calculated
  // Each chi is an extinction coefficient [cm^-1]; density is already included when it is calculated, so a
  // homogeneous path contributes an optical depth chi * pathlength.
  double chi_escatter{0.};  // Thomson scattering on free electrons (stay rpacket) contribution to the opacity
  double chi_freefree_heat{0.};  // free-free heating (become kpacket) contribution to the opacity
  double chi_boundfree{0.};  // bound-free (photoionization) contribution to the opacity
  int nonemptymgi{-1};
  int timestep{-1};
  Phixslist phixslist;

  constexpr explicit ContinuumOpacity(const bool alloc_phixslist) {
    if (alloc_phixslist) {
// NOLINTBEGIN(*-avoid-c-arrays)
#pragma clang unsafe_buffer_usage begin
      phixslist._groundcont_gamma_contr = std::make_unique<double[]>(globals::nbfcontinua_ground);
      phixslist.groundcont_gamma_contr =
          std::span<double>(phixslist._groundcont_gamma_contr.get(), globals::nbfcontinua_ground);
      const auto bfestimcount = globals::bfestim_nu_edge.size();
      phixslist._gamma_contr = std::make_unique<double[]>(bfestimcount);
      phixslist.gamma_contr = std::span<double>(phixslist._gamma_contr.get(), bfestimcount);
#pragma clang unsafe_buffer_usage end
      // NOLINTEND(*-avoid-c-arrays)
    }
  }

  // default constructor allocates phixslist
  constexpr ContinuumOpacity() : ContinuumOpacity(true) {}

  // total continuum absorption coefficient at nu [cm^-1]
  [[nodiscard]] constexpr auto total() const { return chi_escatter + chi_boundfree + chi_freefree_heat; }
};

DEVICE_FUNC void do_rpkt(Packet& pkt, double t2, ContinuumOpacity& chi_rpkt_cont);
DEVICE_FUNC void emit_rpkt(Packet& pkt);
template <bool USECELLHISTANDUPDATEPHIXSLIST>
void calculate_chi_rpkt_cont(double nu_cmf, ContinuumOpacity& chi_rpkt_cont, int nonemptymgi);
extern template void calculate_chi_rpkt_cont<true>(double nu_cmf, ContinuumOpacity& chi_rpkt_cont, int nonemptymgi);
extern template void calculate_chi_rpkt_cont<false>(double nu_cmf, ContinuumOpacity& chi_rpkt_cont, int nonemptymgi);
[[nodiscard]] DEVICE_FUNC auto sample_planck_times_expansion_opacity(int nonemptymgi, rngstate_type& rngstate)
    -> double;
void allocate_expansionopacities();
// Convert Sobolev line optical depths in each wavelength bin into an expansion mass opacity. When requested, also
// construct the Planck-weighted cumulative distribution used to sample thermal re-emission frequencies.
// Eastman & Pinto (1993), doi:10.1086/172957; Karp et al. (1977), doi:10.1086/155241.
void calculate_expansion_opacities(int nonemptymgi);
void MPI_Bcast_binned_opacities(ptrdiff_t nonemptymgi, int root_node_id);
auto calculate_chi_ffheat_nnionpart(int nonemptymgi) -> double;

[[nodiscard]] constexpr auto get_linedistance(const double prop_time, const double nu_cmf, const double nu_trans,
                                              const double dnu_on_dl) -> double {
  // distance from packet position to redshifting into line at frequency nu_trans

  if (nu_cmf <= nu_trans) {
    return 0.;  // photon was propagated too far, make sure that we don't miss a line
  }
  const double delta_nu = nu_cmf - nu_trans;  // positive number

  if constexpr (USE_RELATIVISTIC_DOPPLER_SHIFT) {
    // With special relativity, the Doppler shift formula has an extra factor of 1/gamma in it,
    // which changes the distance reach a line resonance and creates a dependence on packet position and direction

    // use linear interpolation of frequency along the path
    return -delta_nu / dnu_on_dl;  // dnu_on_dl is negative, so this is a positive distance
  }

  return CLIGHT * prop_time * delta_nu / nu_trans;
}

static_assert(get_linedistance(100., 1., 2., -0.5) == 0.);  // overshot the line resonance
static_assert(USE_RELATIVISTIC_DOPPLER_SHIFT || get_linedistance(2., 4., 2., -1.) == (CLIGHT * 2. * 2. / 2.));
static_assert(!USE_RELATIVISTIC_DOPPLER_SHIFT || get_linedistance(2., 4., 2., -1.) == 2.);

// find the next transition lineindex redder than nu_cmf
// for the propagation through non empty cells
// return -1 if no transition can be reached
[[gnu::pure]] [[nodiscard]] constexpr auto closest_transition(const double nu_cmf, const int next_trans,
                                                              const std::span<const double> linelistnu) -> int {
  const int nlines = static_cast<int>(linelistnu.size());
  if (next_trans > (nlines - 1)) {
    // packet is tagged as having no more line interactions
    return -1;
  }
  // if nu_cmf is smaller than the lowest frequency in the linelist,
  // no line interaction is possible: return negative value as a flag
  if (nu_cmf < linelistnu.back()) {
    return -1;
  }

  if (next_trans > 0) [[likely]] {
    // if next_trans > 0 we know the next line we should interact with, independent of the packets
    // current nu_cmf which might be smaller than globals::linelist[left].nu due to propagation errors
    return next_trans;
  }
  if (nu_cmf >= linelistnu[0]) {
    // if nu_cmf is larger than the highest frequency in the the linelist,
    // interaction with the first line occurs - no search
    return 0;
  }
  // otherwise go through the list until nu_cmf is located between two
  // entries in the line list and get the index of the closest line to lower frequencies

  // will find the highest frequency (lowest index) line with nu_line <= nu_cmf
  // lower_bound matches the first element where the comparison function is false
  const int matchindex =
      static_cast<int>(std::ranges::lower_bound(linelistnu, nu_cmf, std::ranges::greater{}) - linelistnu.begin());

  return matchindex;
}

// the linelist is sorted by descending frequency
inline constexpr std::array<double, 4> test_closest_transition_linelist_nu{9., 7., 5., 3.};
static_assert(closest_transition(10., -1, test_closest_transition_linelist_nu) == 0);  // bluer than the whole list
static_assert(closest_transition(8., -1, test_closest_transition_linelist_nu) == 1);  // binary search case
static_assert(closest_transition(5., -1, test_closest_transition_linelist_nu) == 2);  // exact frequency match
static_assert(closest_transition(2., -1, test_closest_transition_linelist_nu) == -1);  // redder than the whole list
static_assert(closest_transition(8., 2, test_closest_transition_linelist_nu) == 2);  // known next transition
static_assert(closest_transition(8., 4, test_closest_transition_linelist_nu) == -1);  // no more line interactions

// Whether a bound-free continuum contributes to the opacity of a cell, i.e. whether the cell contains enough
// of the involved atomic species. Continua that fail this are given a zero level population (see
// CellCache::allcont_nnlevel), which is the single condition that calculate_chi_bf_gammacontr() skips on.
[[gnu::pure]] [[nodiscard]] inline auto keep_this_cont(int element, const int ion, const int level,
                                                       const int nonemptymgi, const float nnetot) -> bool {
  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    return grid::get_elem_massfrac(nonemptymgi, element) > 0;
  }
  return ((get_nnion(nonemptymgi, element, ion) / nnetot > 1.e-6) || (level == 0));
}

#endif  // RPKT_H
