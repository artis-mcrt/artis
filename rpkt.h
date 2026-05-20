#ifndef RPKT_H
#define RPKT_H

#include <algorithm>
#include <cstddef>
#include <ctime>
#include <functional>
#include <memory>
#include <span>

#include "artisoptions.h"
#include "constants.h"
#include "grid.h"
#include "ltepop.h"
#include "mpi_logging.h"
#include "packet.h"

constexpr float expopac_lambdamin = 534.5;
constexpr float expopac_lambdamax = 35000.;
constexpr float expopac_deltalambda = 35.5;
constexpr auto expopac_nbins = static_cast<ptrdiff_t>((expopac_lambdamax - expopac_lambdamin) / expopac_deltalambda);

// wavelength bins are ordered by ascending wavelength (descending frequency)

constexpr auto get_expopac_bin_nu_upper(const ptrdiff_t binindex) -> double {
  const auto binindex_f32 = static_cast<decltype(expopac_deltalambda)>(binindex);
  const auto lambda_lower = expopac_lambdamin + (binindex_f32 * expopac_deltalambda);
  return 1e8 * CLIGHT / lambda_lower;
}

constexpr auto get_expopac_bin_nu_lower(const ptrdiff_t binindex) -> double {
  const auto binindexplusone_f32 = static_cast<decltype(expopac_deltalambda)>(binindex + 1);
  const auto lambda_upper = expopac_lambdamin + (binindexplusone_f32 * expopac_deltalambda);
  return 1e8 * CLIGHT / lambda_upper;
}

// kappa in cm^2/g for each bin of each non-empty cell
inline MPI_shared_array<float> expansionopacities{};

class Phixslist {
  // NOLINTBEGIN(*-avoid-c-arrays)
 public:
  // for either USE_LUT_PHOTOION = true or USE_ION_BFHEATING_ESTIMATORS = true. Size =
  // nbfcontinua_ground
  std::span<double> groundcont_gamma_contr;
  // cumulative sum of all bound-free continua absorption coefficients. Size = nbfcontinua
  std::span<double> chi_bf_sum;
  // needed for DETAILED_BF_ESTIMATORS_ON. Size = bfestimcount
  std::span<double> gamma_contr;
  int allcontend{-1};
  int allcontbegin{0};
  int bfestimend{-1};
  int bfestimbegin{0};

#pragma clang unsafe_buffer_usage begin
  constexpr Phixslist(const int nbfcontinua_ground, const int nbfcontinua, const int bfestimcount)
      : _groundcont_gamma_contr{std::make_unique<double[]>(nbfcontinua_ground)},
        _chi_bf_sum{std::make_unique<double[]>(nbfcontinua)},
        _gamma_contr{std::make_unique<double[]>(bfestimcount)} {
    groundcont_gamma_contr = {_groundcont_gamma_contr.get(), static_cast<std::size_t>(nbfcontinua_ground)};
    chi_bf_sum = {_chi_bf_sum.get(), static_cast<std::size_t>(nbfcontinua)};
    gamma_contr = {_gamma_contr.get(), static_cast<std::size_t>(bfestimcount)};
  }
#pragma clang unsafe_buffer_usage end

  constexpr Phixslist() = default;

 private:
  // unique ptrs are used instead of vectors for nvc++ compatibility,
  // since std::bad_alloc exceptions are not supported on device
  // (still true as of NVC++ 25.11)
  std::unique_ptr<double[]> _groundcont_gamma_contr;
  std::unique_ptr<double[]> _chi_bf_sum;
  std::unique_ptr<double[]> _gamma_contr;
  // NOLINTEND(*-avoid-c-arrays)
};

struct ContinuumOpacity {
  double nu{-1.};  // frequency at which opacity was calculated
  double chi_freefree_scatter{0.};
  double chi_freefree_heat{0.};
  double chi_boundfree{0.};
  int nonemptymgi{-1};
  int timestep{-1};
  Phixslist phixslist;

  constexpr ContinuumOpacity(const int nbfcontinua_ground, const int nbfcontinua, const int bfestimcount)
      : phixslist{nbfcontinua_ground, nbfcontinua, bfestimcount} {}

  constexpr ContinuumOpacity() = default;
  [[nodiscard]] constexpr auto total() const { return chi_freefree_scatter + chi_boundfree + chi_freefree_heat; }
};

DEVICE_FUNC void do_rpkt(Packet& pkt, double t2);
DEVICE_FUNC void emit_rpkt(Packet& pkt);
template <bool USECELLHISTANDUPDATEPHIXSLIST>
void calculate_chi_rpkt_cont(double nu_cmf, ContinuumOpacity& chi_rpkt_cont, int nonemptymgi);
extern template void calculate_chi_rpkt_cont<true>(double nu_cmf, ContinuumOpacity& chi_rpkt_cont, int nonemptymgi);
extern template void calculate_chi_rpkt_cont<false>(double nu_cmf, ContinuumOpacity& chi_rpkt_cont, int nonemptymgi);
[[nodiscard]] DEVICE_FUNC auto sample_planck_times_expansion_opacity(int nonemptymgi) -> double;
void allocate_expansionopacities();
void calculate_expansion_opacities(int nonemptymgi);
void MPI_Bcast_binned_opacities(ptrdiff_t nonemptymgi, int root_node_id);
auto calculate_chi_ffheat_nnionpart(int nonemptymgi) -> double;

[[nodiscard]] constexpr auto get_linedistance(const double prop_time, const double nu_cmf, const double nu_trans,
                                              const double dnu_on_dl) -> double {
  // distance from packet position to redshifting into line at frequency nu_trans

  if (nu_cmf <= nu_trans) {
    return 0;  // photon was propagated too far, make sure that we don't miss a line
  }

  if constexpr (USE_RELATIVISTIC_DOPPLER_SHIFT) {
    // With special relativity, the Doppler shift formula has an extra factor of 1/gamma in it,
    // which changes the distance reach a line resonance and creates a dependence
    // on packet position and direction

    // use linear interpolation of frequency along the path
    return (nu_trans - nu_cmf) / dnu_on_dl;
  }

  return CLIGHT * prop_time * ((nu_cmf / nu_trans) - 1);
}

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
  // entries in the line list and get the index of the closest line
  // to lower frequencies

  // will find the highest frequency (lowest index) line with nu_line <= nu_cmf
  // lower_bound matches the first element where the comparison function is false
  const int matchindex =
      static_cast<int>(std::ranges::lower_bound(linelistnu, nu_cmf, std::ranges::greater{}) - linelistnu.begin());

  return matchindex;
}

[[gnu::pure]] [[nodiscard]] inline auto keep_this_cont(int element, const int ion, const int level,
                                                       const int nonemptymgi, const float nnetot) -> bool {
  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    return grid::get_elem_abundance(nonemptymgi, element) > 0;
  }
  return ((get_nnion(nonemptymgi, element, ion) / nnetot > 1.e-6) || (level == 0));
}

#endif  // RPKT_H
