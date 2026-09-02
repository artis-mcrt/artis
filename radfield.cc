// The radiation field model: accumulates Monte Carlo estimators of the radiation field
// (full-spectrum and binned J_nu, plus detailed bound-bound Jb_lu estimators) and fits diluted
// blackbodies (W, T_R) per cell and per frequency bin for use in the photoionisation and heating rates.
//
// The multibin radiation field model and the photoionisation, bound-bound, and heating estimators
// are described by Shingles et al. (2020), MNRAS, 492, 2029, section 2.2, doi:10.1093/mnras/stz3412.

#include "radfield.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <fstream>
#include <iterator>
#include <print>
#include <span>
#include <utility>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>

#include <cstdint>

#include "toms748.h"
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "mpi_logging.h"
#include "rpkt.h"
#include "sn3d.h"

namespace radfield {

namespace {

constexpr double bins_T_R_min = 500;
constexpr double bins_T_R_max = 250000;

static_assert(RADFIELDBINS_T_E_SUPERBIN_NU_MAX >= RADFIELDBINS_NU_MAX,
              "The T_e superbin upper boundary must be greater than or equal to the upper boundary of the other bins");

static_assert(!DETAILED_BF_ESTIMATORS_ON || !USE_LUT_PHOTOION,
              "USE_LUT_PHOTOION must be false when DETAILED_BF_ESTIMATORS_ON is true");

std::vector<double> J_normfactor;

struct RadFieldBinSolution {
  // these two parameters are used in the current timestep, but were calculated
  // from the values of J and nuJ in the previous timestep
  float W;  // dilution (scaling) factor
  float T_R;  // radiation temperature
};

struct RadFieldBins {
  std::vector<double> J_raw;
  std::vector<double> nuJ_raw;

  void resize(const ptrdiff_t nonempty_npts_model) {
    reserve_resize(J_raw, nonempty_npts_model * RADFIELDBINCOUNT);
    reserve_resize(nuJ_raw, nonempty_npts_model * RADFIELDBINCOUNT);
  }
};

constexpr double radfieldbins_delta_nu =
    (RADFIELDBINS_NU_MAX - RADFIELDBINS_NU_MIN) / (RADFIELDBINCOUNT - 1);  // - 1 for the top super bin

RadFieldBins radfieldbins;

MPI_shared_array<float> radfieldbin_solutions_W;

MPI_shared_array<float> radfieldbin_solutions_T_R;

struct Jb_lu_estimator {
  double value = 0.;
  int contribcount = 0;
};

int detailed_linecount = 0;

// array of indices into the linelist[] array for selected lines
std::vector<int> detailed_lineindices;

std::vector<std::vector<Jb_lu_estimator>> prev_Jb_lu_normed{};  // value from the previous timestep
std::vector<std::vector<Jb_lu_estimator>> Jb_lu_raw{};  // unnormalised estimator for the current timestep

MPI_shared_array<float> prev_bfrate_normed{};  // values from the previous timestep
std::vector<double> bfrate_raw;  // unnormalised estimators for the current timestep

// for DETAILED_BF_ESTIMATORS_ON: maps an index into the allphixstargets arrays (the
// element/ion/level/phixstargetindex ordering) to the same continuum's index in the nu_edge-sorted
// globals::allcont arrays. Built in init() as the inverse of the sort permutation from setup_phixs_list().
MPI_shared_array<int> allcontindex_of_allphixstargetindex{};

// J and nuJ are accumulated and then normalised in-place
// i.e. be sure the normalisation has been applied (exactly once) before using the values here!

std::vector<double> J;  // after normalisation: [ergs/s/sr/cm2/Hz]
#ifdef DO_TITER
std::vector<double> J_reduced_save;
#endif

std::vector<double> nuJ;  // after normalisation: [ergs/s/sr/cm2]
#ifdef DO_TITER
std::vector<double> nuJ_reduced_save;
#endif

std::fstream radfieldfile;

constexpr auto get_bin_nu_upper(const int binindex) -> double {
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  if (binindex == RADFIELDBINCOUNT - 1) {
    return RADFIELDBINS_T_E_SUPERBIN_NU_MAX;
  }
  return RADFIELDBINS_NU_MIN + ((binindex + 1) * radfieldbins_delta_nu);
}

constexpr auto get_bin_nu_lower(const int binindex) -> double {
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);

  if (binindex > 0) {
    return get_bin_nu_upper(binindex - 1);
  }
  return RADFIELDBINS_NU_MIN;
}

// find the left-closed bin [nu_lower, nu_upper) that nu belongs to
constexpr auto select_bin(const double nu) -> int {
  if (nu < RADFIELDBINS_NU_MIN) {
    return -2;  // out of range, nu lower than lowest bin's lower boundary
  }
  if (nu >= RADFIELDBINS_T_E_SUPERBIN_NU_MAX) {
    // out of range, nu higher than highest bin's upper boundary
    return -1;
  }
  if (nu >= RADFIELDBINS_NU_MAX) {
    // in the superbin. separate case because the delta_nu is different to the other bins
    return RADFIELDBINCOUNT - 1;
  }

  const auto binindex = static_cast<int>(get_linearbinindex(nu, RADFIELDBINS_NU_MIN, radfieldbins_delta_nu));

  if (nu == get_bin_nu_upper(binindex)) {
    // exactly on the upper boundary of the bin, so add 1 to ensure we get the left-closed bin
    return binindex + 1;
  }
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < (RADFIELDBINCOUNT - 1));  // -1 because the superbin is a special case

  return binindex;
}

static_assert(get_bin_nu_lower(1) == get_bin_nu_upper(0));  // bins are contiguous
static_assert(select_bin(RADFIELDBINS_NU_MIN / 2) == -2);  // below the lowest bin
static_assert(select_bin(RADFIELDBINS_NU_MIN) == 0);  // bins are left-closed
static_assert(select_bin(RADFIELDBINS_NU_MAX) == RADFIELDBINCOUNT - 1);  // start of the T_e superbin
static_assert(select_bin(RADFIELDBINS_T_E_SUPERBIN_NU_MAX) == -1);  // at and above the superbin upper bound

// associate a Jb_lu estimator with a particular lineindex to be used instead of the general radiation field model
void add_detailed_line(const int lineindex) {
  detailed_linecount++;
  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    prev_Jb_lu_normed[nonemptymgi].push_back({.value = 0, .contribcount = 0});
    assert_always(detailed_linecount == std::ssize(prev_Jb_lu_normed[nonemptymgi]));

    // zero_estimators should do the next part anyway, but just to be sure:
    Jb_lu_raw[nonemptymgi].push_back({.value = 0, .contribcount = 0});
    assert_always(detailed_linecount == std::ssize(Jb_lu_raw[nonemptymgi]));
  }
  detailed_lineindices.push_back(lineindex);
  assert_always(detailed_linecount == std::ssize(detailed_lineindices));
}

// get the normalised J value for a bin
auto get_bin_J(const std::ptrdiff_t nonemptymgi, const int binindex) -> double {
  assert_testmodeonly(J_normfactor[nonemptymgi] > 0.0);
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  return radfieldbins.J_raw[(nonemptymgi * RADFIELDBINCOUNT) + binindex] * J_normfactor[nonemptymgi];
}

// get the normalised nuJ value for a bin
auto get_bin_nuJ(const std::ptrdiff_t nonemptymgi, const int binindex) -> double {
  assert_testmodeonly(J_normfactor[nonemptymgi] > 0.0);
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  return radfieldbins.nuJ_raw[(nonemptymgi * RADFIELDBINCOUNT) + binindex] * J_normfactor[nonemptymgi];
}

// get <nuJ> / <J> for a bin
auto get_bin_nu_bar(const std::ptrdiff_t nonemptymgi, const int binindex) -> double {
  const double nuJ_bin = get_bin_nuJ(nonemptymgi, binindex);
  const double J_bin = get_bin_J(nonemptymgi, binindex);
  return nuJ_bin / J_bin;
}

auto get_bin_W(const std::ptrdiff_t nonemptymgi, const int binindex) -> float {
  return radfieldbin_solutions_W[(nonemptymgi * RADFIELDBINCOUNT) + binindex];
}

auto get_bin_T_R(const std::ptrdiff_t nonemptymgi, const int binindex) -> float {
  return radfieldbin_solutions_T_R[(nonemptymgi * RADFIELDBINCOUNT) + binindex];
}

void update_bfestimators(const ptrdiff_t nonemptymgi, const double distance_e_cmf, const double nu_cmf,
                         const Phixslist& phixslist) {
  assert_testmodeonly(DETAILED_BF_ESTIMATORS_ON);

  // No doppler factor: like J, the photoionisation rate coefficient is a comoving-frame volume estimator, and the
  // frame conversion is already fully accounted for by e_cmf and evaluating sigma_bf at nu_cmf (consistent with the
  // LUT gammaestimator in rpkt.cc)
  const double distance_e_cmf_over_nu = distance_e_cmf / nu_cmf;

  // The packet has propagated since the phixslist was built, so nu_cmf has drifted (downwards) from the value the
  // list was truncated at. Re-deriving the bounds from the current nu_cmf can therefore narrow the range further
  // than the stored phixslist.bfestimbegin/bfestimend, and the contributions are only made over that narrower range.
  const auto bfestimcount = std::ssize(globals::bfestim_nu_edge);

  assert_testmodeonly(phixslist.bfestimend <= bfestimcount);
  const auto bfestimend =
      std::ranges::distance(globals::bfestim_nu_edge.begin(),
                            std::ranges::upper_bound(globals::bfestim_nu_edge.first(phixslist.bfestimend), nu_cmf));
  assert_testmodeonly(bfestimend <= bfestimcount);
  assert_testmodeonly(phixslist.bfestimbegin >= 0);
  // bfestimend was recalculated from the current (lower) nu_cmf, so it can fall below the bfestimbegin that
  // was stored for the higher nu_cmf. Clamp so that the subspan count below can never go negative (which
  // would wrap to a huge size_t), which just leaves an empty range with no contributions to make.
  const auto bfestimbegin_stored = std::min<ptrdiff_t>(phixslist.bfestimbegin, bfestimend);
  const auto bfestimbegin = std::ranges::distance(
      globals::bfestim_nu_edge.begin(),
      std::ranges::lower_bound(globals::bfestim_nu_edge.subspan(bfestimbegin_stored, bfestimend - bfestimbegin_stored),
                               nu_cmf / last_phixs_nuovernuedge));

  for (auto bfestimindex = bfestimbegin; bfestimindex < bfestimend; bfestimindex++) {
    atomicadd(bfrate_raw[(nonemptymgi * bfestimcount) + bfestimindex],
              phixslist.gamma_contr[bfestimindex] * distance_e_cmf_over_nu);
  }
}

// Series expansion used for computing the integral of the Planck function between x and infinity,
// where x = H * nu / (KB * T_R). epsrel is the relative error for convergence of the sum, i.e. the termination
// condition is term < sum * epsrel.
auto partial_planck_integral_x_to_inf(const double x, const double epsrel) -> double {
  assert_testmodeonly(x >= 0);
  if (x > 700) {
    return 0.0;  // e^-x underflows double precision anyway
  }

  double sum = 0.0;
  for (int n = 1; n < 1000; ++n) {
    const double n2 = n * n;
    const double n3 = n2 * n;
    const double n4 = n3 * n;

    const double term = std::exp(-n * x) * ((pow3(x) / n) + (3.0 * pow2(x) / n2) + (6.0 * x / n3) + (6.0 / n4));
    sum += term;

    if (term < sum * epsrel) {
      break;  // Convergence check
    }
  }
  return sum;
}

// Similar to partial_planck_integral_x_to_inf except with an extra factor of nu in the integrand
auto partial_nu_planck_integral_x_to_inf(const double x, const double epsrel) -> double {
  assert_testmodeonly(x >= 0);
  if (x > 700) {
    return 0.0;  // e^-x underflows double precision
  }

  double sum = 0.0;
  for (int n = 1; n < 1000; ++n) {
    const double n2 = n * n;
    const double n3 = n2 * n;
    const double n4 = n3 * n;
    const double n5 = n4 * n;

    const double term = std::exp(-n * x) *
                        ((pow4(x) / n) + (4.0 * pow3(x) / n2) + (12.0 * pow2(x) / n3) + (24.0 * x / n4) + (24.0 / n5));
    sum += term;

    if (term < sum * epsrel) {
      break;
    }
  }
  return sum;
}

// Calculate the intensity-weighted mean frequency without allowing the common Wien-tail exponential to underflow.
auto calculate_planck_mean_frequency(const double temperature, const double nu_low, const double nu_high) -> double {
  assert_testmodeonly(temperature > 0.);
  assert_testmodeonly(nu_low >= 0.);
  assert_testmodeonly(nu_high > nu_low);

  const double x_low = (H * nu_low) / (KB * temperature);
  const double x_high = (H * nu_high) / (KB * temperature);

  // For x >= 100, replacing 1 / (exp(x) - 1) with exp(-x) has a relative error below 4e-44 (and an absolute error
  // proportional to exp(-2x)). Factoring exp(-x_low) out of both moments prevents the individual Planck integrals
  // from underflowing. Using expm1 also keeps the difference between the two incomplete-gamma polynomials accurate
  // for narrow bins.
  constexpr double wien_tail_threshold = 100.;
  if (x_low >= wien_tail_threshold) {
    const auto moment3_tail_polynomial = [](const double x) { return pow3(x) + (3 * pow2(x)) + (6 * x) + 6; };
    const auto moment4_tail_polynomial = [](const double x) {
      return pow4(x) + (4 * pow3(x)) + (12 * pow2(x)) + (24 * x) + 24;
    };
    const auto scaled_bin_moment = [x_low, x_high](const auto tail_polynomial) {
      const double bin_width = x_high - x_low;
      const double exp_minus_bin_width = std::exp(-bin_width);
      const double polynomial_low = tail_polynomial(x_low);
      const double polynomial_high = tail_polynomial(x_high);
      return (-std::expm1(-bin_width) * polynomial_low) - (exp_minus_bin_width * (polynomial_high - polynomial_low));
    };

    const double planck_moment = scaled_bin_moment(moment3_tail_polynomial);
    const double nu_planck_moment = scaled_bin_moment(moment4_tail_polynomial);
    assert_testmodeonly(planck_moment > 0.);
    assert_testmodeonly(nu_planck_moment > 0.);

    return (KB * temperature / H) * (nu_planck_moment / planck_moment);
  }

  const double nu_planck_integral = calculate_planck_integral(temperature, nu_low, nu_high, true);
  const double planck_integral = calculate_planck_integral(temperature, nu_low, nu_high, false);
  return nu_planck_integral / planck_integral;
}

// nu_bar_planck_minus_estimator = nu_bar_planck(T_R) - nu_bar_estimator, where nu_bar is the intensity-weighted mean
// frequency in a bin, which is given by the ratio of nuJ and J estimators.
auto nu_bar_planck_minus_estimator(const double T_R, const int nonemptymgi, const int binindex) -> double {
  const double nu_lower = get_bin_nu_lower(binindex);
  const double nu_upper = get_bin_nu_upper(binindex);
  const double nu_bar_estimator = get_bin_nu_bar(nonemptymgi, binindex);

  const double nu_bar_planck_T_R = calculate_planck_mean_frequency(T_R, nu_lower, nu_upper);

  const double delta_nu_bar = nu_bar_planck_T_R - nu_bar_estimator;

  if (!std::isfinite(delta_nu_bar)) {
    printlnlog(
        "[warning] nu_bar_planck_minus_estimator: cell {} bin {}: delta_nu_bar is {:g}. T_R {:g} [K] nu_lower {:g} "
        "[Hz] nu_upper {:g} [Hz] nu_bar_planck_T_R {:g} [Hz] nu_bar_estimator {:g} [Hz]",
        grid::get_mgi_of_nonemptymgi(nonemptymgi), binindex, delta_nu_bar, T_R, nu_lower, nu_upper, nu_bar_planck_T_R,
        nu_bar_estimator);
  }

  return delta_nu_bar;
}

// Find the radiation temperature of a bin by matching the Planck mean frequency to the bin's nuJ/J estimator.
// nu_bar_planck_minus_estimator() increases with T_R, so a root outside [bins_T_R_min, bins_T_R_max] is clamped to
// whichever bound it lies beyond.
auto find_bin_T_R(const int nonemptymgi, const int binindex) -> float {
  const auto f_deltanubar = [nonemptymgi, binindex](const double T_R) {
    return nu_bar_planck_minus_estimator(T_R, nonemptymgi, binindex);
  };

  const double f_Tmin = f_deltanubar(bins_T_R_min);
  const double f_Tmax = f_deltanubar(bins_T_R_max);

  const bool invalid_values = (!std::isfinite(f_Tmin) || !std::isfinite(f_Tmax));

  // a sign change over the interval guarantees a root that the bracketing solver can find
  if (!invalid_values && f_Tmin * f_Tmax < 0) {
    constexpr double epsrel = 1e-4;
    const auto maxit = 100U;

    uintmax_t iteration_num = maxit;
    const auto result =
        toms748_solve(f_deltanubar, bins_T_R_min, bins_T_R_max, f_Tmin, f_Tmax, ftol<epsrel>, iteration_num);
    const auto T_R_solution = static_cast<float>(0.5 * (result.first + result.second));
    if (iteration_num >= maxit) {
      printlnlog(
          "[warning] find_bin_T_R: cell {} bin {}: T_R did not converge within {} iterations. interval [{:g}, "
          "{:g}] [K]",
          grid::get_mgi_of_nonemptymgi(nonemptymgi), binindex, iteration_num, result.first, result.second);
    }
    return T_R_solution;
  }
  if (invalid_values || f_Tmax < 0) {
    // the Planck mean frequency is still below the estimator at the top of the range (or the residual is not
    // finite), so any root lies above bins_T_R_max
    return bins_T_R_max;
  }
  // the residual is positive at both ends, so the Planck mean frequency already exceeds the estimator at the
  // bottom of the range and any root lies below bins_T_R_min
  return bins_T_R_min;
}

void set_params_fullspec(const int nonemptymgi, const int timestep) {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double nubar = nuJ[nonemptymgi] / J[nonemptymgi];
  if (!std::isfinite(nubar) || nubar == 0.) {
    printlnlog("[warning] T_R estimator infinite in cell {}, keep T_R, T_J, W of last timestep. J = {:g}. nuJ = {:g}",
               modelgridindex, J[nonemptymgi], nuJ[nonemptymgi]);
  } else {
    auto T_J = static_cast<float>(pow(J[nonemptymgi] * PI / STEBO, 1 / 4.));
    if (T_J > MAXTEMP) {
      printlnlog("[warning] temperature estimator T_J = {:g} exceeds T_max {:g} in cell {}. Setting T_J = T_max!", T_J,
                 MAXTEMP, modelgridindex);
      T_J = MAXTEMP;
    } else if (T_J < MINTEMP) {
      printlnlog("[warning] temperature estimator T_J = {:g} below T_min {:g} in cell {}. Setting T_J = T_min!", T_J,
                 MINTEMP, modelgridindex);
      T_J = MINTEMP;
    }
    grid::TJ_allcells[nonemptymgi] = T_J;

    auto T_R = static_cast<float>(H * nubar / KB / 3.832229494);
    if (T_R > MAXTEMP) {
      printlnlog("[warning] temperature estimator T_R = {:g} exceeds T_max {:g} in cell {}. Setting T_R = T_max!", T_R,
                 MAXTEMP, modelgridindex);
      T_R = MAXTEMP;
    } else if (T_R < MINTEMP) {
      printlnlog("[warning] temperature estimator T_R = {:g} below T_min {:g} in cell {}. Setting T_R = T_min!", T_R,
                 MINTEMP, modelgridindex);
      T_R = MINTEMP;
    }
    grid::TR_allcells[nonemptymgi] = T_R;

    const auto W = static_cast<float>(J[nonemptymgi] * PI / STEBO / pow4(T_R));
    grid::W_allcells[nonemptymgi] = W;

    printlnlog(
        "Full-spectrum fit radfield for cell {} at timestep {}: J {:g}, lambda_bar {:5.1f} [Angstrom], T_J {:g} [K], "
        "T_R {:g} [K], W {:g}",
        modelgridindex, timestep, J[nonemptymgi], 1e8 * CLIGHT / nubar, T_J, T_R, W);
  }
}

auto get_allcontindex(const int element, const int lowerion, const int lower, const int phixstargetindex) -> int {
  // direct lookup via the precomputed inverse of the nu_edge sort permutation (built in init())
  if (allcontindex_of_allphixstargetindex.empty()) {
    // there are no continua (or init() has not been called), so nothing can be found
    return -1;
  }
  const auto allphixstargetindex =
      get_allphixstargetindex(get_uniquelevelindex(element, lowerion, lower), phixstargetindex);
  assert_testmodeonly(allphixstargetindex >= 0);
  assert_testmodeonly(allphixstargetindex < allcontindex_of_allphixstargetindex.ssize());
  return allcontindex_of_allphixstargetindex[allphixstargetindex];
}

void write_to_file(const int nonemptymgi, const int timestep) {
  assert_always(MULTIBIN_RADFIELD_MODEL_ON);
  const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
#ifdef _OPENMP
#pragma omp critical(out_file)
  {
#endif

    for (int binindex = -1 - detailed_linecount; binindex < RADFIELDBINCOUNT; binindex++) {
      double nu_lower = 0.;
      double nu_upper = 0.;
      double nuJ_out = 0.;
      double J_out = 0.;
      float T_R = 0.;
      float W = 0.;
      double J_nu_bar = 0.;

      if (binindex >= 0) {
        nu_lower = get_bin_nu_lower(binindex);
        nu_upper = get_bin_nu_upper(binindex);
        nuJ_out = get_bin_nuJ(nonemptymgi, binindex);
        J_out = get_bin_J(nonemptymgi, binindex);
        T_R = get_bin_T_R(nonemptymgi, binindex);
        W = get_bin_W(nonemptymgi, binindex);
        J_nu_bar = J_out / (nu_upper - nu_lower);
      } else if (binindex == -1) {  // bin -1 is the full spectrum fit
        nuJ_out = nuJ[nonemptymgi];
        J_out = J[nonemptymgi];
        T_R = grid::TR_allcells[nonemptymgi];
        W = grid::W_allcells[nonemptymgi];
      } else {  // use binindex < -1 for detailed line Jb_lu estimators
        const int jblueindex = -2 - binindex;  // -2 is the first detailed line, -3 is the second, etc
        const int lineindex = detailed_lineindices[jblueindex];
        const double nu_trans = globals::linelist.nu[lineindex];
        nu_lower = nu_trans;
        nu_upper = nu_trans;
        nuJ_out = -1.;
        J_out = -1.;
        T_R = -1.;
        W = -1.;
        J_nu_bar = prev_Jb_lu_normed[nonemptymgi][jblueindex].value;
      }

      std::println(radfieldfile, "{:d} {:d} {:d} {:.5e} {:.5e} {:.3e} {:.3e} {:.3e} {:.1f} {:.5e}", timestep,
                   modelgridindex, binindex, nu_lower, nu_upper, nuJ_out, J_out, J_nu_bar, T_R, W);
    }
    radfieldfile.flush();
#ifdef _OPENMP
  }
#endif
}

}  // anonymous namespace

// Compute the integral of the Planck function (or nu times the Planck function) between frequency nu=nu_low to
// nu=nu_high. Units are ergs/s/sr/cm2 for the integral of the Planck function, and ergs/s2/sr/cm2 for the integral of
// nu times the Planck function.
auto calculate_planck_integral(const double temperature, const double nu_low, const double nu_high, const bool times_nu)
    -> double {
  if (temperature <= 0) {
    return 0.0;
  }

  constexpr double epsrel = 1e-15;  // relative error for convergence of the series expansion

  // integration variable x = H * nu / (KB * T_R)
  const double x_low = (H * nu_low) / (KB * temperature);
  const double x_high = (H * nu_high) / (KB * temperature);

  if (times_nu) {
    const double constant_factor = (2.0 * pow5(KB) * pow5(temperature)) / (pow4(H) * pow2(CLIGHT));
    const auto low_to_inf = partial_nu_planck_integral_x_to_inf(x_low, epsrel);
    const auto high_to_inf = partial_nu_planck_integral_x_to_inf(x_high, epsrel);

    return constant_factor * (low_to_inf - high_to_inf);
  }

  const double constant_factor = (2.0 * pow4(KB) * pow4(temperature)) / (pow3(H) * pow2(CLIGHT));
  const auto low_to_inf = partial_planck_integral_x_to_inf(x_low, epsrel);
  const auto high_to_inf = partial_planck_integral_x_to_inf(x_high, epsrel);

  return constant_factor * (low_to_inf - high_to_inf);
}

void init() {
  // this should be called only after the atomic data is in memory

  const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();

  reserve_resize(J_normfactor, nonempty_npts_model + 1);
  reserve_resize(J, nonempty_npts_model + 1);

  // J and nuJ are accumulated and then normalised in-place
  // i.e. be sure the normalisation has been applied (exactly once) before using the values here!
  reserve_resize(nuJ, nonempty_npts_model + 1);

#ifdef DO_TITER
  reserve_resize(J_reduced_save, nonempty_npts_model + 1);
  reserve_resize(nuJ_reduced_save, nonempty_npts_model + 1);
#endif

  reserve_resize(prev_Jb_lu_normed, nonempty_npts_model);
  reserve_resize(Jb_lu_raw, nonempty_npts_model);

  detailed_linecount = 0;

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    for (int i = 0; i < globals::nlines; i++) {
      const int element = globals::linelist.elementindex[i];
      const int Z = get_atomicnumber(element);
      if (Z == 26) {
        const int ion = globals::linelist.ionindex[i];
        const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
        const auto lower_uniquelevelindex = globals::linelist.uniquelevelindex_lower[i];
        const int lowerlevel = lower_uniquelevelindex - ionuniquelevelindexstart;

        bool addline = false;
        if (lowerlevel <= 15) {  // selection for detailed estimators (by lower-level index)
          addline = true;
        }

        if (addline) {
          add_detailed_line(i);
        }
      }
    }
    printlnlog("There are {} lines with detailed Jblue_lu estimators.", detailed_linecount);
  }

  printlog("DETAILED_BF_ESTIMATORS {}", DETAILED_BF_ESTIMATORS_ON ? "ON" : "OFF");
  if (DETAILED_BF_ESTIMATORS_ON) {
    printlnlog(" from timestep {}", DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP);
  } else {
    printlnlog("");
  }

  if (MULTIBIN_RADFIELD_MODEL_ON) {
    printlnlog("The multibin radiation field is being used from timestep {} onwards.", FIRST_NLTE_RADFIELD_TIMESTEP);

    printlnlog(
        "Initialising multibin radiation field with {} bins from ({:.2f} [eV], {:6.1f} [A]) to ({:.2f} [eV], "
        "{:6.1f} [A]) and a T_e superbin up to ({:.2f} [eV], {:6.1f} [A]).",
        RADFIELDBINCOUNT - 1, H * RADFIELDBINS_NU_MIN / EV, 1e8 * CLIGHT / RADFIELDBINS_NU_MIN,
        H * RADFIELDBINS_NU_MAX / EV, 1e8 * CLIGHT / RADFIELDBINS_NU_MAX, H * RADFIELDBINS_T_E_SUPERBIN_NU_MAX / EV,
        1e8 * CLIGHT / RADFIELDBINS_T_E_SUPERBIN_NU_MAX);
    if (grid::get_ndo_nonempty(globals::my_rank) > 0) {
      assert_always(!radfieldfile.is_open());
      radfieldfile = open_rank_outfile("radfield");
      std::println(radfieldfile, "timestep modelgridindex bin_num nu_lower nu_upper nuJ J J_nu_avg T_R W");
      radfieldfile.flush();
    }

    const size_t mem_usage_bins = nonempty_npts_model * RADFIELDBINCOUNT * ((2 * sizeof(double)) + sizeof(int));
    radfieldbins.resize(nonempty_npts_model);

    printlnlog("[info] mem_usage: radiation field bin accumulators for non-empty cells occupy {:.3f} MB",
               mem_usage_bins / 1024. / 1024.);

    radfieldbin_solutions_W = MPI_shared_array<float>(nonempty_npts_model * RADFIELDBINCOUNT);

    radfieldbin_solutions_T_R = MPI_shared_array<float>(nonempty_npts_model * RADFIELDBINCOUNT);

    const size_t mem_usage_bin_solutions = nonempty_npts_model * RADFIELDBINCOUNT * sizeof(RadFieldBinSolution);
    printlnlog(
        "[info] mem_usage: radiation field bin solutions for non-empty cells occupy {:.3f} MB (node shared memory)",
        mem_usage_bin_solutions / 1024. / 1024.);
  } else {
    printlnlog("The radiation field model is a full-spectrum fit to a single dilute blackbody TR & W.");
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    if (globals::nbfcontinua > 0) {
      // build the O(1) lookup used by get_allcontindex(): for each continuum in the nu_edge-sorted
      // allcont list, record its position under the allphixstargets (element/ion/level/target) ordering
      auto temp_allcontindex = MPI_shared_array<int>(globals::nbfcontinua, -1);
      if (globals::rank_in_node == 0) {
        for (int allcontindex = 0; allcontindex < globals::nbfcontinua; allcontindex++) {
          const auto allphixstargetindex = get_allphixstargetindex(globals::allcont.uniquelevelindex[allcontindex],
                                                                   globals::allcont.phixstargetindex[allcontindex]);
          temp_allcontindex[allphixstargetindex] = allcontindex;
        }
      }
      MPI_Barrier_node();
      allcontindex_of_allphixstargetindex = std::move(temp_allcontindex);
    }

    const auto bfestimcount = std::ssize(globals::bfestim_nu_edge);
    prev_bfrate_normed = MPI_shared_array<float>(nonempty_npts_model * bfestimcount);
    if (globals::rank_in_node == 0) {
      std::ranges::fill(prev_bfrate_normed, 0.);
    }
    MPI_Barrier_node();
    printlnlog("[info] mem_usage: detailed bf estimators for non-empty cells occupy {:.3f} MB (node shared memory)",
               nonempty_npts_model * bfestimcount * sizeof(float) / 1024. / 1024.);

    reserve_resize(bfrate_raw, nonempty_npts_model * bfestimcount);

    printlnlog("[info] mem_usage: detailed bf estimator acculumators for non-empty cells occupy {:.3f} MB",
               nonempty_npts_model * bfestimcount * sizeof(double) / 1024. / 1024.);
  }

  zero_estimators();

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    MPI_Barrier_node();
    if (globals::rank_in_node == 0) {
      std::ranges::fill(radfieldbin_solutions_W, -1.);
      std::ranges::fill(radfieldbin_solutions_T_R, -1.);
    }
    MPI_Barrier_node();
  }
}

// Initialise estimator arrays which hold the last time steps values (used to damp out
// fluctuations over timestep iterations if DO_TITER is defined) to -1.
void initialise_prev_titer_photoionestimators() {
#ifdef DO_TITER
  std::ranges::fill(globals::ffheatingestimator_save, -1.);
  std::ranges::fill(globals::colheatingestimator_save, -1.);
  std::ranges::fill(J_reduced_save, -1.);
  std::ranges::fill(nuJ_reduced_save, -1.);
  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
        const int groundcontindex = get_groundcontindex(element, ion);
        if (groundcontindex < 0) {
          // an ion without a ground photoionisation table has no estimator slot
          continue;
        }
        if constexpr (USE_LUT_PHOTOION) {
          globals::gammaestimator_save[(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) +
                                       groundcontindex] = -1.;
        }
        if constexpr (USE_ION_BFHEATING_ESTIMATORS) {
          globals::bfheatingestimator_save[(static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) +
                                           groundcontindex] = -1.;
        }
      }
    }
  }
#endif
}

auto get_Jblueindex(const int lineindex) -> int {
  // returns -1 if the line does not have a Jblue estimator
  if constexpr (!DETAILED_LINE_ESTIMATORS_ON) {
    return -1;
  }

  // use a binary search, assuming the list is sorted

  int low = 0;
  int high = detailed_linecount - 1;
  while (low <= high) {
    const int mid = low + ((high - low) / 2);
    if (detailed_lineindices[mid] < lineindex) {
      low = mid + 1;
    } else if (detailed_lineindices[mid] > lineindex) {
      high = mid - 1;
    } else {
      assert_always(mid < detailed_linecount);
      return mid;
    }
  }

  return -1;
}

auto get_Jb_lu(const int nonemptymgi, const int jblueindex) -> double {
  assert_always(jblueindex >= 0);
  assert_always(jblueindex < detailed_linecount);
  return prev_Jb_lu_normed[nonemptymgi][jblueindex].value;
}

// set up the new bins and clear the estimators in preparation for a timestep
void zero_estimators() {
  std::ranges::fill(J_normfactor, -1.0);
  std::ranges::fill(J, 0.0);
  std::ranges::fill(nuJ, 0.0);
  std::ranges::fill(bfrate_raw, 0.0);

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    std::ranges::fill(radfieldbins.J_raw, 0.0);
    std::ranges::fill(radfieldbins.nuJ_raw, 0.0);
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      std::fill_n(Jb_lu_raw[nonemptymgi].data(), detailed_linecount, Jb_lu_estimator{.value = 0., .contribcount = 0});
    }
  }
}

DEVICE_FUNC void update_estimators(const ptrdiff_t nonemptymgi, const double distance_e_cmf, const double nu_cmf,
                                   const Phixslist& phixslist, const bool thickcell) {
  if (distance_e_cmf == 0) {
    return;
  }

  atomicadd(J[nonemptymgi], distance_e_cmf);
  atomicadd(nuJ[nonemptymgi], distance_e_cmf * nu_cmf);

  if (thickcell) {
    return;
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    update_bfestimators(nonemptymgi, distance_e_cmf, nu_cmf, phixslist);
  }

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    const int binindex = select_bin(nu_cmf);

    if (binindex >= 0) {
      const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
      atomicadd(radfieldbins.J_raw[mgibinindex], distance_e_cmf);
      atomicadd(radfieldbins.nuJ_raw[mgibinindex], distance_e_cmf * nu_cmf);
    }
  }
}

DEVICE_FUNC void update_lineestimator(const int nonemptymgi, const int lineindex, const double increment) {
  if constexpr (!DETAILED_LINE_ESTIMATORS_ON) {
    return;
  }

  const int jblueindex = get_Jblueindex(lineindex);
  if (jblueindex >= 0) {
    atomicadd(Jb_lu_raw[nonemptymgi][jblueindex].value, increment);
    atomicadd(Jb_lu_raw[nonemptymgi][jblueindex].contribcount, 1);
  }
}

// mean intensity J_nu [ergs/s/sr/cm2/Hz]
DEVICE_FUNC auto radfield(const double nu, const int nonemptymgi) -> double {
  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    if (globals::timestep >= FIRST_NLTE_RADFIELD_TIMESTEP) {
      const int binindex = select_bin(nu);
      if (binindex >= 0) {
        const auto W = get_bin_W(nonemptymgi, binindex);
        if (W >= 0.) {
          return W * planck(nu, get_bin_T_R(nonemptymgi, binindex));
        }
      }
      return 0.;
    }
  }
  // full spectrum fit to a single dilute blackbody
  return grid::W_allcells[nonemptymgi] * planck(nu, grid::TR_allcells[nonemptymgi]);
}

// Fit the radiation field parameters of one cell to the estimators accumulated over the last timestep:
// the full-spectrum diluted blackbody (W, T_R), and with MULTIBIN_RADFIELD_MODEL_ON a separate (W, T_R)
// per frequency bin.
void fit_parameters(const int nonemptymgi, const int timestep) {
  set_params_fullspec(nonemptymgi, timestep);
  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    assert_always(J_normfactor[nonemptymgi] >= 0.);

    double J_bin_sum = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      J_bin_sum += get_bin_J(nonemptymgi, binindex);
    }
    const auto mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);

    printlnlog(
        "timestep {} cell {}: radfield bins sum to J of {:g} ({:.1f}% of total J). Finding parameters for {} bins...",
        timestep, mgi, J_bin_sum, 100. * J_bin_sum / J[nonemptymgi], RADFIELDBINCOUNT);

    double J_bin_max = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const double J_bin = get_bin_J(nonemptymgi, binindex);
      J_bin_max = std::max(J_bin_max, J_bin);
    }

    int count_T_R_min = 0;
    int count_T_R_max = 0;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const double nu_lower = get_bin_nu_lower(binindex);
      const double nu_upper = get_bin_nu_upper(binindex);
      const double J_bin = get_bin_J(nonemptymgi, binindex);
      float T_R_bin = -1.;
      float W_bin = -1.;

      if (J_bin > 0) {
        if (binindex == RADFIELDBINCOUNT - 1) {
          T_R_bin = grid::Te_allcells[nonemptymgi];
        } else {
          T_R_bin = find_bin_T_R(nonemptymgi, binindex);
          if (T_R_bin <= bins_T_R_min) {
            count_T_R_min++;
          } else if (T_R_bin >= bins_T_R_max) {
            count_T_R_max++;
          }
        }
        double planck_integral_result = calculate_planck_integral(T_R_bin, nu_lower, nu_upper, false);

        W_bin = static_cast<float>(J_bin / planck_integral_result);

        if (W_bin > 1e4 || !std::isfinite(W_bin)) {
          const auto W_bin_first = W_bin;
          planck_integral_result = calculate_planck_integral(bins_T_R_max, nu_lower, nu_upper, false);
          W_bin = static_cast<float>(J_bin / planck_integral_result);
          if (W_bin > 1e4) {
            printlnlog(
                "[warning] cell {} ts {} bin {}: W {:g} too high or non-finite (T_R {:7.1f} [K], J_bin {:g}) and W "
                "{:g} still too high after retry at T_R_max {:g} [K]. Zeroing bin (T_R set to sentinel -99)",
                mgi, timestep, binindex, W_bin_first, T_R_bin, J_bin, W_bin, bins_T_R_max);
            T_R_bin = -99.;
            W_bin = 0.;
          } else {
            printlnlog(
                "[warning] cell {} ts {} bin {}: W {:g} too high or non-finite (T_R {:7.1f} [K], J_bin {:g}); "
                "continuing with W {:g} at T_R_max {:g} [K]",
                mgi, timestep, binindex, W_bin_first, T_R_bin, J_bin, W_bin, bins_T_R_max);
            T_R_bin = bins_T_R_max;
          }
        }
      } else {
        T_R_bin = 0.;
        W_bin = 0.;
      }

      const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
      radfieldbin_solutions_T_R[mgibinindex] = T_R_bin;
      radfieldbin_solutions_W[mgibinindex] = W_bin;
    }

    if (count_T_R_min > 0 || count_T_R_max > 0) {
      printlnlog(
          "[warning] timestep {} cell {}: Some bin T_R values were clamped. {} bins at T_R=T_R_min={} and {} bins at "
          "T_R=T_R_max={}",
          timestep, mgi, count_T_R_min, bins_T_R_min, count_T_R_max, bins_T_R_max);
    }

    write_to_file(nonemptymgi, timestep);
  }
}

void set_J_normfactor(const int nonemptymgi, const double normfactor) { J_normfactor[nonemptymgi] = normfactor; }

void normalise_J(const int nonemptymgi, const double estimator_normfactor_over4pi) {
  assert_always(std::isfinite(J[nonemptymgi]));
  J[nonemptymgi] *= estimator_normfactor_over4pi;
  for (int i = 0; i < detailed_linecount; i++) {
    prev_Jb_lu_normed[nonemptymgi][i].value = Jb_lu_raw[nonemptymgi][i].value * estimator_normfactor_over4pi;
    prev_Jb_lu_normed[nonemptymgi][i].contribcount = Jb_lu_raw[nonemptymgi][i].contribcount;
  }
}

void normalise_bf_estimators(const int nts, const int nts_prev, const int titer, const double deltat) {
  // these conditions are the same on every rank, so all node ranks reach the barrier below together
  if (globals::lte_iteration) {
    return;
  }
  if (nts == globals::timestep_initial && titer == 0) {
    return;
  }
  if (globals::rank_in_node == 0) {
    const auto bfestimcount = std::ssize(globals::bfestim_nu_edge);
    const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();
    for (auto nonemptymgi = 0Z; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
      if (grid::thick_allcells[nonemptymgi] == grid::CellThickness::THICK) {
        continue;
      }
      const auto mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
      const double deltaV =
          grid::get_modelcell_assocvolume_tmin(mgi) * pow3(globals::timesteps[nts_prev].mid / globals::tmin);
      const double estimator_normfactor = 1 / deltaV / deltat / globals::nprocs;
      for (int i = 0; i < bfestimcount; i++) {
        const auto mgibfindex = (nonemptymgi * bfestimcount) + i;
        prev_bfrate_normed[mgibfindex] = static_cast<float>(bfrate_raw[mgibfindex] * (estimator_normfactor / H));
      }
    }
  }
  // every rank on the node reads prev_bfrate_normed during the upcoming grid update, so
  // wait here for the node leader to finish writing the node-shared array
  MPI_Barrier_node();
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_bfrate_estimator(const int element, const int lowerion,
                                                                  const int lower, const int phixstargetindex,
                                                                  const int nonemptymgi) -> double {
  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    const int allcontindex = get_allcontindex(element, lowerion, lower, phixstargetindex);
    if (allcontindex >= 0) {
      const auto bfestimindex = globals::allcont.bfestimindex[allcontindex];
      if (bfestimindex >= 0) {
        return prev_bfrate_normed[(nonemptymgi * std::ssize(globals::bfestim_nu_edge)) + bfestimindex];
      }
    }
  }

  return -1.;
}

void normalise_nuJ(const int nonemptymgi, const double estimator_normfactor_over4pi) {
  assert_always(std::isfinite(nuJ[nonemptymgi]));
  nuJ[nonemptymgi] *= estimator_normfactor_over4pi;
}

// Get the radiation temperature of a cell from the mean intensity alone, by equating J to the
// Stefan-Boltzmann law (T_J = (pi J / sigma)^1/4). Non-finite values keep the previous timestep's value,
// and the result is clamped to [MINTEMP, MAXTEMP].
auto get_T_J_from_J(const int nonemptymgi) -> float {
  const auto T_J = static_cast<float>(pow(J[nonemptymgi] * PI / STEBO, 1. / 4.));
  if (!std::isfinite(T_J)) {
    // keep old value of T_J
    const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    printlnlog("[warning] get_T_J_from_J: T_J estimator infinite in cell {}, use value of last timestep",
               modelgridindex);
    return grid::TJ_allcells[nonemptymgi];
  }
  // Make sure that T is in the allowed temperature range.
  if (T_J > MAXTEMP) {
    printlnlog(
        "[warning] get_T_J_from_J: cell {} ts {}: T_J would be {:.1f} [K] > MAXTEMP. Clamping to MAXTEMP = {:.0f} [K]",
        grid::get_mgi_of_nonemptymgi(nonemptymgi), globals::timestep, T_J, MAXTEMP);
    return MAXTEMP;
  }
  if (T_J < MINTEMP) {
    printlnlog(
        "[warning] get_T_J_from_J: cell {} ts {}: T_J would be {:.1f} [K] < MINTEMP. Clamping to MINTEMP = {:.0f} [K]",
        grid::get_mgi_of_nonemptymgi(nonemptymgi), globals::timestep, T_J, MINTEMP);
    return MINTEMP;
  }
  return T_J;
}

#ifdef DO_TITER
void titer_J(const int nonemptymgi) { titer_average(J[nonemptymgi], J_reduced_save[nonemptymgi]); }

void titer_nuJ(const int nonemptymgi) { titer_average(nuJ[nonemptymgi], nuJ_reduced_save[nonemptymgi]); }
#endif

// reduce and broadcast (allreduce) the estimators for J and nuJ in all bins
void reduce_estimators() {
  MPI_Allreduce_safe(J, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce_safe(nuJ, MPI_SUM, MPI_COMM_WORLD);

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    // reduce all ranks on each node first, then reduce the node leaders
    // then broadcast the final result to all ranks on each node
    // this seems necessary to avoid congestion compared to a single MPI_Allreduce when using many ranks per node
    MPI_Reduce_safe(bfrate_raw, MPI_SUM, 0, globals::mpi_comm_node);
    if (globals::rank_in_node == 0) {
      MPI_Allreduce_safe(bfrate_raw, MPI_SUM, globals::mpi_comm_internode);
    }
    MPI_Bcast_safe(bfrate_raw, 0, globals::mpi_comm_node);
  }

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    const auto sys_time_start_reduction = std::chrono::steady_clock::now();
    printlog("Reducing binned radiation field estimators");

    MPI_Allreduce_safe(radfieldbins.J_raw, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce_safe(radfieldbins.nuJ_raw, MPI_SUM, MPI_COMM_WORLD);

    const auto duration_reduction =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_reduction).count();
    printlnlog(" (took {:.1f} s)", duration_reduction);
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    const auto sys_time_start_reduction = std::chrono::steady_clock::now();
    printlog("Reducing detailed line estimators");

    for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
        MPI_Allreduce_safe(Jb_lu_raw[nonemptymgi][jblueindex].value, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce_safe(Jb_lu_raw[nonemptymgi][jblueindex].contribcount, MPI_SUM, MPI_COMM_WORLD);
      }
    }
    const auto duration_reduction =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_reduction).count();
    printlnlog(" (took {:.1f} s)", duration_reduction);
  }
  MPI_Barrier_allranks();
}

// broadcast the radiation field parameters of the cells that belong to the root rank to all ranks. The caller
// puts a barrier after this call, so that the other ranks on a node read the shared arrays after the broadcast.
void do_MPI_Bcast(const ptrdiff_t nstart_nonempty, const ptrdiff_t ndo_nonempty, const int root,
                  const int root_node_id) {
  MPI_Bcast_safe(std::span{J_normfactor}.subspan(nstart_nonempty, ndo_nonempty), root, MPI_COMM_WORLD);

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    if (globals::rank_in_node == 0) {
      MPI_Bcast_safe(
          radfieldbin_solutions_W.subspan(nstart_nonempty * RADFIELDBINCOUNT, ndo_nonempty * RADFIELDBINCOUNT),
          root_node_id, globals::mpi_comm_internode);
      MPI_Bcast_safe(
          radfieldbin_solutions_T_R.subspan(nstart_nonempty * RADFIELDBINCOUNT, ndo_nonempty * RADFIELDBINCOUNT),
          root_node_id, globals::mpi_comm_internode);
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    for (auto nonemptymgi = nstart_nonempty; nonemptymgi < (nstart_nonempty + ndo_nonempty); nonemptymgi++) {
      for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
        MPI_Bcast_safe(prev_Jb_lu_normed[nonemptymgi][jblueindex].value, root, MPI_COMM_WORLD);
        MPI_Bcast_safe(prev_Jb_lu_normed[nonemptymgi][jblueindex].contribcount, root, MPI_COMM_WORLD);
      }
    }
  }
}

void write_restart_data(FILE* gridsave_file) {
  printlog("binned radiation field and detailed lines, ");

  fprintf(gridsave_file, "%d\n", 30490824);  // special number marking the beginning of radfield data

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    fprintf(gridsave_file, "%d %la %la %la %la\n", RADFIELDBINCOUNT, RADFIELDBINS_NU_MIN, RADFIELDBINS_NU_MAX,
            bins_T_R_min, bins_T_R_max);

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      fprintf(gridsave_file, "%d %la\n", binindex, get_bin_nu_upper(binindex));
    }
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    const int nbfcontinua = globals::nbfcontinua;
    fprintf(gridsave_file, "%d\n", nbfcontinua);

    const int bfestimcount = static_cast<int>(globals::bfestim_nu_edge.size());
    fprintf(gridsave_file, "%d\n", bfestimcount);

    for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      fprintf(gridsave_file, "%d\n", nonemptymgi);
      for (int i = 0; i < bfestimcount; i++) {
        fprintf(gridsave_file, "%a ", prev_bfrate_normed[(nonemptymgi * bfestimcount) + i]);
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    fprintf(gridsave_file, "%d\n", detailed_linecount);

    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      fprintf(gridsave_file, "%d ", detailed_lineindices[jblueindex]);
    }
  }

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    assert_testmodeonly(nonemptymgi >= 0);
    fprintf(gridsave_file, "%d %la\n", nonemptymgi, J_normfactor[nonemptymgi]);

    if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
        const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
        fprintf(gridsave_file, "%la %la %a %a\n", radfieldbins.J_raw[mgibinindex], radfieldbins.nuJ_raw[mgibinindex],
                radfieldbin_solutions_W[mgibinindex], radfieldbin_solutions_T_R[mgibinindex]);
      }
    }

    if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
      for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
        fprintf(gridsave_file, "%la %d\n", Jb_lu_raw[nonemptymgi][jblueindex].value,
                Jb_lu_raw[nonemptymgi][jblueindex].contribcount);
      }
    }
  }
  fprintf(gridsave_file, "%d\n", 42809403);  // special number marking the end of radfield data
}

void read_restart_data(FILE* gridsave_file) {
  printlnlog("Reading restart data for radiation field");

  int code_check = 0;
  assert_always(fscanf(gridsave_file, "%d\n", &code_check) == 1);
  assert_always(code_check == 30490824);

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    double T_R_min_in{NAN};
    double T_R_max_in{NAN};
    double nu_min_in{NAN};
    double nu_max_in{NAN};
    int bincount_in = 0;
    assert_always(fscanf(gridsave_file, "%d %la %la %la %la\n", &bincount_in, &nu_min_in, &nu_max_in, &T_R_min_in,
                         &T_R_max_in) == 5);

    double nu_lower_first_ratio = nu_min_in / RADFIELDBINS_NU_MIN;
    if (nu_lower_first_ratio > 1.0) {
      nu_lower_first_ratio = 1 / nu_lower_first_ratio;
    }

    double nu_upper_last_ratio = nu_max_in / RADFIELDBINS_NU_MAX;
    if (nu_upper_last_ratio > 1.0) {
      nu_upper_last_ratio = 1 / nu_upper_last_ratio;
    }

    if (bincount_in != RADFIELDBINCOUNT || T_R_min_in != bins_T_R_min || T_R_max_in != bins_T_R_max ||
        nu_lower_first_ratio < 0.999 || nu_upper_last_ratio < 0.999) {
      printlnlog("[error] gridsave file specifies {} bins, nu_min {} nu_max {} T_R_min {} T_R_max {}", bincount_in,
                 nu_min_in, nu_max_in, T_R_min_in, T_R_max_in);
      fatal_crash("require {} bins, RADFIELDBINS_NU_MIN {:g} RADFIELDBINS_NU_MAX {:g} T_R_min {:g} T_R_max {:g}",
                  RADFIELDBINCOUNT, RADFIELDBINS_NU_MIN, RADFIELDBINS_NU_MAX, bins_T_R_min, bins_T_R_max);
    }

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      int binindex_in = 0;
      double nu_upper_in = NAN;
      assert_always(fscanf(gridsave_file, "%d %la\n", &binindex_in, &nu_upper_in) == 2);
      assert_always(binindex_in == binindex);
      assert_always(nu_upper_in == get_bin_nu_upper(binindex));
    }
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    int gridsave_nbf_in = 0;
    assert_always(fscanf(gridsave_file, "%d\n", &gridsave_nbf_in) == 1);
    assert_always(gridsave_nbf_in == globals::nbfcontinua);

    const auto bfestimcount = std::ssize(globals::bfestim_nu_edge);
    int gridsave_nbfestim_in = 0;
    assert_always(fscanf(gridsave_file, "%d\n", &gridsave_nbfestim_in) == 1);
    assert_always(gridsave_nbfestim_in == bfestimcount);

    for (auto nonemptymgi = 0Z; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      int nonemptymgi_in = 0;
      assert_always(fscanf(gridsave_file, "%d\n", &nonemptymgi_in) == 1);
      assert_always(nonemptymgi_in == nonemptymgi);
      for (int i = 0; i < bfestimcount; i++) {
        float bfrate_normed = 0;
        assert_always(fscanf(gridsave_file, "%a ", &bfrate_normed) == 1);

        if (globals::rank_in_node == 0) {
          prev_bfrate_normed[(nonemptymgi * bfestimcount) + i] = bfrate_normed;
        }
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    int detailed_linecount_in = 0;
    assert_always(fscanf(gridsave_file, "%d\n", &detailed_linecount_in) == 1);

    if (detailed_linecount_in != detailed_linecount) {
      fatal_crash("gridsave file specifies {} detailed lines but this simulation has {}.", detailed_linecount_in,
                  detailed_linecount);
    }

    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      assert_always(fscanf(gridsave_file, "%d ", &detailed_lineindices[jblueindex]) == 1);
    }
  }

  for (auto nonemptymgi = 0Z; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    int nonemptymgi_in = 0;
    assert_always(fscanf(gridsave_file, "%d %la\n", &nonemptymgi_in, &J_normfactor[nonemptymgi]) == 2);
    assert_always(nonemptymgi_in == nonemptymgi);

    if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
        const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
        float W = 0;
        float T_R = 0;
        assert_always(fscanf(gridsave_file, "%la %la %a %a\n", &radfieldbins.J_raw[mgibinindex],
                             &radfieldbins.nuJ_raw[mgibinindex], &W, &T_R) == 4);
        if (globals::rank_in_node == 0) {
          radfieldbin_solutions_W[mgibinindex] = W;
          radfieldbin_solutions_T_R[mgibinindex] = T_R;
        }
      }
    }

    if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
      for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
        assert_always(fscanf(gridsave_file, "%la %d\n", &Jb_lu_raw[nonemptymgi][jblueindex].value,
                             &Jb_lu_raw[nonemptymgi][jblueindex].contribcount) == 2);
        // normalise_J() is skipped on the timestep that a run resumes from, so the normalised values
        // that the macro atom reads have to be rebuilt here. Otherwise every detailed line estimator
        // would be zero for the first timestep after each restart. J_normfactor holds exactly the
        // factor that normalise_J() applies.
        prev_Jb_lu_normed[nonemptymgi][jblueindex].value =
            Jb_lu_raw[nonemptymgi][jblueindex].value * J_normfactor[nonemptymgi];
        prev_Jb_lu_normed[nonemptymgi][jblueindex].contribcount = Jb_lu_raw[nonemptymgi][jblueindex].contribcount;
      }
    }
  }
  assert_always(fscanf(gridsave_file, "%d\n", &code_check) == 1);
  assert_always(code_check == 42809403);
}

}  // namespace radfield
