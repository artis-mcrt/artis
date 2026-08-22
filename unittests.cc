// Unit tests for the pure numeric helpers (geometry, special relativity, sampling, decay chains,
// cross-sections, binning, input parsing, and atomic level structure). Build and run with:
//   make unittests && ./unittests
// The tests only cover functions with header-visible definitions or external linkage; they use no
// input files and no MPI communication, and a non-zero exit code means at least one check failed.
// (Compile-time checks of the constexpr helpers live in static_asserts next to their definitions.)

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <print>
#include <span>
#include <string_view>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "chargetransfer.h"
#include "constants.h"
#include "decay.h"
#include "exspec.h"
#include "gammapkt.h"
#include "globals.h"
#include "input.h"
#include "integrator.h"
#include "macroatom.h"
#include "mpi_logging.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "random.h"
#include "rpkt.h"
#include "sn3d.h"
#include "toms748.h"
#include "vectors.h"

namespace {

int checks_total = 0;
int checks_failed = 0;

void check(const bool pass, const std::string_view description) {
  checks_total++;
  if (!pass) {
    checks_failed++;
    std::println("  FAIL: {}", description);
  }
}

// check that a and b agree to within a relative tolerance
void check_close(const double a, const double b, const double reltol, const std::string_view description) {
  const bool pass = std::isfinite(a) && std::isfinite(b) &&
                    ((a == b) || (std::abs(a - b) <= (reltol * std::max(std::abs(a), std::abs(b)))));
  if (!pass) {
    std::println("  values: {:.15g} vs {:.15g} (reltol {:g})", a, b, reltol);
  }
  check(pass, description);
}

void test_binindex_helpers() {
  std::println("bin index and log-grid helpers...");
  const double minvalue = 1e14;
  const double dlog = 0.05;
  const ptrdiff_t nbins = 100;
  bool all_midpoints_match = true;
  bool edges_are_contiguous = true;
  for (ptrdiff_t i = 0; i < nbins; i++) {
    const double edge_low = get_loggrid_edge(minvalue, dlog, static_cast<double>(i));
    const double edge_high = get_loggrid_edge(minvalue, dlog, static_cast<double>(i + 1));
    edges_are_contiguous = edges_are_contiguous && (edge_high > edge_low);
    // the geometric midpoint of a log-spaced bin must map back to that bin
    const double midvalue = std::sqrt(edge_low * edge_high);
    all_midpoints_match = all_midpoints_match && (get_logbinindex(midvalue, minvalue, dlog, nbins) == i);
  }
  check(edges_are_contiguous, "get_loggrid_edge produces increasing bin edges");
  check(all_midpoints_match, "get_logbinindex returns the bin containing its geometric midpoint");

  check(get_logbinindex(minvalue / 10., minvalue, dlog, nbins) == 0, "get_logbinindex clamps below-range to bin 0");
  check(get_logbinindex(minvalue * 1e10, minvalue, dlog, nbins) == nbins - 1,
        "get_logbinindex clamps above-range to the last bin");
}

void test_range_chunks() {
  std::println("get_range_chunk / get_chunk_count...");
  rngstate_type rngstate{20260729};
  bool all_partitions_valid = true;
  for (int trial = 0; trial < 200; trial++) {
    const auto size = static_cast<ptrdiff_t>(rng_uniform(rngstate) * 10000);
    const auto nchunks = 1 + static_cast<ptrdiff_t>(rng_uniform(rngstate) * 32);
    ptrdiff_t expected_next_start = 0;
    ptrdiff_t maxchunksize = 0;
    ptrdiff_t minchunksize = size + 1;
    for (ptrdiff_t nchunk = 0; nchunk < nchunks; nchunk++) {
      const auto [nstart, nsize] = get_range_chunk(size, nchunks, nchunk);
      all_partitions_valid = all_partitions_valid && (nstart == expected_next_start) && (nsize >= 0);
      expected_next_start = nstart + nsize;
      maxchunksize = std::max(maxchunksize, nsize);
      minchunksize = std::min(minchunksize, nsize);
    }
    all_partitions_valid =
        all_partitions_valid && (expected_next_start == size) && ((maxchunksize - minchunksize) <= 1);
  }
  check(all_partitions_valid, "get_range_chunk contiguously partitions any range into nearly-equal chunks");

  check(get_chunk_count(0, 5) == 0, "get_chunk_count of an empty range is zero");
  check(get_chunk_count(10, 5) == 2, "get_chunk_count with exact division");
  check(get_chunk_count(11, 5) == 3, "get_chunk_count rounds up");
}

void test_vector_geometry() {
  std::println("vector geometry and special relativity...");
  rngstate_type rngstate{81102};

  bool cross_orthogonal = true;
  bool lagrange_identity = true;
  for (int trial = 0; trial < 100; trial++) {
    const auto vec_a = get_rand_isotropic_unitvec(rngstate);
    const auto vec_b = vec_scale(get_rand_isotropic_unitvec(rngstate), 2.5);
    const auto crossprod = cross_prod(vec_a, vec_b);
    cross_orthogonal =
        cross_orthogonal && (std::abs(dot(crossprod, vec_a)) < 1e-12) && (std::abs(dot(crossprod, vec_b)) < 1e-12);
    // |a x b|^2 + (a.b)^2 = |a|^2 |b|^2
    const double lhs = dot(crossprod, crossprod) + pow2(dot(vec_a, vec_b));
    const double rhs = dot(vec_a, vec_a) * dot(vec_b, vec_b);
    lagrange_identity = lagrange_identity && (std::abs(lhs - rhs) <= (1e-12 * rhs));
  }
  check(cross_orthogonal, "cross product is orthogonal to both inputs");
  check(lagrange_identity, "cross and dot products satisfy Lagrange's identity");

  bool aberration_roundtrip = true;
  for (int trial = 0; trial < 100; trial++) {
    const auto dir_frame1 = get_rand_isotropic_unitvec(rngstate);
    const auto vel = vec_scale(get_rand_isotropic_unitvec(rngstate), rng_uniform(rngstate) * 0.3 * CLIGHT);
    const auto dir_frame2 = angle_ab(dir_frame1, vel);
    const auto dir_frame1_roundtrip = angle_ab(dir_frame2, vec_scale(vel, -1.));
    for (int d = 0; d < 3; d++) {
      aberration_roundtrip = aberration_roundtrip && (std::abs(dir_frame1_roundtrip[d] - dir_frame1[d]) < 1e-12);
    }
  }
  check(aberration_roundtrip, "angle_ab aberration with -v inverts aberration with +v");

  // check the Doppler factor against a direct implementation
  const auto pos_rf = Vec3d{1.1e14, -2.4e14, 0.8e14};
  const auto dir_rf = vec_norm(Vec3d{0.3, -0.1, 0.9});
  const double prop_time = 10. * DAY;
  const auto vel_rf = get_velocity(pos_rf, prop_time);
  double doppler_expected = 1. - (dot(dir_rf, vel_rf) / CLIGHT);
  if (USE_RELATIVISTIC_DOPPLER_SHIFT) {
    doppler_expected /= std::sqrt(1 - (dot(vel_rf, vel_rf) / CLIGHTSQUARED));
  }
  check_close(calculate_doppler_nucmf_on_nurf(pos_rf, dir_rf, prop_time), doppler_expected, 1e-14,
              "calculate_doppler_nucmf_on_nurf matches the direct formula");

  // moving a packet must preserve the frame-invariant ratio e_cmf / nu_cmf = e_rf / nu_rf
  // (an outward-moving packet, so that the monotonic-nu_cmf clamp in move_pkt_withtime stays inactive)
  auto pos = vec_scale(dir_rf, 2.0e14);
  double t_current = 10. * DAY;
  const double nu_rf = 3e15;
  const double e_rf = 4e-12;
  const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pos, dir_rf, t_current);
  double nu_cmf = nu_rf * dopplerfactor;
  double e_cmf = e_rf * dopplerfactor;
  move_pkt_withtime(pos, dir_rf, t_current, nu_rf, nu_cmf, e_rf, e_cmf, 3.0e13);
  check_close(e_cmf / nu_cmf, e_rf / nu_rf, 1e-12, "move_pkt_withtime preserves e_cmf / nu_cmf = e_rf / nu_rf");
}

void test_escapedirectionbin() {
  std::println("get_escapedirectionbin...");
  rngstate_type rngstate{5501};
  constexpr int ndirs = 200000;
  std::array<int, MABINS> bincounts{};
  bool all_in_range = true;
  for (int n = 0; n < ndirs; n++) {
    const auto dirbin = get_escapedirectionbin(get_rand_isotropic_unitvec(rngstate));
    if (dirbin >= 0 && dirbin < MABINS) {
      bincounts[dirbin]++;
    } else {
      all_in_range = false;
    }
  }
  check(all_in_range, "escape direction bins are all within [0, MABINS)");

  // every bin covers an equal solid angle, so isotropic directions should fill them equally
  constexpr double expected_count = static_cast<double>(ndirs) / MABINS;
  const double sixsigma = 6. * std::sqrt(expected_count * (1. - (1. / MABINS)));
  const auto [mincount, maxcount] = std::ranges::minmax(bincounts);
  check(std::abs(mincount - expected_count) < sixsigma && std::abs(maxcount - expected_count) < sixsigma,
        "isotropic directions fill all equal-solid-angle bins to within 6 sigma");
}

void test_frame_transform() {
  std::println("frame_transform...");
  rngstate_type rngstate{99001};
  bool roundtrip_ok = true;
  bool polarisation_invariant = true;
  for (int trial = 0; trial < 100; trial++) {
    const auto n_rf = get_rand_isotropic_unitvec(rngstate);
    const double pol_p = rng_uniform(rngstate) * 0.9;
    const double pol_angle = rng_uniform(rngstate) * 2. * PI;
    const double q0 = pol_p * std::cos(pol_angle);
    const double u0 = pol_p * std::sin(pol_angle);
    const auto vel = vec_scale(get_rand_isotropic_unitvec(rngstate), rng_uniform(rngstate) * 0.2 * CLIGHT);

    const auto [n_cmf, q_cmf, u_cmf] = frame_transform(n_rf, q0, u0, vel);
    polarisation_invariant = polarisation_invariant && (std::abs(std::sqrt(pow2(q_cmf) + pow2(u_cmf)) - pol_p) < 1e-10);

    const auto [n_rf2, q2, u2] = frame_transform(n_cmf, q_cmf, u_cmf, vec_scale(vel, -1.));
    for (int d = 0; d < 3; d++) {
      roundtrip_ok = roundtrip_ok && (std::abs(n_rf2[d] - n_rf[d]) < 1e-10);
    }
    roundtrip_ok = roundtrip_ok && (std::abs(q2 - q0) < 1e-8) && (std::abs(u2 - u0) < 1e-8);
  }
  check(polarisation_invariant, "frame_transform preserves the polarisation degree");
  check(roundtrip_ok, "frame_transform with -v inverts frame_transform with +v");
}

void test_random_sampling() {
  std::println("random sampling...");
  rngstate_type rngstate{31415};
  constexpr int nsamples = 1000000;
  double total = 0.;
  bool all_in_range = true;
  for (int n = 0; n < nsamples; n++) {
    const auto zrand = rng_uniform(rngstate);
    all_in_range = all_in_range && (zrand >= 0.) && (zrand < 1.);
    total += zrand;
  }
  check(all_in_range, "rng_uniform values are within [0, 1)");
  // standard error of the mean is 1/sqrt(12 N)
  check(std::abs((total / nsamples) - 0.5) < (6. / std::sqrt(12. * nsamples)), "rng_uniform mean is 0.5");

  double sum_mu = 0.;
  double sum_musquared = 0.;
  bool unit_vectors = true;
  for (int n = 0; n < nsamples; n++) {
    const auto dir = get_rand_isotropic_unitvec(rngstate);
    unit_vectors = unit_vectors && (std::abs(vec_len(dir) - 1.) < 1e-6);
    sum_mu += dir[2];
    sum_musquared += pow2(dir[2]);
  }
  check(unit_vectors, "get_rand_isotropic_unitvec returns unit vectors");
  check(std::abs(sum_mu / nsamples) < (6. / std::sqrt(3. * nsamples)), "isotropic <cos(theta)> is zero");
  check(std::abs((sum_musquared / nsamples) - (1. / 3.)) < 1e-3, "isotropic <cos^2(theta)> is 1/3");
}

void test_planck() {
  std::println("Planck function integrals...");
  const double temperature = 6000.;
  // the frequency range covers x = h nu / k T from ~1e-7 to past the series cutoff, i.e. effectively 0 to infinity
  const double total = radfield::calculate_planck_integral(temperature, 1e8, 1e25, false);
  // tolerance covers the (mixed-vintage) physical constants rather than the numerics
  check_close(total, STEBO * pow4(temperature) / PI, 1e-3,
              "full-range Planck integral matches the Stefan-Boltzmann law");

  const double nu_split = 5.879e10 * temperature;  // Wien peak
  const double part_a = radfield::calculate_planck_integral(temperature, 1e8, nu_split, false);
  const double part_b = radfield::calculate_planck_integral(temperature, nu_split, 1e25, false);
  check_close(part_a + part_b, total, 1e-12, "Planck integrals over adjacent ranges are additive");

  const double nu_bar = radfield::calculate_planck_integral(temperature, 1e8, 1e25, true) / total;
  // full-range intensity-weighted mean frequency: <nu> = (k T / h) * 4 zeta(5) / zeta(4) = 3.832229... k T / h,
  // the same constant that radfield::set_params_fullspec() uses to convert nu_bar into T_R
  constexpr double zeta4 = 1.08232323371114;
  constexpr double zeta5 = 1.03692775514337;
  check_close(nu_bar, KB * temperature / H * 4. * zeta5 / zeta4, 1e-6,
              "intensity-weighted mean frequency of the full Planck spectrum");
}

void test_bateman() {
  std::println("Bateman decay chain solutions...");
  const double lambda_a = 1. / (8.80 * DAY);  // Ni56-like
  const double lambda_b = 1. / (113.7 * DAY);  // Co56-like
  const double time = 30. * DAY;
  const double initabund = 2.5;

  check_close(decay::calculate_decaychain(initabund, std::array{lambda_a}, time, false),
              initabund * std::exp(-lambda_a * time), 1e-14, "single-nuclide chain reduces to exponential decay");

  // daughter abundance of a two-step chain: N_b(t) = N_a(0) lambda_a / (lambda_b - lambda_a) (e^-l_a t - e^-l_b t)
  const double nb_analytic =
      initabund * lambda_a / (lambda_b - lambda_a) * (std::exp(-lambda_a * time) - std::exp(-lambda_b * time));
  check_close(decay::calculate_decaychain(initabund, std::array{lambda_a, lambda_b}, time, false), nb_analytic, 1e-12,
              "two-step chain matches the analytic Bateman solution");

  // with a stable sink at the end of the chain, everything eventually accumulates there
  check_close(decay::calculate_decaychain(initabund, std::array{lambda_a, 0.}, 1000. * 8.80 * DAY, false), initabund,
              1e-12, "a stable daughter eventually accumulates the whole initial abundance");

  // The expansion-factor mode weights each decay by t'/t for the adiabatic loss since it happened,
  // so for a single decay into a stable sink it is the closed form of
  // integral_0^t lambda N_0 exp(-lambda t') (t'/t) dt'. It is positive for every lambda*t > 0:
  // adiabatic losses cannot remove more energy than the decay released.
  for (const double x : {0.5, 3.4, 1e4}) {
    const double weighted_decaycount = (-std::expm1(-x) / x) - std::exp(-x);
    check_close(decay::calculate_decaychain(initabund, std::array{lambda_a, 0.}, x / lambda_a, true),
                initabund * weighted_decaycount, 1e-12,
                "expansion factor is the energy-weighted decay count of a chain into a stable sink");
  }

  // Below x of about 1e-3 the closed form is itself a difference of two quantities near one and
  // loses accuracy as epsilon/x, so compare against its series x/2 - x^2/3 instead. What matters
  // here is that the value is there at all: the (1 + 1/x) exp(-x) - 1/x form this replaced was 0.6%
  // wrong at x = 1e-7 and underflowed to exactly zero below about 1e-8, which silently cost every
  // long-lived nuclide its contribution to the initial temperature.
  for (const double x : {1e-9, 1e-7, 1e-5}) {
    check_close(decay::calculate_decaychain(initabund, std::array{lambda_a, 0.}, x / lambda_a, true),
                initabund * ((x / 2.) - (x * x / 3.)), 1e-5, "expansion factor at small lambda*timediff");
  }

  // and it is finite at the endpoint, where both forms divide by x
  check_close(decay::calculate_decaychain(initabund, std::array{lambda_a, 0.}, 0., true), 0., 1e-12,
              "expansion factor is zero when no time has elapsed");
}

void test_compton() {
  std::println("Compton scattering cross sections...");

  // reference: the standard Klein-Nishina total cross section, sigma = 2 pi r_e^2 * {...}, with 2 pi r_e^2 = 3/4
  // sigma_T
  const auto sigma_kleinnishina_total = [](const double x) -> double {
    return 0.75 * SIGMA_T *
           ((((1. + x) / pow3(x)) * (((2. * x * (1. + x)) / (1. + (2. * x))) - std::log1p(2. * x))) +
            (std::log1p(2. * x) / (2. * x)) - ((1. + (3. * x)) / pow2(1. + (2. * x))));
  };

  bool total_matches = true;
  for (const double x : {0.05, 0.2, 1., 5.}) {
    const double sigma_partial_full = gammapkt::sigma_compton_partial(x, 1. + (2. * x));
    total_matches = total_matches && (std::abs(sigma_partial_full - sigma_kleinnishina_total(x)) <=
                                      (1e-10 * sigma_kleinnishina_total(x)));
  }
  check(total_matches, "sigma_compton_partial over the full energy-loss range matches the Klein-Nishina formula");

  bool inversion_ok = true;
  for (const double x : {0.05, 1., 5.}) {
    for (const double zrand : {0.1, 0.5, 0.9}) {
      const double f = gammapkt::choose_f(x, zrand);
      const double fraction = gammapkt::sigma_compton_partial(x, f) / gammapkt::sigma_compton_partial(x, 1. + (2. * x));
      inversion_ok = inversion_ok && (std::abs(fraction - zrand) < 2e-4);
    }
  }
  check(inversion_ok, "choose_f inverts the partial Compton cross section to the solver tolerance");

  check_close(gammapkt::meanf_sigma(THOMSON_LIMIT * (1. - 1e-6)), gammapkt::meanf_sigma(THOMSON_LIMIT * (1. + 1e-6)),
              1e-5, "meanf_sigma is continuous across the Taylor-series/closed-form crossover");
}

void test_rad_deexcitation() {
  std::println("radiative deexcitation rate coefficient...");
  const double epsilon_trans = 2. * EV;
  const float einstein_A = 1e7;
  const double upperstatweight = 3.;
  const double lowerstatweight = 1.;
  const double t_current = 20. * DAY;

  check(rad_deexcitation_ratecoeff(epsilon_trans, einstein_A, upperstatweight, lowerstatweight, 0., 0., t_current) ==
            einstein_A,
        "radiative deexcitation is the Einstein A coefficient in the optically thin limit");

  // in the self-absorbed case the rate is reduced by the Sobolev escape probability beta = (1 - e^-tau) / tau
  const double nnupper = 1e5;
  const double nnlower = 1e8;
  const double nu_trans = epsilon_trans / H;
  const double b_ul = CLIGHTSQUAREDOVERTWOH / pow3(nu_trans) * einstein_A;
  const double b_lu = upperstatweight / lowerstatweight * b_ul;
  const double tau_sobolev = ((b_lu * nnlower) - (b_ul * nnupper)) * HCLIGHTOVERFOURPI * t_current;
  check(tau_sobolev > 1., "the self-absorbed test point is optically thick");
  check_close(rad_deexcitation_ratecoeff(epsilon_trans, einstein_A, upperstatweight, lowerstatweight, nnupper, nnlower,
                                         t_current),
              einstein_A * (-std::expm1(-tau_sobolev)) / tau_sobolev, 1e-12,
              "radiative deexcitation applies the Sobolev escape probability");
}

void test_phixs_table_lookup() {
  std::println("photoionisation cross-section table lookup...");
  globals::NPHIXSPOINTS = 10;
  globals::NPHIXSNUINCREMENT = 0.1;
  last_phixs_nuovernuedge = 1.0 + (globals::NPHIXSNUINCREMENT * (globals::NPHIXSPOINTS - 1));

  std::vector<float> photoion_xs(globals::NPHIXSPOINTS);
  for (int i = 0; i < globals::NPHIXSPOINTS; i++) {
    photoion_xs[i] = static_cast<float>((i + 1) * 1e-18);
  }
  const double nu_edge = 3e15;

  check(photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_edge * 0.99) == 0.,
        "cross section is zero below the threshold");

  if constexpr (PHIXS_CLASSIC_NO_INTERPOLATION) {
    // classic mode compares nu == nu_edge directly, so the threshold value is exact in every build
    check(photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_edge) == photoion_xs[0],
          "classic mode cross section at the threshold is the first table point");
    check(photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_edge * (1. + 0.25)) == photoion_xs[2],
          "classic mode truncates to the nearest lower table point");
    // regression test for the former out-of-bounds read: scan frequencies approaching the upper limit of
    // the tabulated range from below (the last few representable values fall in the final table cell)
    bool tail_reads_in_table = true;
    const double nu_out_bound = nu_edge * (1 + (globals::NPHIXSNUINCREMENT * globals::NPHIXSPOINTS));
    double nu = nu_out_bound * (1. - 1e-13);
    while (nu < nu_out_bound) {
      const auto sigma = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu);
      tail_reads_in_table = tail_reads_in_table && (std::ranges::find(photoion_xs, sigma) != photoion_xs.end());
      nu = std::nextafter(nu, nu_out_bound);
    }
    check(tail_reads_in_table, "classic mode reads within the table for every nu below the tabulated range's end");
  } else {
    // probe just above the threshold rather than exactly at it: the fast-math production builds can
    // round the interpolation index at exactly nu_edge to either side of zero
    check_close(photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_edge * (1. + 1e-9)), photoion_xs[0],
                1e-6, "cross section just above the threshold is the first table point");
    check_close(photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_edge * (1. + 0.05)),
                0.5 * (photoion_xs[0] + photoion_xs[1]), 1e-6, "cross section is interpolated between table points");
  }

  // past the tabulated range, the Kramers nu^-3 extrapolation continues from the last table point.
  // classic mode anchors the extrapolation at nu_edge * (1 + increment * npoints) instead of the
  // frequency of the last table point
  const double nu_kramers_anchor = PHIXS_CLASSIC_NO_INTERPOLATION
                                       ? nu_edge * (1. + (globals::NPHIXSNUINCREMENT * globals::NPHIXSPOINTS))
                                       : nu_edge * last_phixs_nuovernuedge;
  const double nu_above_table = 4. * nu_edge;
  check_close(photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu_above_table),
              photoion_xs[globals::NPHIXSPOINTS - 1] * pow3(nu_kramers_anchor / nu_above_table), 1e-6,
              "cross section follows the Kramers extrapolation above the tabulated range");
}

void test_closest_transition_randomised() {
  std::println("closest_transition against a reference implementation...");
  rngstate_type rngstate{777};

  // build a descending-frequency linelist
  std::vector<double> linelist_nu(500);
  for (auto& nu : linelist_nu) {
    nu = 1e14 + (rng_uniform(rngstate) * 1e15);
  }
  std::ranges::sort(linelist_nu, std::ranges::greater{});

  bool all_match = true;
  for (int trial = 0; trial < 1000; trial++) {
    const double nu_cmf = 0.5e14 + (rng_uniform(rngstate) * 1.2e15);
    // reference: the first (bluest) line with nu_line <= nu_cmf, or -1 if none
    int expected = -1;
    for (int i = 0; i < std::ssize(linelist_nu); i++) {
      if (linelist_nu[i] <= nu_cmf) {
        expected = i;
        break;
      }
    }
    all_match = all_match && (closest_transition(nu_cmf, -1, linelist_nu) == expected);
  }
  check(all_match, "closest_transition matches a linear-scan reference for random frequencies");
}

void test_input_helpers() {
  std::println("input parsing helpers...");
  check(lineiscommentonly(""), "empty line is comment-only");
  check(lineiscommentonly("   \t"), "whitespace line is comment-only");
  check(lineiscommentonly("# a comment"), "hash line is comment-only");
  check(lineiscommentonly("   # indented comment"), "indented hash line is comment-only");
  check(!lineiscommentonly("1 2 3"), "data line is not comment-only");
  check(!lineiscommentonly("1 2 3 # trailing comment"), "data line with trailing comment is not comment-only");
}

void test_parse_next_token() {
  std::println("parse_next_token...");
  {
    auto remainder = std::string_view{"  42\t-7 +5  "};
    int i = 0;
    check(parse_next_token(remainder, i) && i == 42, "int token with leading whitespace");
    check(parse_next_token(remainder, i) && i == -7, "negative int token after tab");
    check(parse_next_token(remainder, i) && i == 5, "leading plus sign accepted");
    check(!parse_next_token(remainder, i) && i == 5, "trailing whitespace has no token and value is untouched");
  }
  {
    auto remainder = std::string_view{};
    int i = -99;
    check(!parse_next_token(remainder, i), "empty line has no token");
  }
  {
    auto remainder = std::string_view{"12abc 3"};
    int i = -99;
    check(!parse_next_token(remainder, i) && remainder == "12abc 3",
          "token with trailing junk is rejected and not consumed");
  }
  {
    auto remainder = std::string_view{"+-0"};
    int i = -99;
    check(!parse_next_token(remainder, i), "+-0 is rejected");
  }
  {
    auto remainder = std::string_view{"3000000000"};
    int i = -99;
    check(!parse_next_token(remainder, i), "int overflow is rejected");
  }
  {
    auto remainder = std::string_view{"1.5e3 -2.5 +0.25"};
    double d = 0.;
    check(parse_next_token(remainder, d) && d == 1500., "double token in exponent notation");
    check(parse_next_token(remainder, d) && d == -2.5, "negative double token");
    check(parse_next_token(remainder, d) && d == 0.25, "double token with leading plus sign");
  }
  for (const auto* const token : {"nan", "inf", "-inf", "NAN"}) {
    auto remainder = std::string_view{token};
    double d = 0.;
    check(!parse_next_token(remainder, d), "nan and inf spellings are rejected");
  }
  {
    auto remainder = std::string_view{"1e400"};
    double d = -99.;
    check(!parse_next_token(remainder, d), "double overflow is rejected");
  }
  {
    auto remainder = std::string_view{"1e-400"};
    double d = -99.;
    check(parse_next_token(remainder, d) && d == 0., "underflow below the double range reads as zero");
  }
  {
    auto remainder = std::string_view{"1e39"};
    float f = -99.;
    check(!parse_next_token(remainder, f), "float overflow is rejected even when the value fits in a double");
  }
  {
    auto remainder = std::string_view{"1e-60"};
    float f = -99.;
    check(parse_next_token(remainder, f) && f == 0.F, "underflow below the float range reads as zero");
  }
}

void test_count_groundterm_levels() {
  std::println("ground term level count...");
  // The lowest levels of real ions, as (energies [eV], statistical weights). The static_asserts in input.h
  // show the LS coupling rules. These tests add ions with more levels.
  const auto count = [](const std::vector<double>& energies, const std::vector<float>& statweights) {
    return count_groundterm_levels(energies, statweights);
  };
  check(count({}, {}) == 0, "no levels");
  check(count({0.}, {1.F}) == 1, "one level");
  check(count({0., 0.1}, {2.F, 4.F}) == 2, "two-level 2P1/2,3/2 with no third level to compare with");
  check(count({0., 0.1}, {2.F, 2.F}) == 1, "two levels with the same statistical weight");

  // Fe I 5D4..0 (inverted) followed by 5F5 (g = 11)
  check(count({0., 0.0516, 0.0873, 0.1101, 0.1213, 0.8590, 0.9146}, {9.F, 7.F, 5.F, 3.F, 1.F, 11.F, 9.F}) == 5,
        "Fe I 5D4..0");
  // Fe III 5D4..0 (inverted) followed by 3H6 (g = 13)
  check(count({0., 0.0541, 0.0916, 0.1156, 0.1274, 2.4059, 2.4860}, {9.F, 7.F, 5.F, 3.F, 1.F, 5.F, 13.F}) == 5,
        "Fe III 5D4..0");
  // Co I 4F9/2..3/2 (inverted), then 3d8 4s2 4F9/2 (g = 10) only 0.2 eV above it
  check(count({0., 0.1012, 0.1744, 0.2243, 0.4318, 0.5136}, {10.F, 8.F, 6.F, 4.F, 10.F, 8.F}) == 4, "Co I 4F9/2..3/2");
  // Ti II 4F3/2..9/2 (normal) followed by b4F3/2 (g = 4)
  check(count({0., 0.0117, 0.0280, 0.0488, 0.1126, 0.1220}, {4.F, 6.F, 8.F, 10.F, 4.F, 6.F}) == 4, "Ti II 4F3/2..9/2");
  // Ni II 2D5/2,3/2 (inverted) followed by 4F9/2 (g = 10)
  check(count({0., 0.1868, 1.0407, 1.1568}, {6.F, 4.F, 10.F, 8.F}) == 2, "Ni II 2D5/2,3/2");
  // Sc II 3D1,2,3 (normal) followed by 1D2 (g = 5)
  check(count({0., 0.0084, 0.0220, 0.3150, 0.5955}, {3.F, 5.F, 7.F, 5.F, 5.F}) == 3, "Sc II 3D1,2,3");
  // Cr I 7S3, then 5S2, then 5D0..4. The weight changes correctly from 7S3 to 5S2, but the splitting is large
  check(count({0., 0.9414, 0.9610, 0.9684, 0.9829}, {7.F, 5.F, 1.F, 3.F, 5.F}) == 1, "Cr I 7S3");
  // He I 1S0 then 2 3S1 (g = 3) then 2 1S0
  check(count({0., 19.8196, 20.6158, 20.9641}, {1.F, 3.F, 1.F, 5.F}) == 1, "He I 1S0");
  // O II 4S3/2, then 2D5/2,3/2 (almost the same energy), then 2P3/2,1/2
  check(count({0., 3.3241, 3.3266, 5.0173, 5.0174}, {4.F, 6.F, 4.F, 4.F, 2.F}) == 1, "O II 4S3/2");
  // Ca V 3P2,1,0 (inverted, splitting ratio 2.8) then 1D2
  check(count({0., 0.2981, 0.4061, 2.3347, 5.4350}, {5.F, 3.F, 1.F, 5.F, 1.F}) == 3, "Ca V 3P2,1,0");
  // Ar V 3P0,1,2 (normal) then 1D2
  check(count({0., 0.0948, 0.2517, 2.0209, 4.7007}, {1.F, 3.F, 5.F, 5.F, 1.F}) == 3, "Ar V 3P0,1,2");
  // Fe XXI 3P0,1,2. At this charge the interval rule is inaccurate: the second splitting is 0.29 of
  // the prediction and not 2. Only the larger tolerance keeps this term together
  check(count({0., 9.1524, 14.5438, 30.3089, 46.0904}, {1.F, 3.F, 5.F, 5.F, 1.F}) == 3, "Fe XXI 3P0,1,2");
  // Fe XI 3P2,1,0 (inverted, the same condition: error factor 3.85), then 1D2
  check(count({0., 1.5700, 1.7737, 4.6776, 10.0156}, {5.F, 3.F, 1.F, 5.F, 1.F}) == 3, "Fe XI 3P2,1,0");
  // Mg IV data gives 2P1/2 and 2P3/2 both at zero energy. The two levels stay in one term
  check(count({0., 0., 0.2762, 38.6251}, {2.F, 4.F, 2.F, 2.F}) == 2, "degenerate 2P1/2,3/2 pair");
}

void test_calculate_timesteps() {
  std::println("timestep grid calculation...");
  const double tmin = 2. * DAY;
  const double tmax = 50. * DAY;
  const int nts = 40;

  const auto grid_is_valid = [](const std::vector<globals::TimeStep>& ts, const int ntimesteps, const double t_start,
                                const double t_end) {
    bool ok = std::ssize(ts) == ntimesteps + 1;
    ok = ok && ts[0].start == t_start;
    ok = ok && ts[ntimesteps].start == t_end && ts[ntimesteps].mid == t_end && ts[ntimesteps].width == 0.;
    for (int n = 0; ok && n < ntimesteps; n++) {
      ok = ok && ts[n].width > 0.;
      ok = ok && ts[n].mid > ts[n].start && ts[n].mid < ts[n].start + ts[n].width;
      const double ts_end = ts[n].start + ts[n].width;
      const double next_start = (n < ntimesteps - 1) ? ts[n + 1].start : t_end;
      ok = ok && std::abs((ts_end / next_start) - 1.) < 0.001;
    }
    return ok;
  };

  {
    const auto ts = calculate_timesteps(TimeStepSizeMethod::LOGARITHMIC, tmin, tmax, nts, -1., -1.);
    check(grid_is_valid(ts, nts, tmin, tmax), "logarithmic grid is contiguous from tmin to tmax");
    bool constant_ratio = true;
    for (int n = 1; n < nts; n++) {
      constant_ratio =
          constant_ratio && std::abs(((ts[n].start / ts[n - 1].start) / (ts[1].start / ts[0].start)) - 1.) < 1e-10;
    }
    check(constant_ratio, "logarithmic timesteps have a constant start-time ratio");
  }
  {
    const auto ts = calculate_timesteps(TimeStepSizeMethod::CONSTANT, tmin, tmax, nts, -1., -1.);
    check(grid_is_valid(ts, nts, tmin, tmax), "constant grid is contiguous from tmin to tmax");
    bool equal_widths = true;
    for (int n = 0; n < nts; n++) {
      equal_widths = equal_widths && std::abs((ts[n].width / ts[0].width) - 1.) < 1e-10;
    }
    check(equal_widths, "constant timesteps have equal widths");
  }
  {
    // 2 to 30 days logarithmic, then 1-day fixed steps to 50 days
    const auto ts = calculate_timesteps(TimeStepSizeMethod::LOGARITHMIC_THEN_CONSTANT, tmin, tmax, nts, 1., 30.);
    check(grid_is_valid(ts, nts, tmin, tmax), "logarithmic-then-constant grid is contiguous from tmin to tmax");
    // ceil((50 - 30) / 1) = 20 fixed steps at the end
    bool fixed_tail = true;
    for (int n = nts - 20; n < nts; n++) {
      fixed_tail = fixed_tail && std::abs((ts[n].width / DAY) - 1.) < 1e-10;
    }
    check(fixed_tail, "logarithmic-then-constant grid ends with the fixed-width steps");
  }
  {
    // 1-day fixed steps from 2 to 30 days, then logarithmic to 50 days
    const auto ts = calculate_timesteps(TimeStepSizeMethod::CONSTANT_THEN_LOGARITHMIC, tmin, tmax, nts, 1., 30.);
    check(grid_is_valid(ts, nts, tmin, tmax), "constant-then-logarithmic grid is contiguous from tmin to tmax");
    // ceil((30 - 2) / 1) = 28 fixed steps at the start
    bool fixed_head = true;
    for (int n = 0; n < 28; n++) {
      fixed_head = fixed_head && std::abs((ts[n].width / DAY) - 1.) < 1e-10;
    }
    check(fixed_head, "constant-then-logarithmic grid starts with the fixed-width steps");
  }
  {
    const auto ts = calculate_timesteps(TimeStepSizeMethod::LOGARITHMIC, tmin, tmax, 1, -1., -1.);
    check(grid_is_valid(ts, 1, tmin, tmax), "a single logarithmic timestep covers the whole time range");
  }
}

void test_rank_outfile_name() {
  std::println("per-rank output filename matching...");
  bool match_all = true;
  for (const auto* const name : {
           "output_0-0.txt",
           "output_123-45.txt",
           "estimators_0000.out",
           "nlte_0003.out",
           "radfield_0100.out",
           "macroatom_0000.out",
           "output_0-0.txt.zst",
           "estimators_0000.out.gz",
           "radfield_0000.out.xz",
       }) {
    match_all = match_all && is_rank_outfile_name(name);
  }
  check(match_all, "generated per-rank filenames are matched, with and without compression extensions");

  bool match_none = false;
  for (const auto* const name : {
           "output_0.txt",
           "output_a-0.txt",
           "output_-0.txt",
           "output_0-.txt",
           "output_0-0.log",
           "estimators_00x0.out",
           "estimators_.out",
           "estimators_0000.dat",
           "nlte_0000.OUT",
           "spec.out",
           "packets00_0000.out",
           "deposition.out",
           "model.txt",
           "input.txt",
           "output_0-0.txt.bz2",
           "estimators_0000.out.zst.zst",
       }) {
    match_none = match_none || is_rank_outfile_name(name);
  }
  check(!match_none, "other filenames are never matched, so they cannot be deleted from an output folder");
}

void test_gth_solver() {
  std::println("GTH stationary distribution solver...");

  // build matrices in the NLTE solver layout A[(to * n) + from]: each off-diagonal holds a transition rate and each
  // diagonal element holds the negated column sum, as the NLTE rate matrix assembly does (GTH must never read it)
  const auto set_rate = [](std::vector<double>& matrix, const ptrdiff_t n, const ptrdiff_t from, const ptrdiff_t to,
                           const double rate) {
    matrix[(to * n) + from] += rate;
    matrix[(from * n) + from] -= rate;
  };

  // three-state chain with the exact stationary distribution {53, 23, 18} / 94 from the Markov chain tree theorem
  // (each state's weight is the sum over its rooted spanning trees of the product of the tree's rates)
  const auto make_threestate_matrix = [&set_rate] {
    std::vector<double> matrix(3 * 3, 0.);
    set_rate(matrix, 3, 0, 1, 1.);
    set_rate(matrix, 3, 0, 2, 2.);
    set_rate(matrix, 3, 1, 0, 3.);
    set_rate(matrix, 3, 1, 2, 4.);
    set_rate(matrix, 3, 2, 0, 5.);
    set_rate(matrix, 3, 2, 1, 6.);
    return matrix;
  };

  const auto check_threestate_result = [](const std::vector<double>& vec_x, const std::string_view description) {
    constexpr std::array expected{53. / 94., 23. / 94., 18. / 94.};
    for (ptrdiff_t k = 0; k < std::ssize(expected); k++) {
      check_close(vec_x[k], expected[k], 1e-14, description);
    }
  };

  // A birth-death chain of adjacent states. rate_up[k] is the rate from state k to k + 1, and rate_down[k]
  // is the rate in the opposite direction. The function fills the matrix of the caller and does not return
  // one, because gcc 16.0.1 with LTO gives an incorrect maybe-uninitialized warning for a returned vector
  // whose size comes from the data. make_threestate_matrix below returns a vector and gives no warning.
  // Return a vector here again when the gcc version in CI changes.
  const auto make_chain_matrix = [&set_rate](std::vector<double>& matrix, const std::span<const double> rate_up,
                                             const std::span<const double> rate_down) {
    const auto n = std::ssize(rate_up) + 1;
    matrix.assign(n * n, 0.);
    for (ptrdiff_t k = 0; k < n - 1; k++) {
      set_rate(matrix, n, k, k + 1, rate_up[k]);
      set_rate(matrix, n, k + 1, k, rate_down[k]);
    }
  };

  {
    auto matrix = make_threestate_matrix();
    std::vector<double> vec_x(3, 0.);
    const auto result = gth_stationary_distribution(matrix, vec_x);
    check(!result.has_value(), "GTH solves an irreducible three-state chain");
    check_threestate_result(vec_x, "GTH three-state stationary distribution");
    check_close(vec_x[0] + vec_x[1] + vec_x[2], 1., 1e-14, "GTH stationary distribution sums to one");
  }

  {
    // corrupting the diagonal cannot change the result because GTH never reads it
    auto matrix = make_threestate_matrix();
    matrix[(0 * 3) + 0] = 12345.;
    matrix[(1 * 3) + 1] = -67890.;
    matrix[(2 * 3) + 2] = 0.5;
    std::vector<double> vec_x(3, 0.);
    const auto result = gth_stationary_distribution(matrix, vec_x);
    check(!result.has_value(), "GTH solves with a corrupted diagonal");
    check_threestate_result(vec_x, "GTH ignores the diagonal");
  }

  {
    // birth-death chain with extreme dynamic range: GTH must reproduce the exact detailed-balance ratios of
    // neighbouring populations (the componentwise relative accuracy that motivates the algorithm, where an LU solve
    // with a normalisation row loses the small components to absolute rounding)
    constexpr std::array<double, 3> rate_up = {2e10, 3e-4, 5e2};
    constexpr std::array<double, 3> rate_down = {7e-6, 1e8, 4e-1};
    std::vector<double> matrix;
    make_chain_matrix(matrix, rate_up, rate_down);
    std::vector<double> vec_x(rate_up.size() + 1, 0.);
    const auto result = gth_stationary_distribution(matrix, vec_x);
    check(!result.has_value(), "GTH solves a birth-death chain");
    for (ptrdiff_t k = 0; k < std::ssize(rate_up); k++) {
      check_close(vec_x[k + 1] / vec_x[k], rate_up[k] / rate_down[k], 1e-13,
                  "GTH birth-death chain preserves a detailed-balance population ratio");
    }
  }

  {
    // seven-state chain whose raw stationary weights span 1e360 relative to state zero, which would overflow the
    // double range without the subtraction-free rescaling in the back-substitution
    const std::vector<double> rate_up(6, 1e30);
    const std::vector<double> rate_down(6, 1e-30);
    std::vector<double> matrix;
    make_chain_matrix(matrix, rate_up, rate_down);
    std::vector<double> vec_x(7, 0.);
    const auto result = gth_stationary_distribution(matrix, vec_x);
    check(!result.has_value(), "GTH survives raw weights beyond the double range");
    check(std::ranges::all_of(vec_x, [](const double x) { return std::isfinite(x) && x >= 0.; }),
          "GTH extreme-ratio distribution is finite and non-negative");
    check_close(vec_x[6], 1., 1e-12, "GTH extreme-ratio distribution is dominated by the top state");
    check_close(vec_x[5] / vec_x[6], 1e-60, 1e-12, "GTH extreme-ratio detailed-balance ratio survives the rescaling");
  }

  {
    // a single up/down rate ratio of 1e400 exceeds the double range outright: the pre-division rescaling must keep
    // the dominant state finite instead of letting the weight overflow to infinity and decay into NaNs
    std::vector<double> matrix;
    make_chain_matrix(matrix, std::array{1e200}, std::array{1e-200});
    std::vector<double> vec_x(2, 0.);
    const auto result = gth_stationary_distribution(matrix, vec_x);
    check(!result.has_value(), "GTH solves a two-state chain with a weight ratio beyond the double range");
    check(std::ranges::all_of(vec_x, [](const double x) { return std::isfinite(x) && x >= 0.; }),
          "GTH beyond-range-ratio distribution is finite and non-negative");
    check_close(vec_x[1], 1., 1e-12, "GTH beyond-range-ratio distribution is dominated by the top state");
  }

  {
    // state 2 has no departure rate: the chain is reducible, so the solve must fail, identify state 2, and leave the
    // output vector untouched
    std::vector<double> matrix(3 * 3, 0.);
    set_rate(matrix, 3, 0, 1, 1.);
    set_rate(matrix, 3, 1, 0, 2.);
    set_rate(matrix, 3, 1, 2, 3.);
    std::vector<double> vec_x(3, -42.);
    const auto result = gth_stationary_distribution(matrix, vec_x);
    check(result.has_value() && result.value() == 2, "GTH reports a state with no departure rate");
    check(std::ranges::all_of(vec_x, [](const double x) { return x == -42.; }),
          "GTH leaves the output vector untouched on failure");
  }

  {
    // two decoupled two-state blocks: eliminating state 3 succeeds (it connects to state 2), then state 2 has no
    // rate into states 0 or 1 and the solve must fail there
    std::vector<double> matrix(4 * 4, 0.);
    set_rate(matrix, 4, 0, 1, 1.);
    set_rate(matrix, 4, 1, 0, 2.);
    set_rate(matrix, 4, 2, 3, 3.);
    set_rate(matrix, 4, 3, 2, 4.);
    std::vector<double> vec_x(4, 0.);
    const auto result = gth_stationary_distribution(matrix, vec_x);
    check(result.has_value() && result.value() == 2, "GTH detects decoupled blocks of states");
  }

  {
    // an infinite rate reaching a weight that the overflow rescaling has already driven to zero gives an inflow of
    // 0 * inf = NaN, which no further rescaling can clear. The rescue must give up after a bounded number of
    // passes and hand the non-finite result to the caller's population check: an unbounded loop hangs here rather
    // than returning a wrong answer, so reaching the checks below at all is most of what this pins down.
    // States 0 and 1 carry the beyond-range ratio that forces a rescale, and the rate from 0 to 2 is infinite.
    // Row 2 is never rewritten by the elimination (it is eliminated first), so that infinity survives into the
    // back-substitution for state 2.
    std::vector<double> matrix(3 * 3, 0.);
    set_rate(matrix, 3, 0, 1, 1e200);
    set_rate(matrix, 3, 1, 0, 1e-200);
    set_rate(matrix, 3, 2, 0, 1.);
    set_rate(matrix, 3, 0, 2, std::numeric_limits<double>::infinity());
    std::vector<double> vec_x(3, 0.);
    const auto result = gth_stationary_distribution(matrix, vec_x);
    check(!result.has_value(), "GTH terminates on an unrescuable non-finite inflow");
    check(!std::isfinite(vec_x[2]), "GTH passes the unrescuable inflow through as a non-finite population");
  }

  {
    // a single state has stationary distribution {1}, and the value in the (never-read) diagonal is irrelevant
    std::vector<double> matrix{12345.};
    std::vector<double> vec_x(1, 0.);
    const auto result = gth_stationary_distribution(matrix, vec_x);
    check(!result.has_value(), "GTH solves the single-state chain");
    check_close(vec_x[0], 1., 1e-15, "GTH single-state distribution is {1}");
  }
}

void test_nonthermal_solve_upper_triangular() {
  std::println("Spencer-Fano compacted upper triangular solver...");
  constexpr int sfpts = nonthermal::SFPTS;

  // deterministic 64-bit LCG (Knuth MMIX constants) giving uniform doubles in [0, 1)
  uint64_t lcgstate = 20260808;
  const auto next_uniform = [&lcgstate] {
    lcgstate = (lcgstate * 6364136223846793005ULL) + 1442695040888963407ULL;
    return std::ldexp(static_cast<double>(lcgstate >> 11U), -53);
  };

  // element (i, j) of the compacted upper triangular matrix (same layout as nonthermal.cc's uppertriangular())
  const auto uppertriangular = [](const int i, const int j) { return (sfpts * i) - (i * (i + 1) / 2) + j; };

  // diagonally dominant upper triangular system: diagonal in [1, 2), off-diagonal magnitude < 2^-13
  std::vector<double> uppertri(static_cast<size_t>(sfpts) * (sfpts + 1) / 2);
  for (int i = 0; i < sfpts; i++) {
    uppertri[uppertriangular(i, i)] = 1. + next_uniform();
    for (int j = i + 1; j < sfpts; j++) {
      uppertri[uppertriangular(i, j)] = std::ldexp(next_uniform() - 0.5, -12);
    }
  }
  std::array<double, sfpts> bvec{};
  for (auto& b : bvec) {
    b = 2. * (next_uniform() - 0.5);
  }
  // exact zeros in the last rows (solved first) exercise the skip-division-on-exact-zero path
  for (int i = sfpts - 8; i < sfpts; i++) {
    bvec[i] = 0.;
  }

  std::array<double, sfpts> xvec{};
  nonthermal::solve_upper_triangular(uppertri, bvec, xvec);

  // back substitution on this diagonally dominant system is accurate to a small multiple of machine epsilon
  double maxresidual = 0.;
  for (int i = 0; i < sfpts; i++) {
    double matvecprod = 0.;
    for (int j = i; j < sfpts; j++) {
      matvecprod += uppertri[uppertriangular(i, j)] * xvec[j];
    }
    maxresidual = std::max(maxresidual, std::fabs(bvec[i] - matvecprod));
  }
  check(maxresidual < 1e-11, "small residual");
  check(std::bit_cast<uint64_t>(xvec[sfpts - 1]) == std::bit_cast<uint64_t>(0.),
        "exactly-zero rhs rows give exactly +0. solution values");
}

// the TOMS 748 root finder and Gauss-Kronrod quadrature extracted from Boost.Math
void test_chargetransfer_helpers() {
  // fit form of Kingdon & Ferland (1996): k = a * (T/1e4)^b * (1 + c * exp(d * T/1e4)) * exp(-eexp/T)
  check(chargetransfer::evaluate_ctfit(1e-9, 0., 0., 0., 0., 1e3, 1e5, 1e4) == 1e-9,
        "evaluate_ctfit gives the coefficient a at T = 1e4 K for a flat fit");
  check_close(chargetransfer::evaluate_ctfit(1e-9, 1., 0., 0., 0., 1e3, 1e5, 2e4), 2e-9, 1e-12,
              "evaluate_ctfit applies the power law term");
  check_close(chargetransfer::evaluate_ctfit(1e-9, 1., 0., 0., 0., 1e3, 1e5, 1e7), 10e-9, 1e-12,
              "evaluate_ctfit clamps the fit temperature at tmax");
  check_close(chargetransfer::evaluate_ctfit(1e-9, 0., 1., -1., 0., 1e3, 1e5, 1e4), 1e-9 * (1. + std::exp(-1.)), 1e-12,
              "evaluate_ctfit applies the exponential correction term");
  check_close(chargetransfer::evaluate_ctfit(1e-9, 0., 0., 0., 2e4, 1e3, 1e5, 1e4), 1e-9 * std::exp(-2.), 1e-12,
              "evaluate_ctfit applies the Boltzmann factor of an endothermic fit");
  check(chargetransfer::evaluate_ctfit(1e-9, 0., -2., 0., 0., 1e3, 1e5, 1e4) == 0.,
        "evaluate_ctfit floors a negative fit value at zero");

  // Landau-Zener cross section for capture by a doubly charged ion from a neutral donor
  const double deltae = 2. * EV;
  const double ip_donor = 13. * EV;
  const double sigma_mid = chargetransfer::sigma_lz_channel(2, deltae, ip_donor, 1e7);
  check(sigma_mid > 0., "sigma_lz_channel is positive for a favourable energy defect");

  // the cross section stays below the geometric limit pi * rx^2 (the transfer probability tops out at 0.5)
  const double rx_cm = 5.29177211e-9 * (27.211386 / 2.);
  check(sigma_mid < PI * rx_cm * rx_cm, "sigma_lz_channel stays below the geometric limit");

  // the crossing is traversed adiabatically at a very low velocity and diabatically at a high
  // velocity, so the cross section vanishes at both ends
  check(chargetransfer::sigma_lz_channel(2, deltae, ip_donor, 1e2) < 1e-3 * sigma_mid,
        "sigma_lz_channel vanishes in the adiabatic low velocity limit");
  check(chargetransfer::sigma_lz_channel(2, deltae, ip_donor, 5e8) < sigma_mid,
        "sigma_lz_channel decreases in the diabatic high velocity limit");

  // a large energy defect moves the crossing inward and suppresses the rate
  check(chargetransfer::sigma_lz_channel(2, 15. * EV, ip_donor, 1e7) < sigma_mid,
        "sigma_lz_channel is smaller for a large energy defect");

  // the near-resonant estimate: an acceptor charge of one to three and an energy defect in (0, 4 eV]
  check(chargetransfer::is_near_resonant(1, 2. * EV), "is_near_resonant accepts a singly charged ion");
  check(chargetransfer::is_near_resonant(3, 4. * EV), "is_near_resonant accepts the upper edges of charge and defect");
  check(!chargetransfer::is_near_resonant(4, 2. * EV), "is_near_resonant rejects a charge above three");
  check(!chargetransfer::is_near_resonant(2, 4.5 * EV), "is_near_resonant rejects a large energy defect");
  check(!chargetransfer::is_near_resonant(2, 0.), "is_near_resonant rejects a reaction with no energy release");
  check(!chargetransfer::is_near_resonant(2, -1. * EV), "is_near_resonant rejects an endothermic reaction");
}

void test_toms748_and_gauss_kronrod() {
  const auto f_cosfixedpoint = [](const double x) { return std::cos(x) - x; };
  std::uintmax_t iterations = 50;
  const auto [rootlow, roothigh] =
      toms748_solve(f_cosfixedpoint, 0., 1., f_cosfixedpoint(0.), f_cosfixedpoint(1.), ftol<1e-12>, iterations);
  const double root = 0.5 * (rootlow + roothigh);
  check(std::fabs(root - 0.7390851332151607) < 1e-10, "toms748_solve finds the fixed point of cos(x)");
  check(iterations < 50, "toms748_solve converges in fewer than the maximum iterations");

  double abserr{NAN};
  const double integral_sin =
      integrator([](const double x) { return std::sin(x); }, 0., std::numbers::pi, 1e-10, &abserr);
  check(std::fabs(integral_sin - 2.) < 1e-12, "gauss_kronrod_integrate of sin over [0, pi] is 2");
  check(abserr < 1e-10, "gauss_kronrod_integrate error estimate is small for a smooth integrand");

  // sqrt(pi / 100) / 2 * (erf(117) + erf(3))
  const double integral_bump =
      integrator<31>([](const double x) { return std::exp(-100. * (x - 0.3) * (x - 0.3)); }, 0., 12., 1e-10, &abserr);
  check(std::fabs(integral_bump - 0.17724342737116647) < 1e-9,
        "31-point adaptive gauss_kronrod_integrate resolves a narrow Gaussian bump");
}

}  // anonymous namespace

auto main() -> int {
  test_binindex_helpers();
  test_range_chunks();
  test_vector_geometry();
  test_escapedirectionbin();
  test_frame_transform();
  test_random_sampling();
  test_planck();
  test_bateman();
  test_compton();
  test_rad_deexcitation();
  test_phixs_table_lookup();
  test_closest_transition_randomised();
  test_input_helpers();
  test_parse_next_token();
  test_count_groundterm_levels();
  test_calculate_timesteps();
  test_rank_outfile_name();
  test_gth_solver();
  test_chargetransfer_helpers();
  test_nonthermal_solve_upper_triangular();
  // the solver/integrator throw std::domain_error only on precondition violations that this test does not trigger
  // cppcheck-suppress throwInEntryPoint
  test_toms748_and_gauss_kronrod();

  std::println("unit tests: {} of {} checks passed", checks_total - checks_failed, checks_total);

  return (checks_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
