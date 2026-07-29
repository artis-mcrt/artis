// Unit tests for the pure numeric helpers (geometry, special relativity, sampling, decay chains,
// cross-sections, and binning). Build and run with:
//   make unittests && ./unittests
// The tests only cover functions with header-visible definitions or external linkage; they use no
// input files and no MPI communication, and a non-zero exit code means at least one check failed.
// (Compile-time checks of the constexpr helpers live in static_asserts next to their definitions.)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <print>
#include <string_view>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "exspec.h"
#include "gammapkt.h"
#include "globals.h"
#include "input.h"
#include "macroatom.h"
#include "mpi_logging.h"
#include "radfield.h"
#include "random.h"
#include "rpkt.h"
#include "sn3d.h"
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
    for (double nu = nu_out_bound * (1. - 1e-13); nu < nu_out_bound; nu = std::nextafter(nu, nu_out_bound)) {
      const auto sigma = photoionisation_crosssection_fromtable(photoion_xs, nu_edge, nu);
      tail_reads_in_table = tail_reads_in_table && (std::ranges::find(photoion_xs, sigma) != photoion_xs.end());
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

  std::println("unit tests: {} of {} checks passed", checks_total - checks_failed, checks_total);

  return (checks_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
