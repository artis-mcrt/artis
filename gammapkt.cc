// Gamma-ray packet transport: radioactive decay line spectra, Compton scattering, photoelectric
// absorption, and pair production, propagating gamma packets through the ejecta until they
// deposit their energy or escape.

#include "gammapkt.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <limits>
#include <print>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "mpi_logging.h"
#include "packet.h"
#include "random.h"
#include "stats.h"
#include "vectors.h"

namespace gammapkt {
// Code for handing gamma rays - creation and propagation

namespace {
struct GammaLine {
  double energy{};  // in erg
  double probability{};
};

std::vector<std::vector<GammaLine>> gamma_spectra;

struct ElementPhotoionData {
  double energy;  // energy in MeV
  double sigma_xcom;  // cross section in barns/atom
};

// the XCOM table covers Z = 1..xcom_max_atomic_number (indexed by Z - 1)
constexpr int xcom_max_atomic_number = USE_XCOM_GAMMAPHOTOION ? 100 : 0;

std::array<std::vector<ElementPhotoionData>, xcom_max_atomic_number> photoion_data;

// Gamma-ray energies expressed as frequencies [Hz]. NB: these are very slightly inconsistent with MEV / H
// from constants.h (which gives 2.41805e+20 for 1 MeV); the historical values are kept here so that
// results are unchanged.
constexpr double nu_100kev = 2.41326e+19;
constexpr double nu_1mev = 2.41326e+20;
constexpr double nu_1p022mev = 2.46636e+20;  // electron-positron pair rest mass energy (pair production threshold)
constexpr double nu_1p5mev = 3.61990e+20;

void read_gamma_spectrum(const int nucindex, const std::string& filename) {
  // read the gamma ray lines and store the average energy in gamma rays per nuclear decay
  auto gammafile = fstream_required(filename, std::ios::in);
  std::string line;
  get_noncommentline(gammafile, line);
  std::istringstream ssline(line);
  int nlines = 0;
  ssline >> nlines;

  gamma_spectra[nucindex].reserve(nlines);
  gamma_spectra[nucindex].clear();

  double E_gamma_avg = 0.;
  for (int n = 0; n < nlines; n++) {
    get_noncommentline(gammafile, line);
    double en_mev = 0.;
    double prob = 0.;
    ssline.clear();
    ssline.str(line);
    assert_always(ssline >> en_mev >> prob);
    gamma_spectra[nucindex].push_back({.energy = en_mev * MEV, .probability = prob});
    E_gamma_avg += en_mev * MEV * prob;
  }

  decay::set_nucdecayenergygamma(nucindex, E_gamma_avg);
}

void set_trivial_gamma_spectrum(const int nucindex) {
  // there is no gamma-ray table, so just set a single gamma-ray line with 100% probability
  const int nlines = 1;
  gamma_spectra[nucindex].resize(nlines, {});
  gamma_spectra[nucindex][0].energy = decay::nucdecayenergygamma(nucindex);
  gamma_spectra[nucindex][0].probability = 1.;
}

void read_decaydata() {
  // migrate from old filenames that didn't specify the nuclide mass number
  if (!std::filesystem::exists("gamma_ni56.txt") && std::filesystem::exists("ni_lines.txt")) {
    std::error_code rename_error;
    std::filesystem::rename("ni_lines.txt", "gamma_ni56.txt", rename_error);
    if (rename_error) {
      printlnlog("[error] failed to move ni_lines.txt to gamma_ni56.txt: {}", rename_error.message());
    } else {
      printlnlog("Moved ni_lines.txt to gamma_ni56.txt");
    }
  }

  if (!std::filesystem::exists("gamma_co56.txt") && std::filesystem::exists("co_lines.txt")) {
    std::error_code rename_error;
    std::filesystem::rename("co_lines.txt", "gamma_co56.txt", rename_error);
    if (rename_error) {
      printlnlog("[error] failed to move co_lines.txt to gamma_co56.txt: {}", rename_error.message());
    } else {
      printlnlog("Moved co_lines.txt to gamma_co56.txt");
    }
  }

  gamma_spectra.resize(decay::get_num_nuclides(), {});
  int tables_found = 0;
  int nuclides_without_table = 0;
  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    gamma_spectra[nucindex].clear();
    const int z = decay::get_nuc_z(nucindex);
    const int a = decay::get_nuc_a(nucindex);
    if (z < 1) {
      continue;
    }

    auto strelname = decay::get_elname(z);
    std::ranges::transform(strelname, strelname.begin(), [](unsigned char c) { return std::tolower(c); });

    const auto striso = std::format("{}{}", strelname, a);

    // look in the current folder first, then in the data/ subfolder
    const std::array filenames = {std::format("gamma_{}.txt", striso), std::format("{}_lines.txt", striso)};

    bool tablefound = false;
    for (const auto& datadir : datafolders) {
      for (const auto& filename : filenames) {
        const auto filepath = std::format("{}{}", datadir, filename);
        if (std::filesystem::exists(filepath)) {
          tablefound = true;
          tables_found++;
          read_gamma_spectrum(nucindex, filepath);
          break;
        }
      }
      if (tablefound) {
        break;
      }
    }

    if (!tablefound && decay::nucdecayenergygamma(nucindex) > 0.) {
      assert_always(z != 28 || a != 56);  // Ni-56 must have a gamma spectrum
      assert_always(z != 27 || a != 56);  // Co-56 must have a gamma spectrum
      assert_always(z != 23 || a != 48);  // V-48 must have a gamma spectrum
      assert_always(z != 24 || a != 48);  // Cr-48 must have a gamma spectrum
      assert_always(z != 28 || a != 57);  // Ni-57 must have a gamma spectrum if present in list of nuclides
      assert_always(z != 27 || a != 57);  // Co-57 must have a gamma spectrum if present in list of nuclides
      set_trivial_gamma_spectrum(nucindex);
      nuclides_without_table++;
    }
  }

  if (decay::nuc_exists(26, 52)) {
    decay::set_nucdecayenergygamma(decay::get_nucindex(26, 52), 0.86 * MEV);  // Fe52
  }
  if (decay::nuc_exists(25, 52)) {
    decay::set_nucdecayenergygamma(decay::get_nucindex(25, 52), 3.415 * MEV);  // Mn52
  }

  printlnlog("[info] read gamma-ray line tables for {} nuclides", tables_found);
  if (nuclides_without_table > 0) {
    printlnlog(
        "[info] no gamma-ray line table found for {} nuclides with non-zero gamma decay energy, so a single line "
        "carrying the mean gamma energy per decay is used for each",
        nuclides_without_table);
  }
}

// construct an energy ordered gamma ray line list.
void init_gamma_linelist() {
  read_decaydata();

  // make a combined list of all gamma lines, sorted by energy, and write it to a file

  const ptrdiff_t total_lines = std::ranges::fold_left(
      gamma_spectra, 0, [](const ptrdiff_t sum, const auto& lines) { return sum + std::ssize(lines); });
  printlnlog("total gamma-ray lines {}", total_lines);

  if (globals::my_rank == 0) {
    struct NucGammaLine {
      int nucindex;  // is it a Ni56, Co56, a fake line, etc
      int nucgammaindex;  // which of the lines of that nuclide is it
      double energy;  // in erg
    };
    std::vector<NucGammaLine> allnuc_gamma_line_list;
    allnuc_gamma_line_list.reserve(total_lines);

    for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
      for (int j = 0; j < std::ssize(gamma_spectra[nucindex]); j++) {
        allnuc_gamma_line_list.push_back(
            {.nucindex = nucindex, .nucgammaindex = j, .energy = gamma_spectra[nucindex][j].energy});
      }
    }

    assert_always(std::ssize(allnuc_gamma_line_list) == total_lines);
    std::ranges::SORT_OR_STABLE_SORT(allnuc_gamma_line_list, [](const NucGammaLine& g1, const NucGammaLine& g2) {
      return std::tie(g1.energy, g1.nucindex, g1.nucgammaindex) < std::tie(g2.energy, g2.nucindex, g2.nucgammaindex);
    });

    auto gammalinelist = fstream_required("gammalinelist.out", std::ofstream::out | std::ofstream::trunc);
    std::println(gammalinelist, "#index nucindex Z A nucgammmaindex en_gamma_mev gammaline_probability");

    for (auto i = 0Z; i < total_lines; i++) {
      const int nucindex = allnuc_gamma_line_list[i].nucindex;
      const int index = allnuc_gamma_line_list[i].nucgammaindex;
      std::println(gammalinelist, "{} {} {} {} {} {:g} {:g}", static_cast<int>(i), allnuc_gamma_line_list[i].nucindex,
                   decay::get_nuc_z(allnuc_gamma_line_list[i].nucindex),
                   decay::get_nuc_a(allnuc_gamma_line_list[i].nucindex), allnuc_gamma_line_list[i].nucgammaindex,
                   gamma_spectra[nucindex][index].energy / MEV, gamma_spectra[nucindex][index].probability);
    }
    printlnlog("Wrote combined gamma-ray line list to gammalinelist.out");
  }
}

void init_xcom_photoion_data() {
  printlnlog("reading XCOM photoionisation data...");
  for (int Z = 0; Z < xcom_max_atomic_number; Z++) {
    photoion_data[Z].reserve(100);
  }

  auto data_fs = fstream_required("xcom_photoion_data.txt", std::ios::in);
  std::string line_str;
  while (get_noncommentline(data_fs, line_str)) {
    int Z = 0;
    double E = 0;
    double sigma = 0;
    std::stringstream(line_str) >> Z >> E >> sigma;
    assert_always(Z > 0);
    assert_always(Z <= xcom_max_atomic_number);
    // convert XCOM data to cgs units already here
    photoion_data[Z - 1].push_back({.energy = E, .sigma_xcom = sigma * 1e-24});
  }
}

// the absorption coefficient [cm^-1] for Compton scattering in the co-moving frame
[[nodiscard]] auto get_chi_compton_cmf(const int nonemptymgi, const double nu_cmf) -> double {
  if constexpr (GAMMA_USE_KAPPA_GREY.has_value()) {
    return 0.;
  }

  const double xx = H * nu_cmf / ME / CLIGHT / CLIGHT;

  // Use this to decide whether the Thompson limit is acceptable.

  const double sigma_cmf = (xx < THOMSON_LIMIT) ? SIGMA_T : sigma_compton_partial(xx, 1 + (2 * xx));

  // Now need to multiply by the electron number density.
  const double chi_cmf = sigma_cmf * grid::get_nnetot(nonemptymgi);

  assert_testmodeonly(std::isfinite(chi_cmf));

  return chi_cmf;
}

// For Thomson scattering we can get the new angle from a random number very easily.
auto thomson_angle(rngstate_type& rngstate) -> double {
  const double B_coeff = (8. * rng_uniform(rngstate)) - 4.;

  const double t_coeff = std::cbrt((std::sqrt(pow2(B_coeff) + 4) - B_coeff) / 2);

  const double mu = (1 / t_coeff) - t_coeff;

  assert_testmodeonly(fabs(mu) <= 1);

  return mu;
}

// scattering a direction through angle theta.
[[nodiscard]] auto scatter_dir(const Vec3d& dir_in, const double cos_theta, rngstate_type& rngstate) -> Vec3d {
  // begin with setting the direction in coordinates where original direction is parallel to z-hat.

  assert_testmodeonly(std::fabs(vec_len(dir_in) - 1.) < 1e-10);  // dir_in must be a unit vector

  const double phi = rng_uniform(rngstate) * 2 * PI;

  const double sin_theta_sq = 1. - pow2(cos_theta);
  const double sin_theta = std::sqrt(sin_theta_sq);
  const double zprime = cos_theta;
  const double xprime = sin_theta * cos(phi);
  const double yprime = sin_theta * sin(phi);

  // When dir_in is (anti)parallel to the z-axis the rotation below is singular (norm1 -> inf, giving
  // 0*inf = NaN). Handle it directly: the scattering frame's z-axis is dir_in, so just (anti)align the
  // result along z (matching the pole handling in electron_scatter_rpkt).
  if (std::fabs(dir_in[2]) > 0.999999999) {
    const auto dir_out = Vec3d{xprime, yprime, (dir_in[2] > 0) ? zprime : -zprime};
    assert_testmodeonly(std::fabs(vec_len(dir_out) - 1.) < 1e-10);
    return dir_out;
  }

  // Now need to derotate the coordinates back to real x,y,z. Rotation matrix is determined by dir_in.

  const double norm1 = 1. / std::sqrt(pow2(dir_in[0]) + pow2(dir_in[1]));
  const double norm2 = 1. / vec_len(dir_in);

  const double r11 = dir_in[1] * norm1;
  const double r12 = -dir_in[0] * norm1;
  const double r13 = 0.;
  const double r21 = dir_in[0] * dir_in[2] * norm1 * norm2;
  const double r22 = dir_in[1] * dir_in[2] * norm1 * norm2;
  const double r23 = -norm2 / norm1;
  const double r31 = dir_in[0] * norm2;
  const double r32 = dir_in[1] * norm2;
  const double r33 = dir_in[2] * norm2;

  const auto dir_out = Vec3d{
      (r11 * xprime) + (r21 * yprime) + (r31 * zprime),
      (r12 * xprime) + (r22 * yprime) + (r32 * zprime),
      (r13 * xprime) + (r23 * yprime) + (r33 * zprime),
  };

  assert_testmodeonly(std::fabs(vec_len(dir_out) - 1.) < 1e-10);

  return dir_out;
}

// handle physical Compton scattering event
void compton_scatter(Packet& pkt) {
  const double xx = H * pkt.nu_cmf / ME / CLIGHT / CLIGHT;

  // It is known that a Compton scattering event is going to take place.
  // We need to do two things - (1) decide whether to convert energy
  // to electron or leave as gamma (2) decide properties of new packet.

  // The probability of giving energy to electron is related to the
  // energy change of the gamma ray. This is equivalent to the choice of
  // scattering angle. Probability of scattering into particular angle
  // (i.e. final energy) is related to the partial cross-section.

  // Choose a random number to get the energy. Want to find the factor by which the energy changes "f" such that
  // sigma_partial/sigma_tot = zrand

  // initialise with Thomson limit case (no energy loss)
  double f = 1.;
  bool stay_gamma = true;
  if (xx >= THOMSON_LIMIT) {
    f = choose_f(xx, rng_uniform(get_rngstate(pkt)));

    // Prob of keeping gamma ray is...
    const double prob_gamma = 1. / f;

    stay_gamma = (rng_uniform(get_rngstate(pkt)) < prob_gamma);
  }

  if (stay_gamma) {
    // It stays as a gamma ray. Change frequency and direction in co-moving frame then transfer back to rest frame.

    pkt.nu_cmf = pkt.nu_cmf / f;  // reduce frequency

    // The packet has stored the direction in the rest frame.
    // Use aberration of angles to get this into the co-moving frame.

    const auto vel_vec = get_velocity(pkt.pos, pkt.prop_time);

    const auto cmf_dir = angle_ab(pkt.dir, vel_vec);

    // Now change the direction through the scattering angle.

    const double cos_theta = (xx < THOMSON_LIMIT) ? thomson_angle(get_rngstate(pkt)) : 1. - ((f - 1) / xx);

    const auto new_dir = scatter_dir(cmf_dir, cos_theta, get_rngstate(pkt));

    assert_testmodeonly(fabs(1. - dot(new_dir, new_dir)) < 1e-8);
    assert_testmodeonly(fabs(dot(new_dir, cmf_dir) - cos_theta) < 1e-8);

    // Now convert back again.

    pkt.dir = angle_ab(new_dir, vec_scale(vel_vec, -1.));

    assert_testmodeonly(std::fabs(vec_len(pkt.dir) - 1.) < 1e-10);

    // It now has a rest frame direction and a co-moving frequency. Just need to set the rest frame energy.
    set_pkt_restframe_from_cmf(pkt);
  } else {
    // energy loss of the gamma becomes energy of the electron (needed to calculate time-dependent thermalisation rate)
    if constexpr (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::TIMEDEPENDENTWITHGAMMAPRODUCTS) {
      pkt.nu_cmf = pkt.nu_cmf * (1 - (1 / f));
      pkt.type = TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS;
    } else {
      pkt.type = TYPE_NTLEPTON_DEPOSITED;
    }
    pkt.absorptiontype = ABSTYPE_GAMMA_COMPTON;
    stats::increment(stats::Counter::NT_STAT_FROM_GAMMA);
  }
}

// calculate the absorption coefficient [cm^-1] for photo electric effect scattering in the co-moving frame
[[nodiscard]] auto get_chi_photo_electric_cmf(const int nonemptymgi, const double ffegrp, const double nu_cmf)
    -> double {
  const double rho = grid::get_rho(nonemptymgi);

  if constexpr (GAMMA_USE_KAPPA_GREY.has_value()) {
    return GAMMA_USE_KAPPA_GREY.value() * rho;
  }

  if constexpr (!USE_XCOM_GAMMAPHOTOION) {
    // Cross sections from Equation 2 of Ambwani & Sutherland (1988), attributed to Veigele (1973)

    const double hnu_over_100kev = nu_cmf / nu_100kev;

    const double sigma_cmf_si = 1.16e-24 * pow(hnu_over_100kev, -3.13);

    const double sigma_cmf_fe = 25.7e-24 * pow(hnu_over_100kev, -3.0);

    // Now need to multiply by the particle number density. The composition is approximated as a mix of just
    // two species: an iron-group one with mass fraction ffegrp, and silicon for the remainder. The nuclide
    // number densities therefore use A = 56 for the iron group and A = 28 for silicon.

    const double chi_cmf_si = sigma_cmf_si * (rho / MH / 28);

    const double chi_cmf_fe = sigma_cmf_fe * (rho / MH / 56);

    return (chi_cmf_fe * ffegrp) + (chi_cmf_si * (1. - ffegrp));
  }

  const double hnu_over_1MeV = nu_cmf / nu_1mev;
  const double log10_hnu_over_1MeV = log10(hnu_over_1MeV);
  double chi_cmf{0.};
  for (int i = 0; i < get_nelements(); i++) {
    // determine charge number:
    const int Z = get_atomicnumber(i);
    if (Z > xcom_max_atomic_number) {
      continue;  // no XCOM photoionisation data available (table only covers Z = 1..xcom_max_atomic_number)
    }
    const auto numb_energies = std::ssize(photoion_data[Z - 1]);
    if (numb_energies == 0) {
      continue;
    }
    const double n_i = grid::get_elem_numberdens(nonemptymgi, i);  // number density in the current cell
    if (n_i == 0) {
      continue;
    }
    // get indices of lower and upper boundary
    int idx_above = -1;

    for (int j = 0; j < numb_energies; j++) {
      if (photoion_data[Z - 1][j].energy > hnu_over_1MeV) {
        idx_above = j;
        break;
      }
    }
    if (idx_above == 0) {  // packet energy smaller than all tabulated values
      chi_cmf += photoion_data[Z - 1][0].sigma_xcom * n_i;
      continue;
    }
    if (idx_above == -1) {  // packet energy greater than all tabulated values
      chi_cmf += photoion_data[Z - 1][numb_energies - 1].sigma_xcom * n_i;
      continue;
    }
    assert_always(idx_above > 0);
    assert_always(idx_above < numb_energies);
    const int idx_below = idx_above - 1;
    assert_always(idx_below >= 0);
    const double log10_E = log10_hnu_over_1MeV;
    const double log10_E_above = log10(photoion_data[Z - 1][idx_above].energy);
    const double log10_E_below = log10(photoion_data[Z - 1][idx_below].energy);
    const double log10_sigma_below = log10(photoion_data[Z - 1][idx_below].sigma_xcom);
    const double log10_sigma_above = log10(photoion_data[Z - 1][idx_above].sigma_xcom);
    // interpolate or extrapolate, both linear in log10-log10 space
    const double log10_sigma_interp = log10_sigma_below + ((log10_sigma_above - log10_sigma_below) /
                                                           (log10_E_above - log10_E_below) * (log10_E - log10_E_below));
    const double sigma_interp = pow(10., log10_sigma_interp);
    const double chi_cmf_contrib = sigma_interp * n_i;
    assert_always(sigma_interp >= 0.);
    chi_cmf += chi_cmf_contrib;
  }
  return chi_cmf;
}

// energy-dependent factor of the pair-production cross section, from Equation 2 of Ambwani &
// Sutherland (1988), attributed to Hubbell (1969). Multiply by Z^2 * 1e-27 for a cross section in
// cm^2. Only valid above the 1.022 MeV threshold.
[[nodiscard]] constexpr auto get_sigma_pair_prod_factor(const double nu_cmf) -> double {
  const double hnu_over_mev = nu_cmf / nu_1mev;
  if (nu_cmf > nu_1p5mev) {
    return 0.0481 + (0.301 * (hnu_over_mev - 1.5));
  }
  // the coefficient below 1.5 MeV is 0.10063, not the 1.0063 printed in the paper: the latter makes
  // the fit jump by a factor of ten at the junction between the two branches
  return 0.10063 * (hnu_over_mev - 1.022);
}

// the two branches of the fit must agree where they meet at 1.5 MeV
static_assert((get_sigma_pair_prod_factor(nu_1p5mev * (1. + 1e-12)) - get_sigma_pair_prod_factor(nu_1p5mev)) < 1e-6);
static_assert((get_sigma_pair_prod_factor(nu_1p5mev) - get_sigma_pair_prod_factor(nu_1p5mev * (1. + 1e-12))) < 1e-6);

// calculate the absorption coefficient [cm^-1] for pair production in the comoving frame
[[nodiscard]] auto get_chi_pair_prod_cmf(const int nonemptymgi, const double ffegrp, const double nu_cmf) -> double {
  if constexpr (GAMMA_USE_KAPPA_GREY.has_value()) {
    return 0.;
  }
  const double rho = grid::get_rho(nonemptymgi);

  // below the pair production threshold there is no pair production
  if (nu_cmf <= nu_1p022mev) {
    return 0.;
  }

  const double sigma_factor = get_sigma_pair_prod_factor(nu_cmf);

  const double sigma_cmf_si = sigma_factor * 196.e-27;

  const double sigma_cmf_fe = sigma_factor * 784.e-27;

  // multiply by the particle number density. As in get_chi_photo_electric_cmf(), the composition is
  // approximated as an iron-group species (A = 56) with mass fraction ffegrp plus silicon (A = 28).

  const double chi_cmf_si = sigma_cmf_si * (rho / MH / 28);

  const double chi_cmf_fe = sigma_cmf_fe * (rho / MH / 56);

  const double chi_cmf = (chi_cmf_fe * ffegrp) + (chi_cmf_si * (1. - ffegrp));

  return std::max(chi_cmf, 0.);
}

// get the comoving-frame gamma-ray absorption coefficient (with the expected energy loss fraction per interaction
// factor included). All three terms are comoving-frame coefficients so that they can be combined and converted to
// the deposition estimator with a single frame factor in update_gamma_dep().
[[nodiscard]] auto get_chi_cmf_loss_weighted(const int nonemptymgi, const double nu_cmf) -> double {
  assert_testmodeonly(nonemptymgi >= 0);
  const double ffegrp = grid::get_ffegrp(grid::get_mgi_of_nonemptymgi(nonemptymgi));
  const auto chi_photo_electric_cmf = get_chi_photo_electric_cmf(nonemptymgi, ffegrp, nu_cmf);

  if constexpr (GAMMA_USE_KAPPA_GREY.has_value()) {
    // grey-mode transport interacts only via the grey photoabsorption opacity (the Compton and
    // pair-production coefficients are zero), and every interaction deposits the full packet
    // energy, so the estimator must not include the Compton energy-loss term
    return chi_photo_electric_cmf;
  }

  const double xx = H * nu_cmf / ME / CLIGHT / CLIGHT;
  const auto chi_pair_prod_cmf = get_chi_pair_prod_cmf(nonemptymgi, ffegrp, nu_cmf);

  return ((meanf_sigma(xx) * grid::get_nnetot(nonemptymgi)) + chi_photo_electric_cmf +
          (chi_pair_prod_cmf * (1. - (nu_1p022mev / nu_cmf))));
}

// update the energy deposition estimator for gamma ray path increment
void update_gamma_dep(const Packet& pkt, const int nonemptymgi, const double dist) {
  if (!(dist > 0)) {
    return;
  }
  if constexpr (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::TIMEDEPENDENTWITHGAMMAPRODUCTS) {
    return;  // don't instantly deposit energy from gamma rays, handle the particles they produce instead
  }

  if (nonemptymgi < 0) {
    return;  // empty cell
  }

  // Comoving-frame energy deposited along a rest-frame path increment dist is
  //   chi_cmf * e_cmf * doppler * dist = chi_cmf * e_rf * doppler^2 * dist,
  // since e_cmf = e_rf * doppler and the interaction probability along the rest-frame path is
  // chi_rf * dist = chi_cmf * doppler * dist.
  const double doppler_sq = pow2(calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time));
  const double heating_cont = get_chi_cmf_loss_weighted(nonemptymgi, pkt.nu_cmf) * pkt.e_rf * dist * doppler_sq;

  // The terms in the above are for Compton, photoelectric and pair production. The pair production one
  // assumes that a fraction (1. - (1.022 MeV / nu)) of the gamma's energy is thermalised.
  // The remaining 1.022 MeV is made into gamma rays

  // For normalisation this needs to be
  //  1) divided by volume
  //  2) divided by the length of the time step
  //  3) divided by 4 pi sr
  //  This will all be done later
  assert_testmodeonly(heating_cont >= 0.);
  assert_testmodeonly(std::isfinite(heating_cont));
  atomicadd(globals::dep_estimator_gamma[nonemptymgi], heating_cont);
}

// Give a gamma packet an isotropic direction in the comoving frame, then set its rest-frame
// direction, frequency and energy. Assumes pkt.nu_cmf and pkt.e_cmf are already set.
DEVICE_FUNC void emit_gamma_isotropic(Packet& pkt) {
  const auto dir_cmf = get_rand_isotropic_unitvec(get_rngstate(pkt));

  // This direction is in the cmf - convert it to the rest frame using aberration of angles.
  // Converting cmf -> rest requires the -ve velocity (negative time for the backwards transform).
  const auto vel_vec = get_velocity(pkt.pos, -pkt.prop_time);

  pkt.dir = angle_ab(dir_cmf, vel_vec);

  set_pkt_restframe_from_cmf(pkt);

  pkt.type = TYPE_GAMMA;
}

// handle gamma to electron-positron pair production event
void pair_production(Packet& pkt) {
  // In pair production, the original gamma makes an electron-positron pair with total kinetic energy equal to the
  // gamma-ray energy minus the 1.022 MeV rest-mass energy. The indivisible-packet branch selection assigns either this
  // kinetic-energy share or the rest-mass energy to the packet. The latter produces 0.511 MeV annihilation gamma rays.

  constexpr double pair_rest_mass_energy = 1.022 * MEV;
  const double gamma_energy = H * pkt.nu_cmf;
  const double prob_gamma = pair_rest_mass_energy / gamma_energy;

  assert_always(prob_gamma >= 0);

  if (rng_uniform(get_rngstate(pkt)) > prob_gamma) {
    if constexpr (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::TIMEDEPENDENTWITHGAMMAPRODUCTS) {
      // Represent one of the two particles, assuming that they share the available kinetic energy equally. e_cmf stays
      // unchanged because it is the indivisible macro-packet energy; nu_cmf stores the physical particle energy used by
      // the time-dependent thermalisation calculation.
      const double particle_kinetic_energy = (gamma_energy - pair_rest_mass_energy) / 2;
      assert_testmodeonly(particle_kinetic_energy > 0.);
      pkt.nu_cmf = particle_kinetic_energy / H;
      pkt.type = (rng_uniform(get_rngstate(pkt)) > 0.5) ? TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS
                                                        : TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS;
    } else {
      pkt.type = TYPE_NTLEPTON_DEPOSITED;
    }

    pkt.absorptiontype = ABSTYPE_GAMMA_PAIRPRODUCTION;
    stats::increment(stats::Counter::NT_STAT_FROM_GAMMA);
  } else {
    // The energy goes into emission at 511 keV.
    pkt.nu_cmf = 0.511 * MEV / H;

    // Now give the gamma ray an isotropic direction and set the rest-frame quantities.
    emit_gamma_isotropic(pkt);
  }
}

// move a gamma packet until time t2
void transport_gamma(Packet& pkt, const double t2) {
  // Assign optical depth to next physical event. And start counter of optical depth for this path.
  const double tau_next = -std::log(static_cast<double>(rng_uniform_pos(get_rngstate(pkt))));

  // Start by finding the distance to the crossing of the grid cell
  // boundaries. boundarydist is the boundary distance and next_cellindex is the grid cell into which we pass.

  const auto [boundarydist, next_cellindex] = grid::boundary_distance(pkt.dir, pkt.pos, pkt.prop_time, pkt.cellindex);

  // Now consider the scattering/destruction processes.
  // Compton scattering - need to determine the scattering co-efficient.
  // Routine returns the value in the rest frame.
  const int mgi = grid::get_propcell_modelgridindex(pkt.cellindex);
  const auto nonemptymgi = (mgi >= 0) ? grid::get_nonemptymgi_of_mgi(mgi) : -1;

  const auto doppler = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
  const double ffegrp = (mgi >= 0) ? grid::get_ffegrp(mgi) : 0.;
  const double chi_compton = (mgi >= 0) ? get_chi_compton_cmf(nonemptymgi, pkt.nu_cmf) * doppler : 0.;
  const double chi_photo_electric =
      (mgi >= 0) ? get_chi_photo_electric_cmf(nonemptymgi, ffegrp, pkt.nu_cmf) * doppler : 0.;
  const double chi_pair_prod = (mgi >= 0) ? get_chi_pair_prod_cmf(nonemptymgi, ffegrp, pkt.nu_cmf) * doppler : 0.;
  const double chi_tot = chi_compton + chi_photo_electric + chi_pair_prod;

  assert_testmodeonly(std::isfinite(chi_compton));
  assert_testmodeonly(std::isfinite(chi_photo_electric));
  assert_testmodeonly(std::isfinite(chi_pair_prod));

  // So distance before physical event is...

  const double edist = chi_tot > 0. ? tau_next / chi_tot : std::numeric_limits<double>::max();

  assert_always(edist >= 0);

  // Find how far it can travel during the time interval.

  const double tdist = (t2 - pkt.prop_time) * CLIGHT_PROP;

  assert_always(tdist >= 0);

  // boundarydist can be zero when rounding error has put the packet just past a cell
  // boundary that it is moving towards (an immediate crossing). This cannot repeat
  // endlessly: the just-crossed face is excluded by the velocity conditions in the new
  // cell, so at most ndim consecutive zero-distance crossings can occur (e.g. at a corner)
  if ((boundarydist <= tdist) && (boundarydist <= edist)) {
    move_pkt_withtime(pkt, boundarydist / 2.);

    // Move it into the new cell.
    if (chi_tot > 0) {
      update_gamma_dep(pkt, nonemptymgi, boundarydist);
    }

    move_pkt_withtime(pkt, boundarydist / 2.);

    if (next_cellindex != pkt.cellindex) {
      grid::change_cell_or_escape(pkt, next_cellindex);
    }
  } else if ((tdist < boundarydist) && (tdist <= edist)) {
    // Doesn't reach boundary.
    move_pkt_withtime(pkt, tdist / 2.);

    if (chi_tot > 0) {
      update_gamma_dep(pkt, nonemptymgi, tdist);
    }
    move_pkt_withtime(pkt, tdist / 2.);
    pkt.prop_time = t2;  // prevent roundoff error
  } else if ((edist < boundarydist) && (edist <= tdist)) {
    move_pkt_withtime(pkt, edist / 2.);
    if (chi_tot > 0) {
      update_gamma_dep(pkt, nonemptymgi, edist);
    }
    move_pkt_withtime(pkt, edist / 2.);

    // event occurs. Choose which event and call the appropriate subroutine.
    const double chi_rnd = rng_uniform(get_rngstate(pkt)) * chi_tot;
    if (chi_compton > chi_rnd) {
      // Compton scattering.
      compton_scatter(pkt);
    } else if ((chi_compton + chi_photo_electric) > chi_rnd) {
      // Photo electric effect
      if constexpr (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::TIMEDEPENDENTWITHGAMMAPRODUCTS) {
        pkt.type = TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS;
        // nu_cmf stays the same as the gamma-ray energy becomes the kinetic energy of the electron (minus ionisation
        // energy but this is neglected here)
      } else {
        pkt.type = TYPE_NTLEPTON_DEPOSITED;
      }

      pkt.absorptiontype = ABSTYPE_GAMMA_PHOTOELECTRIC;
      stats::increment(stats::Counter::NT_STAT_FROM_GAMMA);
    } else {
      // It's a pair production
      pair_production(pkt);
    }
  } else {
    assert_always(false);
  }
}

// With probability f_gamma, absorb the packet locally (depositing its energy as a k-packet);
// otherwise let it escape. Shared tail of the gamma-ray thermalisation schemes.
void absorb_or_escape_gamma(Packet& pkt, const double f_gamma) {
  assert_always(f_gamma >= 0.);
  assert_always(f_gamma <= 1.);

  if (rng_uniform(get_rngstate(pkt)) < f_gamma) {
    // packet is absorbed and contributes to the heating as a k-packet.
    // These parameterised thermalisation schemes don't resolve which process absorbed the gamma ray, so the
    // absorption is recorded under the photoelectric type (the dominant absorber at the relevant energies).
    pkt.type = TYPE_NTLEPTON_DEPOSITED;
    pkt.absorptiontype = ABSTYPE_GAMMA_PHOTOELECTRIC;
  } else {
    // let packet escape, i.e. make it inactive. change_cell_or_escape() records pkt.escape_type = pkt.type
    // before setting the type to TYPE_ESCAPE, so leave pkt.type as TYPE_GAMMA here so the escaped gamma is
    // correctly identified (e.g. for the gamma light curve).
    grid::change_cell_or_escape(pkt, -99);
  }
}

void barnes_thermalisation(Packet& pkt)
// Barnes treatment: packet is either getting absorbed immediately and locally
// creating a k-packet or it escapes. The absorption probability matches the
// Barnes thermalization efficiency, for expressions see the original paper:
// https://ui.adsabs.harvard.edu/abs/2016ApJ...829..110B
{
  // compute thermalization efficiency (= absorption probability) using a mean gamma-ray opacity
  // of 0.1, an average value chosen to fit the analytic approximations from the paper.
  // Alternative: distinguish between low-E (kappa = 1) and high-E (kappa = 0.05) packets.

  // determine average initial density via kinetic energy
  const double E_kin = grid::get_ejecta_kinetic_energy();
  const double v_ej = sqrt(E_kin * 2 / grid::mtot_input);

  // t_ineff = sqrt(rho_0 * R_0 * t_0^2 * mean_gamma_opac), expressed via the paper's scaling relation.
  const double t_ineff = 1.4 * DAY * sqrt(grid::mtot_input / (5.e-3 * MSUN)) * ((0.2 * CLIGHT) / v_ej);
  const double tau = pow2(t_ineff / pkt.prop_time);
  const double f_gamma = 1. - exp(-tau);

  absorb_or_escape_gamma(pkt, f_gamma);
}

void wollaeger_thermalisation(Packet& pkt) {
  // corresponds to a local version of the Barnes scheme, i.e. it takes into account the local mass
  // density rather than a value averaged over the ejecta
  constexpr double mean_gamma_opac = 0.1;
  // integration: requires distances within single cells in radial direction and the corresponding densities
  // need to create a packet copy which is moved during the integration
  Packet pkt_copy = pkt;
  pkt_copy.dir = vec_norm(pkt_copy.pos);  // integrate the optical depth radially outwards
  const double t_current = pkt.prop_time;
  double tau = 0.;
  bool end_packet = false;
  while (!end_packet) {
    // distance to the next cell
    const auto [boundarydist, next_cellindex] =
        grid::boundary_distance(pkt_copy.dir, pkt_copy.pos, pkt_copy.prop_time, pkt_copy.cellindex);
    const double s_cont = boundarydist * pow3(t_current / pkt_copy.prop_time);
    const int mgi = grid::get_propcell_modelgridindex(pkt_copy.cellindex);
    if (mgi >= 0) {
      const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
      tau += grid::get_rho(nonemptymgi) * s_cont * mean_gamma_opac;  // contribution to the integral
    }
    // move packet copy now
    move_pkt_withtime(pkt_copy, boundarydist);

    grid::change_cell_or_escape(pkt_copy, next_cellindex, false);
    end_packet = (pkt_copy.type == TYPE_ESCAPE);
  }
  const double f_gamma = 1. - std::exp(-tau);

  absorb_or_escape_gamma(pkt, f_gamma);
}

void guttman_thermalisation(Packet& pkt) {
  // Guttman et al. (2024), doi:10.1093/mnras/stae1795.
  // Extension of the Wollaeger scheme that averages the deposition probability over random emission directions.

  // Mean gamma opacity from section 3.2, using the lower value for nearly symmetric matter at late times.
  constexpr double mean_gamma_opac = 0.03;

  constexpr int num_directions = 100;
  double deposition_probability_sum = 0.;
  for (int i = 0; i < num_directions; i++) {
    Packet pkt_copy = pkt;
    pkt_copy.dir = get_rand_isotropic_unitvec(get_rngstate(pkt));

    double tau = 0.;
    while (pkt_copy.type != TYPE_ESCAPE) {
      const auto [boundarydist, next_cellindex] =
          grid::boundary_distance(pkt_copy.dir, pkt_copy.pos, pkt_copy.prop_time, pkt_copy.cellindex);
      const int mgi = grid::get_propcell_modelgridindex(pkt_copy.cellindex);
      if (mgi >= 0) {
        const double rho = grid::get_rho_tmin(mgi) * pow3(globals::tmin / pkt_copy.prop_time);
        tau += mean_gamma_opac * rho * boundarydist;
      }
      move_pkt_withtime(pkt_copy, boundarydist);
      grid::change_cell_or_escape(pkt_copy, next_cellindex, false);
    }

    deposition_probability_sum -= std::expm1(-tau);
  }

  absorb_or_escape_gamma(pkt, deposition_probability_sum / num_directions);
}

}  // anonymous namespace

void init_gamma_data() {
  init_gamma_linelist();
  if constexpr (USE_XCOM_GAMMAPHOTOION) {
    init_xcom_photoion_data();
  }
}

[[nodiscard]] auto choose_gamma_ray(const int nucindex, rngstate_type& rngstate) -> double {
  // Get the frequency [Hz] of a random gamma ray from the decay spectrum of a given nucleus.
  // If no gamma spectrum is known, return -1.

  if (gamma_spectra[nucindex].empty()) {
    // for historical consistency, Fe52 and Mn52 decay directly to k-packets, which is signalled by a negative nu_cmf
    const auto nuc_z = decay::get_nuc_z(nucindex);
    const auto nuc_a = decay::get_nuc_a(nucindex);
    assert_always((nuc_z == 26 && nuc_a == 52) || (nuc_z == 25 && nuc_a == 52));
    return -1;
  }

  const double E_gamma = decay::nucdecayenergygamma(nucindex);  // Average energy per gamma line of a decay

  const double zrand = rng_uniform(rngstate);
  double runtot = 0.;
  for (auto n = 0Z; n < std::ssize(gamma_spectra[nucindex]); n++) {
    runtot += gamma_spectra[nucindex][n].probability * gamma_spectra[nucindex][n].energy / E_gamma;
    // Strict comparison ensures that zrand == 0 does not select a leading zero-probability line.
    if (zrand < runtot) {
      return gamma_spectra[nucindex][n].energy / H;
    }
  }

  // the probabilities sum to one only up to floating-point rounding, so attribute the remainder to the last line
  return gamma_spectra[nucindex].back().energy / H;
}

// convert a pellet to a gamma ray (or kpkt if no gamma spec loaded)
DEVICE_FUNC void pellet_gamma_decay(Packet& pkt) {
  // Start by getting the position of the pellet at the point of decay. Pellet is moving with the matter.

  // if no gamma spectra is known, then convert straight to kpkts (e.g., Fe52, Mn52)
  if (pkt.nu_cmf < 0) {
    pkt.type = TYPE_KPKT;
    pkt.absorptiontype = ABSTYPE_PELLET_NOGAMMASPEC;
    return;
  }

  assert_testmodeonly(pkt.prop_time == pkt.tdecay);

  // Give the gamma ray an isotropic direction (isotropic emission in the cmf) and set the
  // rest-frame energy and frequency, recording that it's now a gamma ray.
  emit_gamma_isotropic(pkt);
}

DEVICE_FUNC void do_gamma(Packet& pkt, const int nts, const double t2) {
  if constexpr (GAMMA_THERMALISATION_SCHEME == GammaThermalisationScheme::FREQUENCYDEPENDENT) {
    transport_gamma(pkt, t2);
  } else if constexpr (GAMMA_THERMALISATION_SCHEME == GammaThermalisationScheme::BARNES) {
    barnes_thermalisation(pkt);
  } else if constexpr (GAMMA_THERMALISATION_SCHEME == GammaThermalisationScheme::WOLLAEGER) {
    wollaeger_thermalisation(pkt);
  } else if constexpr (GAMMA_THERMALISATION_SCHEME == GammaThermalisationScheme::GUTTMAN) {
    guttman_thermalisation(pkt);
  } else {
    assert_always(false);  // invalid thermalisation scheme
  }

  if (pkt.type != TYPE_GAMMA && pkt.type != TYPE_ESCAPE) {
    if constexpr (PARTICLE_THERMALISATION_SCHEME != ParticleThermalisationScheme::TIMEDEPENDENTWITHGAMMAPRODUCTS) {
      atomicadd(globals::timesteps[nts].gamma_dep_discrete, pkt.e_cmf);
    }

    if constexpr (GAMMA_THERMALISATION_SCHEME != GammaThermalisationScheme::FREQUENCYDEPENDENT) {
      // no transport, so the path-based gamma deposition estimator won't get updated unless we do it here
      const int mgi = grid::get_propcell_modelgridindex(pkt.cellindex);
      if (mgi >= 0) {  // empty cells have no deposition estimator
        const int nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
        atomicadd(globals::dep_estimator_gamma[nonemptymgi], pkt.e_cmf);
      }
    }
  }
}

}  // namespace gammapkt
