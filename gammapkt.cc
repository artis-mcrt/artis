#include "gammapkt.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <span>
#include <string>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "packet.h"
#include "sn3d.h"
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

struct el_photoion_data {
  double energy;      // energy in MeV
  double sigma_xcom;  // cross section in barns/atom
};

constexpr int numb_xcom_elements = USE_XCOM_GAMMAPHOTOION ? 100 : 0;

std::array<std::vector<el_photoion_data>, numb_xcom_elements> photoion_data;

struct NucGammaLine {
  int nucindex;       // is it a Ni56, Co56, a fake line, etc
  int nucgammaindex;  // which of the lines of that nuclide is it
  double energy;      // in erg
};

std::vector<NucGammaLine> allnuc_gamma_line_list;

void read_gamma_spectrum(const int nucindex, const char filename[50])
// reads in gamma_spectra and returns the average energy in gamma rays per nuclear decay
{
  printout("reading gamma spectrum for Z=%d A=%d from %s...", decay::get_nuc_z(nucindex), decay::get_nuc_a(nucindex),
           filename);

  FILE *filein = fopen_required(filename, "r");
  int nlines = 0;
  assert_always(fscanf(filein, "%d", &nlines) == 1);

  gamma_spectra[nucindex].resize(nlines, {});

  double E_gamma_avg = 0.;
  for (int n = 0; n < nlines; n++) {
    double en_mev = 0.;
    double prob = 0.;
    assert_always(fscanf(filein, "%lg %lg", &en_mev, &prob) == 2);
    gamma_spectra[nucindex][n].energy = en_mev * MEV;
    gamma_spectra[nucindex][n].probability = prob;
    E_gamma_avg += en_mev * MEV * prob;
  }
  fclose(filein);

  decay::set_nucdecayenergygamma(nucindex, E_gamma_avg);

  printout("nlines %d avg_en_gamma %g MeV\n", nlines, E_gamma_avg / MEV);
}

void set_trivial_gamma_spectrum(const int nucindex) {
  // printout("Setting trivial gamma spectrum for z %d a %d engamma %g\n", z, a, decay::nucdecayenergygamma(z, a));
  const int nlines = 1;
  gamma_spectra[nucindex].resize(nlines, {});
  gamma_spectra[nucindex][0].energy = decay::nucdecayenergygamma(nucindex);
  gamma_spectra[nucindex][0].probability = 1.;
}

void read_decaydata() {
  // migrate from old filename
  if (!std::filesystem::exists("ni56_lines.txt") && std::filesystem::exists("ni_lines.txt")) {
    printout("Moving ni_lines.txt to ni56_lines.txt\n");
    std::rename("ni_lines.txt", "ni56_lines.txt");
  }

  // migrate from old filename
  if (!std::filesystem::exists("co56_lines.txt") && std::filesystem::exists("co_lines.txt")) {
    printout("Moving co_lines.txt to co56_lines.txt\n");
    std::rename("co_lines.txt", "co56_lines.txt");
  }

  gamma_spectra.resize(decay::get_num_nuclides(), {});

  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    gamma_spectra[nucindex].clear();
    const int z = decay::get_nuc_z(nucindex);
    const int a = decay::get_nuc_a(nucindex);
    if (z < 1) {
      continue;
    }

    auto strelname = decay::get_elname(z);
    std::ranges::transform(strelname, strelname.begin(), [](unsigned char c) { return std::tolower(c); });

    // look in the current folder
    char filename[MAXFILENAMELENGTH];
    snprintf(filename, MAXFILENAMELENGTH, "%s%d_lines.txt", strelname.c_str(), a);

    // look in the 'data' subfolder
    char filename2[MAXFILENAMELENGTH];
    snprintf(filename2, MAXFILENAMELENGTH, "data/%s%d_lines.txt", strelname.c_str(), a);

    if (std::ifstream(filename)) {
      read_gamma_spectrum(nucindex, filename);
    } else if (std::ifstream(filename2)) {
      read_gamma_spectrum(nucindex, filename2);
    } else if (decay::nucdecayenergygamma(nucindex) > 0.) {
      // printout("%s does not exist. Setting 100%% chance of single gamma-line with energy %g MeV\n",
      //   filename, decay::nucdecayenergygamma(z, a) / EV / 1e6);
      set_trivial_gamma_spectrum(nucindex);

      assert_always(z != 28 || a != 56);  // Ni-56 must have a gamma spectrum
      assert_always(z != 27 || a != 56);  // Co-56 must have a gamma spectrum
      assert_always(z != 23 || a != 48);  // V-48 must have a gamma spectrum
      assert_always(z != 24 || a != 48);  // Cr-48 must have a gamma spectrum
      assert_always(z != 28 || a != 57);  // Ni-57 must have a gamma spectrum if present in list of nuclides
      assert_always(z != 28 || a != 57);  // Co-57 must have a gamma spectrum if present in list of nuclides
    } else {
      // printout("%s does not exist. No gamma decay from this nuclide.\n", filename);
    }
  }

  if (decay::nuc_exists(26, 52)) {
    decay::set_nucdecayenergygamma(decay::get_nucindex(26, 52), 0.86 * MEV);  // Fe52
  }
  if (decay::nuc_exists(25, 52)) {
    decay::set_nucdecayenergygamma(decay::get_nucindex(25, 52), 3.415 * MEV);  // Mn52
  }
}

// construct an energy ordered gamma ray line list.
void init_gamma_linelist() {
  read_decaydata();

  // Now do the sorting.

  std::ptrdiff_t total_lines = 0;
  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    total_lines += std::ssize(gamma_spectra[nucindex]);
  }
  printout("total gamma-ray lines %td\n", total_lines);

  allnuc_gamma_line_list = std::vector<NucGammaLine>();
  allnuc_gamma_line_list.reserve(total_lines);

  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    for (std::ptrdiff_t j = 0; j < std::ssize(gamma_spectra[nucindex]); j++) {
      allnuc_gamma_line_list.push_back(
          {.nucindex = nucindex, .nucgammaindex = static_cast<int>(j), .energy = gamma_spectra[nucindex][j].energy});
    }
  }
  allnuc_gamma_line_list.shrink_to_fit();
  assert_always(static_cast<int>(allnuc_gamma_line_list.size()) == total_lines);
  std::ranges::SORT_OR_STABLE_SORT(allnuc_gamma_line_list, [](const NucGammaLine &g1, const NucGammaLine &g2) {
    // true if d1 < d2
    if (g1.energy < g2.energy) {
      return true;
    }
    if (g1.energy == g2.energy && g1.nucindex < g2.nucindex) {
      return true;
    }
    if (g1.energy == g2.energy && g1.nucindex == g2.nucindex && g1.nucgammaindex < g2.nucgammaindex) {
      return true;
    }
    return false;
  });

  FILE *const line_list = fopen_required("gammalinelist.out", "w");

  fprintf(line_list, "#index nucindex Z A nucgammmaindex en_gamma_mev gammaline_probability\n");
  for (std::ptrdiff_t i = 0; i < total_lines; i++) {
    const int nucindex = allnuc_gamma_line_list[i].nucindex;
    const int index = allnuc_gamma_line_list[i].nucgammaindex;
    fprintf(line_list, "%d %d %d %d %d %g %g \n", static_cast<int>(i), allnuc_gamma_line_list[i].nucindex,
            decay::get_nuc_z(allnuc_gamma_line_list[i].nucindex), decay::get_nuc_a(allnuc_gamma_line_list[i].nucindex),
            allnuc_gamma_line_list[i].nucgammaindex, gamma_spectra[nucindex][index].energy / MEV,
            gamma_spectra[nucindex][index].probability);
  }
  fclose(line_list);
}

void init_xcom_photoion_data() {
  // read the file
  printout("reading XCOM photoionization data...\n");
  // reserve memory
  for (int Z = 0; Z < numb_xcom_elements; Z++) {
    photoion_data[Z].reserve(100);
  }
  std::string filepath{"xcom_photoion_data.txt"};
  if (!std::filesystem::exists(filepath)) {
    filepath = "data/xcom_photoion_data.txt";
  }
  assert_always(std::filesystem::exists(filepath));

  std::ifstream data_fs(filepath);
  std::string line_str;
  // now read the file a second time to store the data
  while (getline(data_fs, line_str)) {
    if (line_str[0] != '#') {
      int Z = 0;
      double E = 0;
      double sigma = 0;
      if (3 == std::sscanf(line_str.c_str(), "%d %lg %lg", &Z, &E, &sigma)) {
        assert_always(Z > 0);
        assert_always(Z <= numb_xcom_elements);
        photoion_data[Z - 1].push_back({.energy = E, .sigma_xcom = sigma});
      }
    }
  }
}

__host__ __device__ auto choose_gamma_ray(const int nucindex) -> double {
  // Routine to choose which gamma ray line it'll be.

  const double E_gamma = decay::nucdecayenergygamma(nucindex);  // Average energy per gamma line of a decay

  const double zrand = rng_uniform();
  double runtot = 0.;
  for (ptrdiff_t n = 0; n < std::ssize(gamma_spectra[nucindex]); n++) {
    runtot += gamma_spectra[nucindex][n].probability * gamma_spectra[nucindex][n].energy / E_gamma;
    if (zrand <= runtot) {
      return gamma_spectra[nucindex][n].energy / H;
    }
  }

  printout("Failure to choose line (pellet_nucindex %d). Abort. zrand %g runtot %g\n", nucindex, zrand, runtot);
  assert_always(false);
}

constexpr auto sigma_compton_partial(const double x, const double f_max) -> double
// Routine to compute the partial cross section for Compton scattering.
//   xx is the photon energy (in units of electron mass) and f
//  is the energy loss factor up to which we wish to integrate.
{
  const double term1 = ((x * x) - (2 * x) - 2) * std::log(f_max) / x / x;
  const double term2 = (((f_max * f_max) - 1) / (f_max * f_max)) / 2;
  const double term3 = ((f_max - 1) / x) * ((1 / x) + (2 / f_max) + (1 / (x * f_max)));

  return (3 * SIGMA_T * (term1 + term2 + term3) / (8 * x));
}

auto get_chi_compton_rf(const Packet &pkt) -> double {
  // calculate the absorption coefficient [cm^-1] for Compton scattering in the observer reference frame
  // Start by working out the compton x-section in the co-moving frame.

  const double xx = H * pkt.nu_cmf / ME / CLIGHT / CLIGHT;

  // Use this to decide whether the Thompson limit is acceptable.

  const double sigma_cmf = (xx < THOMSON_LIMIT) ? SIGMA_T : sigma_compton_partial(xx, 1 + (2 * xx));

  // Now need to multiply by the electron number density.
  const double chi_cmf = sigma_cmf * grid::get_nnetot(grid::get_cell_modelgridindex(pkt.where));

  // convert between frames
  const double chi_rf = chi_cmf * calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);

  assert_testmodeonly(std::isfinite(chi_rf));

  return chi_rf;
}

auto choose_f(const double xx, const double zrand) -> double
// To choose the value of f to integrate to - idea is we want
//   sigma_compton_partial(xx,f) = zrand.
{
  double f_max = 1 + (2 * xx);
  double f_min = 1;

  const double norm = zrand * sigma_compton_partial(xx, f_max);

  int count = 0;
  double err = 1e20;

  double ftry = (f_max + f_min) / 2;
  while ((err > 1.e-4) && (count < 1000)) {
    ftry = (f_max + f_min) / 2;
    const double sigma_try = sigma_compton_partial(xx, ftry);
    if (sigma_try > norm) {
      f_max = ftry;
      err = (sigma_try - norm) / norm;
    } else {
      f_min = ftry;
      err = (norm - sigma_try) / norm;
    }
    //      printout("error %g\n",err);
    count++;
    if (count == 1000) {
      printout("Compton hit 1000 tries. %g %g %g %g %g\n", f_max, f_min, ftry, sigma_try, norm);
    }
  }

  return ftry;
}

// For Thomson scattering we can get the new angle from a random number very easily.
auto thomson_angle() -> double {
  const double B_coeff = (8. * rng_uniform()) - 4.;

  const double t_coeff = std::cbrt((std::sqrt((B_coeff * B_coeff) + 4) - B_coeff) / 2);

  const double mu = (1 / t_coeff) - t_coeff;

  assert_always(fabs(mu) <= 1);

  return mu;
}

// scattering a direction through angle theta.
[[nodiscard]] auto scatter_dir(const std::array<double, 3> &dir_in, const double cos_theta) -> std::array<double, 3> {
  // begin with setting the direction in coordinates where original direction
  // is parallel to z-hat.

  const double phi = rng_uniform() * 2 * PI;

  const double sin_theta_sq = 1. - (cos_theta * cos_theta);
  const double sin_theta = std::sqrt(sin_theta_sq);
  const double zprime = cos_theta;
  const double xprime = sin_theta * cos(phi);
  const double yprime = sin_theta * sin(phi);

  // Now need to derotate the coordinates back to real x,y,z.
  // Rotation matrix is determined by dir_in.

  const double norm1 = 1. / std::sqrt((dir_in[0] * dir_in[0]) + (dir_in[1] * dir_in[1]));
  const double norm2 = 1. / vec_len(dir_in);

  const double r11 = dir_in[1] * norm1;
  const double r12 = -1 * dir_in[0] * norm1;
  const double r13 = 0.;
  const double r21 = dir_in[0] * dir_in[2] * norm1 * norm2;
  const double r22 = dir_in[1] * dir_in[2] * norm1 * norm2;
  const double r23 = -1 * norm2 / norm1;
  const double r31 = dir_in[0] * norm2;
  const double r32 = dir_in[1] * norm2;
  const double r33 = dir_in[2] * norm2;

  const auto dir_out = std::array<double, 3>{(r11 * xprime) + (r21 * yprime) + (r31 * zprime),
                                             (r12 * xprime) + (r22 * yprime) + (r32 * zprime),
                                             (r13 * xprime) + (r23 * yprime) + (r33 * zprime)};

  assert_testmodeonly(std::fabs(vec_len(dir_out) - 1.) < 1e-10);

  return dir_out;
}

// handle physical Compton scattering event
void compton_scatter(Packet &pkt) {
  //  printout("Compton scattering.\n");

  const double xx = H * pkt.nu_cmf / ME / CLIGHT / CLIGHT;

  // It is known that a Compton scattering event is going to take place.
  // We need to do two things - (1) decide whether to convert energy
  // to electron or leave as gamma (2) decide properties of new packet.

  // The probability of giving energy to electron is related to the
  // energy change of the gamma ray. This is equivalent to the choice of
  // scattering angle. Probability of scattering into particular angle
  // (i.e. final energy) is related to the partial cross-section.

  // Choose a random number to get the energy. Want to find the
  // factor by which the energy changes "f" such that
  // sigma_partial/sigma_tot = zrand

  // initialise with Thomson limit case (no energy loss)
  double f = 1.;
  bool stay_gamma = true;
  if (xx >= THOMSON_LIMIT) {
    f = choose_f(xx, rng_uniform());

    assert_always(f >= 1.);
    assert_always(f <= (2 * xx + 1.));

    // Prob of keeping gamma ray is...
    const double prob_gamma = 1. / f;

    stay_gamma = (rng_uniform() < prob_gamma);
  }

  if (stay_gamma) {
    // It stays as a gamma ray. Change frequency and direction in
    // co-moving frame then transfer back to rest frame.

    pkt.nu_cmf = pkt.nu_cmf / f;  // reduce frequency

    // The packet has stored the direction in the rest frame.
    // Use aberration of angles to get this into the co-moving frame.

    auto vel_vec = get_velocity(pkt.pos, pkt.prop_time);

    auto cmf_dir = angle_ab(pkt.dir, vel_vec);

    // Now change the direction through the scattering angle.

    const double cos_theta = (xx < THOMSON_LIMIT) ? thomson_angle() : 1. - ((f - 1) / xx);

    const auto new_dir = scatter_dir(cmf_dir, cos_theta);

    assert_testmodeonly(fabs(1. - dot(new_dir, new_dir)) < 1e-8);
    assert_testmodeonly(fabs(dot(new_dir, cmf_dir) - cos_theta) < 1e-8);

    // Now convert back again.

    pkt.dir = angle_ab(new_dir, vec_scale(vel_vec, -1.));

    assert_testmodeonly(std::fabs(vec_len(pkt.dir) - 1.) < 1e-10);

    // It now has a rest frame direction and a co-moving frequency.
    //  Just need to set the rest frame energy.
    const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
    pkt.nu_rf = pkt.nu_cmf / dopplerfactor;
    pkt.e_rf = pkt.e_cmf / dopplerfactor;

    pkt.last_cross = BOUNDARY_NONE;  // allow it to re-cross a boundary
  } else {
    // energy loss of the gamma becomes energy of the electron (needed to calculate time-dependent thermalisation rate)
    if constexpr (PARTICLE_THERMALISATION_SCHEME == ThermalisationScheme::DETAILEDWITHGAMMAPRODUCTS) {
      pkt.nu_cmf = pkt.nu_cmf * (1 - 1 / f);
      pkt.type = TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS;
    } else {
      pkt.type = TYPE_NTLEPTON_DEPOSITED;
    }
    pkt.absorptiontype = -3;
    stats::increment(stats::COUNTER_NT_STAT_FROM_GAMMA);
  }
}

// calculate the absorption coefficient [cm^-1] for photo electric effect scattering in the observer reference frame
auto get_chi_photo_electric_rf(const Packet &pkt) -> double {
  double chi_cmf{NAN};
  // Start by working out the x-section in the co-moving frame.

  const int mgi = grid::get_cell_modelgridindex(pkt.where);

  if (mgi >= grid::get_npts_model()) {
    // empty cell
    return 0.;
  }

  const double rho = grid::get_rho(mgi);

  if (globals::gamma_kappagrey < 0) {
    chi_cmf = 0.;
    if (!USE_XCOM_GAMMAPHOTOION) {
      // Cross sections from Equation 2 of Ambwani & Sutherland (1988), attributed to Veigele (1973)

      // 2.41326e19 Hz = 100 keV / H
      const double hnu_over_100kev = pkt.nu_cmf / 2.41326e+19;

      // double sigma_cmf_cno = 0.0448e-24 * pow(hnu_over_100kev, -3.2);

      const double sigma_cmf_si = 1.16e-24 * pow(hnu_over_100kev, -3.13);

      const double sigma_cmf_fe = 25.7e-24 * pow(hnu_over_100kev, -3.0);

      // Now need to multiply by the particle number density.

      const double chi_cmf_si = sigma_cmf_si * (rho / MH / 28);
      // Assumes Z = 14. So mass = 28.

      const double chi_cmf_fe = sigma_cmf_fe * (rho / MH / 56);
      // Assumes Z = 28. So mass = 56.

      const double f_fe = grid::get_ffegrp(mgi);

      chi_cmf = (chi_cmf_fe * f_fe) + (chi_cmf_si * (1. - f_fe));
    } else {
      const double hnu_over_1MeV = pkt.nu_cmf / 2.41326e+20;
      const double log10_hnu_over_1MeV = log10(hnu_over_1MeV);
      for (int i = 0; i < get_nelements(); i++) {
        // determine charge number:
        const int Z = get_atomicnumber(i);
        auto numb_energies = std::ssize(photoion_data[Z - 1]);
        if (numb_energies == 0) {
          continue;
        }
        const double n_i = grid::get_elem_numberdens(mgi, i);  // number density in the current cell
        if (n_i == 0) {
          continue;
        }
        // get indices of lower and upper boundary
        int E_gtr_idx = -1;

        for (int j = 0; j < numb_energies; j++) {
          if (photoion_data[Z - 1][j].energy > hnu_over_1MeV) {
            E_gtr_idx = j;
            break;
          }
        }
        if (E_gtr_idx == 0) {  // packet energy smaller than all tabulated values
          chi_cmf += photoion_data[Z - 1][0].sigma_xcom * n_i;
          continue;
        }
        if (E_gtr_idx == -1) {  // packet energy greater than all tabulated values
          chi_cmf += photoion_data[Z - 1][numb_energies - 1].sigma_xcom * n_i;
          continue;
        }
        assert_always(E_gtr_idx > 0);
        assert_always(E_gtr_idx < numb_energies);
        const int E_smaller_idx = E_gtr_idx - 1;
        assert_always(E_smaller_idx >= 0);
        const double log10_E = log10_hnu_over_1MeV;
        const double log10_E_gtr = log10(photoion_data[Z - 1][E_gtr_idx].energy);
        const double log10_E_smaller = log10(photoion_data[Z - 1][E_smaller_idx].energy);
        const double log10_sigma_lower = log10(photoion_data[Z - 1][E_smaller_idx].sigma_xcom);
        const double log10_sigma_gtr = log10(photoion_data[Z - 1][E_gtr_idx].sigma_xcom);
        // interpolate or extrapolate, both linear in log10-log10 space
        const double log10_intpol = log10_E_smaller + ((log10_sigma_gtr - log10_sigma_lower) /
                                                       (log10_E_gtr - log10_E_smaller) * (log10_E - log10_E_smaller));
        const double sigma_intpol = pow(10., log10_intpol) * 1.0e-24;  // now in cm^2
        const double chi_cmf_contrib = sigma_intpol * n_i;
        chi_cmf += chi_cmf_contrib;
      }
    }
  } else {
    chi_cmf = globals::gamma_kappagrey * rho;
  }

  // Now convert between frames.

  const double chi_rf = chi_cmf * calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
  return chi_rf;
}

// calculate the absorption coefficient [cm^-1] for pair production in the observer reference frame
auto get_chi_pair_prod_rf(const Packet &pkt) -> double {
  const int mgi = grid::get_cell_modelgridindex(pkt.where);
  const double rho = grid::get_rho(mgi);

  if (globals::gamma_kappagrey >= 0.) {
    return 0.;
  }

  // 2.46636e+20 Hz = 1022 keV / H
  if (pkt.nu_cmf <= 2.46636e+20) {
    return 0.;
  }

  // double sigma_cmf_cno;
  double sigma_cmf_si{NAN};
  double sigma_cmf_fe{NAN};
  const double f_fe = grid::get_ffegrp(mgi);

  // Cross sections from Equation 2 of Ambwani & Sutherland (1988), attributed to Hubbell (1969)

  // 3.61990e+20 = 1500 keV in frequency / H
  const double hnu_over_mev = pkt.nu_cmf / 2.41326e+20;
  if (pkt.nu_cmf > 3.61990e+20) {
    // sigma_cmf_cno = (0.0481 + (0.301 * (hnu_over_mev - 1.5))) * 49.e-27;

    sigma_cmf_si = (0.0481 + (0.301 * (hnu_over_mev - 1.5))) * 196.e-27;

    sigma_cmf_fe = (0.0481 + (0.301 * (hnu_over_mev - 1.5))) * 784.e-27;
  } else {
    // sigma_cmf_cno = 1.0063 * (hnu_over_mev - 1.022) * 49.e-27;

    sigma_cmf_si = 1.0063 * (hnu_over_mev - 1.022) * 196.e-27;

    sigma_cmf_fe = 1.0063 * (hnu_over_mev - 1.022) * 784.e-27;
  }

  // multiply by the particle number density.

  // sigma_cmf_cno *= rho * (1. - f_fe) / MH / 14;
  // Assumes Z = 7. So mass = 14.

  const double chi_cmf_si = sigma_cmf_si * (rho / MH / 28);
  // Assumes Z = 14. So mass = 28.

  const double chi_cmf_fe = sigma_cmf_fe * (rho / MH / 56);
  // Assumes Z = 28. So mass = 56.

  const double chi_cmf = (chi_cmf_fe * f_fe) + (chi_cmf_si * (1. - f_fe));

  // Now need to convert between frames.

  double chi_rf = chi_cmf * calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);

  if (chi_rf < 0) {
    printout("Negative pair production sigma. Setting to zero. Abort? %g\n", chi_rf);
    chi_rf = 0.;
  }

  return chi_rf;
}

// Routine to compute the mean energy converted to non-thermal electrons times the Klein-Nishina cross section.
constexpr auto meanf_sigma(const double x) -> double {
  const double f = 1 + (2 * x);

  const double term0 = 2 / x;
  const double term1 = (1 - (2 / x) - (3 / (x * x))) * log(f);
  const double term2 = ((4 / x) + (3 / (x * x)) - 1) * 2 * x / f;
  const double term3 = (1 - (2 / x) - (1 / (x * x))) * 2 * x * (1 + x) / f / f;
  const double term4 = -2. * x * ((4 * x * x) + (6 * x) + 3) / 3 / f / f / f;

  const double tot = 3 * SIGMA_T * (term0 + term1 + term2 + term3 + term4) / (8 * x);

  return tot;
}

// Subroutine to record the heating rate in a cell due to gamma rays.
// By heating rate I mean, for now, really the rate at which the code is making
// k-packets in that cell which will then convert into r-packets. This is (going
// to be) used for the new light_curve syn-style calculation.
// The intention is that dep_estimator_gamma will contain the emissivity of r-packets
// in the co-moving frame (which is going to be isotropic).
void update_gamma_dep(const Packet &pkt, const double dist, const int mgi, const int nonemptymgi) {
  if (!(dist > 0)) {
    return;
  }
  if constexpr (PARTICLE_THERMALISATION_SCHEME == ThermalisationScheme::DETAILEDWITHGAMMAPRODUCTS) {
    return;  // don't instantly deposit energy from gamma rays, handle the particles they produce instead
  }

  const double doppler_sq = doppler_squared_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);

  const double xx = H * pkt.nu_cmf / ME / CLIGHT / CLIGHT;
  double heating_cont = ((meanf_sigma(xx) * grid::get_nnetot(mgi)) + get_chi_photo_electric_rf(pkt) +
                         (get_chi_pair_prod_rf(pkt) * (1. - (2.46636e+20 / pkt.nu_cmf))));
  heating_cont = heating_cont * pkt.e_rf * dist * doppler_sq;

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

// handle physical pair production event
//
//  In pair production, the original gamma makes an electron positron pair - kinetic energy equal to
//  gamma ray energy - 1.022 MeV. We assume that the electron deposits any kinetic energy directly to
//  the thermal pool. The positron annihilates with an electron locally making a pair of gamma rays
//  at 0.511 MeV in the local cmf (isotropic). So all the thermal energy goes to the thermal pool
//  immediately and the remainder goes into gamma-rays at 0.511 MeV.
void pair_prod(Packet &pkt) {
  const double prob_gamma = 1.022 * MEV / (H * pkt.nu_cmf);

  assert_always(prob_gamma >= 0);

  const auto zrand = rng_uniform();
  if (zrand > prob_gamma) {
    if constexpr (PARTICLE_THERMALISATION_SCHEME == ThermalisationScheme::DETAILEDWITHGAMMAPRODUCTS) {
      // Convert it to an e-minus or positron kinetic energy packet
      pkt.type = (rng_uniform() > 0.5) ? TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS : TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS;
    } else {
      pkt.type = TYPE_NTLEPTON_DEPOSITED;
    }

    // nu_cmf stays the same as the gamma energy becomes the kinetic energy of the electron
    pkt.absorptiontype = -5;
    stats::increment(stats::COUNTER_NT_STAT_FROM_GAMMA);
  } else {
    // The energy goes into emission at 511 keV.
    pkt.nu_cmf = 0.511 * MEV / H;

    // Now let's give the gamma ray a direction.

    const auto dir_cmf = get_rand_isotropic_unitvec();

    // This direction is in the cmf - we want to convert it to the rest
    // frame - use aberration of angles. We want to convert from cmf to
    // rest so need -ve velocity.

    const auto vel_vec = get_velocity(pkt.pos, -1. * pkt.prop_time);
    // negative time since we want the backwards transformation here

    pkt.dir = angle_ab(dir_cmf, vel_vec);

    const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
    pkt.nu_rf = pkt.nu_cmf / dopplerfactor;
    pkt.e_rf = pkt.e_cmf / dopplerfactor;

    pkt.type = TYPE_GAMMA;
    pkt.last_cross = BOUNDARY_NONE;
  }
}

// move a gamma packet until time t2
void transport_gamma(Packet &pkt, const double t2) {
  // Assign optical depth to next physical event. And start counter of
  // optical depth for this path.
  const double zrand = rng_uniform_pos();
  const double tau_next = -1. * log(zrand);
  const double tau_current = 0.;

  // Start by finding the distance to the crossing of the grid cell
  // boundaries. sdist is the boundary distance and snext is the
  // grid cell into which we pass.

  auto [sdist, snext] = grid::boundary_distance(pkt.dir, pkt.pos, pkt.prop_time, pkt.where, &pkt.last_cross);

  const double maxsdist = (GRID_TYPE == GridType::CARTESIAN3D)
                              ? globals::rmax * pkt.prop_time / globals::tmin
                              : 2 * globals::rmax * (pkt.prop_time + sdist / CLIGHT_PROP) / globals::tmin;
  if (sdist > maxsdist) {
    printout("Unreasonably large sdist (gamma). Abort. %g %g %g\n", globals::rmax, pkt.prop_time / globals::tmin,
             sdist);
    assert_always(false);
  }

  if (sdist < 0) {
    printout("Negative distance (sdist). Abort?\n");
    sdist = 0;
  }

  if (((snext < 0) && (snext != -99)) || (snext >= grid::ngrid)) {
    printout("Heading for inappropriate grid cell. Abort.\n");
    printout("Current cell %d, target cell %d.\n", pkt.where, snext);
    assert_always(false);
  }

  if (sdist > globals::max_path_step) {
    sdist = globals::max_path_step;
    snext = pkt.where;
  }

  // Now consider the scattering/destruction processes.
  // Compton scattering - need to determine the scattering co-efficient.
  // Routine returns the value in the rest frame.

  double chi_compton = 0.;
  if (globals::gamma_kappagrey < 0) {
    chi_compton = get_chi_compton_rf(pkt);
  }

  const double chi_photo_electric = get_chi_photo_electric_rf(pkt);
  const double chi_pair_prod = get_chi_pair_prod_rf(pkt);
  const double chi_tot = chi_compton + chi_photo_electric + chi_pair_prod;

  assert_testmodeonly(std::isfinite(chi_compton));
  assert_testmodeonly(std::isfinite(chi_photo_electric));
  assert_testmodeonly(std::isfinite(chi_pair_prod));

  // So distance before physical event is...

  const double edist = chi_tot > 0. ? (tau_next - tau_current) / chi_tot : std::numeric_limits<double>::max();

  assert_always(edist >= 0);

  // Find how far it can travel during the time interval.

  const double tdist = (t2 - pkt.prop_time) * CLIGHT_PROP;

  assert_always(tdist >= 0);

  if ((sdist < tdist) && (sdist < edist)) {
    move_pkt_withtime(pkt, sdist / 2.);

    // Move it into the new cell.
    const int mgi = grid::get_cell_modelgridindex(pkt.where);
    const int nonemptymgi = (mgi < grid::get_npts_model()) ? grid::get_modelcell_nonemptymgi(mgi) : -1;
    if (chi_tot > 0 && nonemptymgi >= 0) {
      update_gamma_dep(pkt, sdist, mgi, nonemptymgi);
    }

    move_pkt_withtime(pkt, sdist / 2.);

    if (snext != pkt.where) {
      grid::change_cell(pkt, snext);
    }
  } else if ((tdist < sdist) && (tdist < edist)) {
    // Doesn't reach boundary.
    move_pkt_withtime(pkt, tdist / 2.);
    const int mgi = grid::get_cell_modelgridindex(pkt.where);
    const int nonemptymgi = (mgi < grid::get_npts_model()) ? grid::get_modelcell_nonemptymgi(mgi) : -1;

    if (chi_tot > 0 && nonemptymgi >= 0) {
      update_gamma_dep(pkt, tdist, mgi, nonemptymgi);
    }
    move_pkt_withtime(pkt, tdist / 2.);
    pkt.prop_time = t2;  // prevent roundoff error
  } else if ((edist < sdist) && (edist < tdist)) {
    move_pkt_withtime(pkt, edist / 2.);
    const int mgi = grid::get_cell_modelgridindex(pkt.where);
    const int nonemptymgi = (mgi < grid::get_npts_model()) ? grid::get_modelcell_nonemptymgi(mgi) : -1;
    if (chi_tot > 0 && nonemptymgi >= 0) {
      update_gamma_dep(pkt, edist, mgi, nonemptymgi);
    }
    move_pkt_withtime(pkt, edist / 2.);

    // event occurs. Choose which event and call the appropriate subroutine.
    const double chi_rnd = rng_uniform() * chi_tot;
    if (chi_compton > chi_rnd) {
      // Compton scattering.
      compton_scatter(pkt);
    } else if ((chi_compton + chi_photo_electric) > chi_rnd) {
      // Photo electric effect
      if constexpr (PARTICLE_THERMALISATION_SCHEME == ThermalisationScheme::DETAILEDWITHGAMMAPRODUCTS) {
        pkt.type = TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS;
        // nu_cmf stays the same as the gamma-ray energy becomes the kinetic energy of the electron (minus ionisation
        // energy but this is neglected here)
      } else {
        pkt.type = TYPE_NTLEPTON_DEPOSITED;
      }

      pkt.absorptiontype = -4;
      stats::increment(stats::COUNTER_NT_STAT_FROM_GAMMA);
    } else {
      // It's a pair production
      pair_prod(pkt);
    }
  } else {
    assert_always(false);
  }
}

void barnes_thermalisation(Packet &pkt)
// Barnes treatment: packet is either getting absorbed immediately and locally
// creating a k-packet or it escapes. The absorption probability matches the
// Barnes thermalization efficiency, for expressions see the original paper:
// https://ui.adsabs.harvard.edu/abs/2016ApJ...829..110B
{
  // compute thermalization efficiency (= absorption probability)
  // 0.1 is an average value to fit the analytic approximations from the paper.
  // Alternative: Distinguish between low-E (kappa = 1) or high-E (kappa = 0.05)
  // packets.
  // constexpr double mean_gamma_opac = 0.1;

  // determine average initial density via kinetic energy
  const double E_kin = grid::get_ejecta_kinetic_energy();
  const double v_ej = sqrt(E_kin * 2 / grid::mtot_input);

  // const double t_ineff = sqrt(rho_0 * R_0 * pow(t_0, 2) * mean_gamma_opac);
  const double t_ineff = 1.4 * 86400. * sqrt(grid::mtot_input / (5.e-3 * 1.989 * 1.e33)) * ((0.2 * 29979200000) / v_ej);
  const double tau = pow(t_ineff / pkt.prop_time, 2.);
  const double f_gamma = 1. - exp(-tau);
  assert_always(f_gamma >= 0.);
  assert_always(f_gamma <= 1.);

  // either absorb packet or let it escape
  if (rng_uniform() < f_gamma) {
    // packet is absorbed and contributes to the heating as a k-packet
    pkt.type = TYPE_NTLEPTON_DEPOSITED;
    pkt.absorptiontype = -4;
  } else {
    // let packet escape, i.e. make it inactive
    pkt.type = TYPE_ESCAPE;
    grid::change_cell(pkt, -99);
  }
}

void wollaeger_thermalisation(Packet &pkt) {
  // corresponds to a local version of the Barnes scheme, i.e. it takes into account the local mass
  // density rather than a value averaged over the ejecta
  constexpr double mean_gamma_opac = 0.1;
  // integration: requires distances within single cells in radial direction and the corresponding densities
  // need to create a packet copy which is moved during the integration
  Packet pkt_copy = pkt;
  pkt.dir = vec_norm(pkt_copy.pos);
  const double t_current = pkt.prop_time;
  double tau = 0.;
  bool end_packet = false;
  while (!end_packet) {
    // distance to the next cell
    const auto [sdist, snext] =
        grid::boundary_distance(pkt_copy.dir, pkt_copy.pos, pkt_copy.prop_time, pkt_copy.where, &pkt_copy.last_cross);
    const double s_cont = sdist * t_current * t_current * t_current / std::pow(pkt_copy.prop_time, 3);
    const int mgi = grid::get_cell_modelgridindex(pkt_copy.where);
    if (mgi != grid::get_npts_model()) {
      tau += grid::get_rho(mgi) * s_cont * mean_gamma_opac;  // contribution to the integral
    }
    // move packet copy now
    move_pkt_withtime(pkt_copy, sdist);

    grid::change_cell(pkt_copy, snext);
    end_packet = (pkt_copy.type == TYPE_ESCAPE);
  }
  const double f_gamma = 1. - std::exp(-tau);
  assert_always(f_gamma >= 0.);
  assert_always(f_gamma <= 1.);

  // either absorb packet or let it escape
  if (rng_uniform() < f_gamma) {
    // packet is absorbed and contributes to the heating as a k-packet
    pkt.type = TYPE_NTLEPTON_DEPOSITED;
    pkt.absorptiontype = -4;

  } else {
    // let packet escape, i.e. make it inactive
    pkt.type = TYPE_ESCAPE;
    grid::change_cell(pkt, -99);
  }
}

void guttman_thermalisation(Packet &pkt) {
  // Guttman+2024, arXiv:2403.08769v1
  // extension of the Wollaeger scheme. Rather than calculating a single optical depth in radial outward
  // direction, it calculates a spherical average in all possible gamma-ray emission directions.

  // toy value for the mean gamma ray opacity for a start, the paper discuss this some more detail, too
  constexpr double mean_gamma_opac =
      0.03;  // mean gamma opacity from section 3.2, taken the lower one due to almost symmetrical matter at late times

  const double t_0 = globals::tmin;
  const double t = pkt.prop_time;

  // calculate average optical w.r.t. to emission direction
  // discretize the two sphere into octants, i.e. four values in phi and two in theta
  constexpr int numb_rnd_dirs = 100;
  auto column_densities = std::array<double, numb_rnd_dirs>{};
  for (int i = 0; i < numb_rnd_dirs; i++) {
    // compute column density by moving an artificial packet outwards and integrating over local density
    // WARNING: This simple implementation relies on a relatively large number of random directions
    // as it assumes that the directions are distributed equally which is likely not the case with low statistics.

    // step 1: draw a random direction
    Packet pkt_copy = pkt;
    // phi rotation: around z-axis
    const auto random_dir = get_rand_isotropic_unitvec();
    pkt_copy.dir = random_dir;  // fix new direction

    // step 2: move packet into the calculated direction and integrate the density
    bool end_packet = false;
    while (!end_packet) {
      // distance to the next cell
      const auto [sdist, snext] =
          grid::boundary_distance(pkt_copy.dir, pkt_copy.pos, pkt_copy.prop_time, pkt_copy.where, &pkt_copy.last_cross);
      const double s_cont = sdist * std::pow(t, 3.) / std::pow(pkt_copy.prop_time, 3.);
      const int mgi = grid::get_cell_modelgridindex(pkt_copy.where);
      if (mgi != grid::get_npts_model()) {
        column_densities[i] += grid::get_rho_tmin(mgi) * s_cont;  // contribution to the integral
      }
      // move packet copy now
      move_pkt_withtime(pkt_copy, sdist);

      grid::change_cell(pkt_copy, snext);
      end_packet = (pkt_copy.type == TYPE_ESCAPE);
    }
  }
  const double avg_column_density =
      std::accumulate(column_densities.cbegin(), column_densities.cend(), 0.) / std::ssize(column_densities);
  const double t_gamma = sqrt(mean_gamma_opac * avg_column_density * t_0 * t_0);

  // compute the (discretized) integral
  double f_gamma = 0.;
  const double width = 4 * PI / numb_rnd_dirs;
  for (int i = 0; i < numb_rnd_dirs; i++) {
    const double summand =
        width * (1 - std::exp(-std::pow(t_gamma, 2.) / std::pow(t, 2.) * column_densities[i] / avg_column_density));
    printout("width: %f t_gamma: %f t: %f column_densities[i]: %f avg_column_density: %f summand: %f", width, t_gamma,
             t, column_densities[i], avg_column_density, summand);
    f_gamma += summand;
  }
  f_gamma /= (4 * PI);

  assert_always(f_gamma >= 0.);
  assert_always(f_gamma <= 1.);

  // either absorb packet or let it escape
  if (rng_uniform() < f_gamma) {
    // packet is absorbed and contributes to the heating as a k-packet
    pkt.type = TYPE_NTLEPTON_DEPOSITED;
    pkt.absorptiontype = -4;
  } else {
    // let packet escape, i.e. make it inactive
    pkt.type = TYPE_ESCAPE;
    grid::change_cell(pkt, -99);
  }
}

}  // anonymous namespace

void init_gamma_data() {
  init_gamma_linelist();
  if constexpr (USE_XCOM_GAMMAPHOTOION) {
    init_xcom_photoion_data();
  }
}

// convert a pellet to a gamma ray (or kpkt if no gamma spec loaded)
__host__ __device__ void pellet_gamma_decay(Packet &pkt) {
  // Start by getting the position of the pellet at the point of decay. Pellet
  // is moving with the matter.

  // if no gamma spectra is known, then covert straight to kpkts (e.g., Fe52, Mn52)
  if (gamma_spectra[pkt.pellet_nucindex].empty()) {
    pkt.type = TYPE_KPKT;
    pkt.absorptiontype = -6;
    return;
  }

  // Now let's give the gamma ray a direction.
  // Assuming isotropic emission in cmf

  const auto dir_cmf = get_rand_isotropic_unitvec();

  // This direction is in the cmf - we want to convert it to the rest
  // frame - use aberration of angles. We want to convert from cmf to
  // rest so need -ve velocity.

  const auto vel_vec = get_velocity(pkt.pos, -1. * pkt.tdecay);
  // negative time since we want the backwards transformation here

  pkt.dir = angle_ab(dir_cmf, vel_vec);

  // Now need to assign the frequency of the packet in the co-moving frame.

  pkt.nu_cmf = choose_gamma_ray(pkt.pellet_nucindex);

  // Finally we want to put in the rest frame energy and frequency. And record
  // that it's now a gamma ray.

  pkt.prop_time = pkt.tdecay;
  const double dopplerfactor = calculate_doppler_nucmf_on_nurf(pkt.pos, pkt.dir, pkt.prop_time);
  pkt.nu_rf = pkt.nu_cmf / dopplerfactor;
  pkt.e_rf = pkt.e_cmf / dopplerfactor;

  pkt.type = TYPE_GAMMA;
  pkt.last_cross = BOUNDARY_NONE;

  // initialise polarisation information
  pkt.stokes = {1., 0., 0.};

  pkt.pol_dir = cross_prod(pkt.dir, std::array<double, 3>{0., 0., 1.});
  if ((dot(pkt.pol_dir, pkt.pol_dir)) < 1.e-8) {
    pkt.pol_dir = cross_prod(pkt.dir, std::array<double, 3>{0., 1., 0.});
  }

  pkt.pol_dir = vec_norm(pkt.pol_dir);
}

__host__ __device__ void do_gamma(Packet &pkt, const int nts, const double t2) {
  if constexpr (GAMMA_THERMALISATION_SCHEME == ThermalisationScheme::DETAILED) {
    transport_gamma(pkt, t2);
  } else if constexpr (GAMMA_THERMALISATION_SCHEME == ThermalisationScheme::BARNES) {
    barnes_thermalisation(pkt);
  } else if constexpr (GAMMA_THERMALISATION_SCHEME == ThermalisationScheme::WOLLAEGER) {
    wollaeger_thermalisation(pkt);
  } else if constexpr (GAMMA_THERMALISATION_SCHEME == ThermalisationScheme::GUTTMAN) {
    guttman_thermalisation(pkt);
  } else {
    __builtin_unreachable();
  }

  if (pkt.type != TYPE_GAMMA && pkt.type != TYPE_ESCAPE) {
    if constexpr (PARTICLE_THERMALISATION_SCHEME != ThermalisationScheme::DETAILEDWITHGAMMAPRODUCTS) {
      atomicadd(globals::timesteps[nts].gamma_dep_discrete, pkt.e_cmf);
    }

    if constexpr (GAMMA_THERMALISATION_SCHEME != ThermalisationScheme::DETAILED &&
                  GAMMA_THERMALISATION_SCHEME != ThermalisationScheme::DETAILEDWITHGAMMAPRODUCTS) {
      // no transport, so the path-based gamma deposition estimator won't get updated unless we do it here
      const int mgi = grid::get_cell_modelgridindex(pkt.where);
      const int nonemptymgi = grid::get_modelcell_nonemptymgi(mgi);
      atomicadd(globals::dep_estimator_gamma[nonemptymgi], pkt.e_cmf);
    }
  }
}

}  // namespace gammapkt
