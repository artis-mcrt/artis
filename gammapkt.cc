#include "gammapkt.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <vector>

#include "decay.h"
#include "grid.h"
#include "nonthermal.h"
#include "packet.h"
#include "sn3d.h"
#include "stats.h"
#include "vectors.h"

namespace gammapkt {
// Code for handing gamma rays - creation and propagation

struct gamma_spec {
  std::unique_ptr<double[]> energy = nullptr;  // in erg
  std::unique_ptr<double[]> probability = nullptr;
  int nlines = 0;
};

static struct gamma_spec *gamma_spectra;

struct gammaline {
  int nucindex;       // is it a Ni56, Co56, a fake line, etc
  int nucgammaindex;  // which of the lines of that nuclide is it
  double energy;      // in erg

  auto operator<(const struct gammaline &g2) const -> bool {
    // true if d1 < d2
    if (energy < g2.energy) {
      return true;
    }
    if (energy == g2.energy && nucindex < g2.nucindex) {
      return true;
    }
    if (energy == g2.energy && nucindex == g2.nucindex && nucgammaindex < g2.nucgammaindex) {
      return true;
    }
    return false;
  }
};

static std::vector<struct gammaline> allnuc_gamma_line_list;

static void read_gamma_spectrum(const int nucindex, const char filename[50])
// reads in gamma_spectra and returns the average energy in gamma rays per nuclear decay
{
  printout("reading gamma spectrum for Z=%d A=%d from %s...", decay::get_nuc_z(nucindex), decay::get_nuc_a(nucindex),
           filename);

  FILE *filein = fopen_required(filename, "r");
  int nlines = 0;
  assert_always(fscanf(filein, "%d", &nlines) == 1);

  gamma_spectra[nucindex].nlines = nlines;

  gamma_spectra[nucindex].energy = std::make_unique<double[]>(nlines);
  gamma_spectra[nucindex].probability = std::make_unique<double[]>(nlines);

  double E_gamma_avg = 0.;
  for (int n = 0; n < nlines; n++) {
    double en_mev = 0.;
    double prob = 0.;
    assert_always(fscanf(filein, "%lg %lg", &en_mev, &prob) == 2);
    gamma_spectra[nucindex].energy[n] = en_mev * MEV;
    gamma_spectra[nucindex].probability[n] = prob;
    E_gamma_avg += en_mev * MEV * prob;
  }
  fclose(filein);

  decay::set_nucdecayenergygamma(nucindex, E_gamma_avg);

  printout("nlines %d avg_en_gamma %g MeV\n", nlines, E_gamma_avg / MEV);
}

static void set_trivial_gamma_spectrum(const int nucindex) {
  // printout("Setting trivial gamma spectrum for z %d a %d engamma %g\n", z, a, decay::nucdecayenergygamma(z, a));
  const int nlines = 1;
  gamma_spectra[nucindex].nlines = nlines;
  gamma_spectra[nucindex].energy = std::make_unique<double[]>(nlines);
  gamma_spectra[nucindex].probability = std::make_unique<double[]>(nlines);
  gamma_spectra[nucindex].energy[0] = decay::nucdecayenergygamma(nucindex);
  gamma_spectra[nucindex].probability[0] = 1.;
}

static void read_decaydata() {
  // migrate from old filename
  if (!std::ifstream("ni56_lines.txt") && std::ifstream("ni_lines.txt")) {
    printout("Moving ni_lines.txt to ni56_lines.txt\n");
    std::rename("ni_lines.txt", "ni56_lines.txt");
  }

  // migrate from old filename
  if (!std::ifstream("co56_lines.txt") && std::ifstream("co_lines.txt")) {
    printout("Moving co_lines.txt to co56_lines.txt\n");
    std::rename("co_lines.txt", "co56_lines.txt");
  }

  gamma_spectra = static_cast<struct gamma_spec *>(calloc(decay::get_num_nuclides(), sizeof(struct gamma_spec)));

  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    gamma_spectra[nucindex].nlines = 0;
    gamma_spectra[nucindex].energy = nullptr;
    gamma_spectra[nucindex].probability = nullptr;
    const int z = decay::get_nuc_z(nucindex);
    const int a = decay::get_nuc_a(nucindex);
    if (z < 1) {
      continue;
    }

    const char *elname = decay::get_elname(z);
    const size_t elnamelen = strlen(elname);  // excluding the NULL terminator
    assert_always(elnamelen < 7);
    char elnamelower[8];
    for (size_t i = 0; i < elnamelen; i++) {
      elnamelower[i] = static_cast<char>(tolower(elname[i]));
    }
    elnamelower[elnamelen] = '\0';

    // look in the current folder
    char filename[MAXFILENAMELENGTH];
    snprintf(filename, MAXFILENAMELENGTH, "%s%d_lines.txt", elnamelower, a);

    // look in the 'data' subfolder
    char filename2[MAXFILENAMELENGTH];
    snprintf(filename2, MAXFILENAMELENGTH, "data/%s%d_lines.txt", elnamelower, a);

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

  int total_lines = 0;
  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    total_lines += gamma_spectra[nucindex].nlines;
  }
  printout("total gamma-ray lines %d\n", total_lines);

  allnuc_gamma_line_list = std::vector<struct gammaline>();
  allnuc_gamma_line_list.reserve(total_lines);

  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    for (int j = 0; j < gamma_spectra[nucindex].nlines; j++) {
      allnuc_gamma_line_list.push_back(
          {.nucindex = nucindex, .nucgammaindex = j, .energy = gamma_spectra[nucindex].energy[j]});
    }
  }
  allnuc_gamma_line_list.shrink_to_fit();
  assert_always(static_cast<int>(allnuc_gamma_line_list.size()) == total_lines);
  std::sort(allnuc_gamma_line_list.begin(), allnuc_gamma_line_list.end());

  FILE *const line_list = fopen_required("gammalinelist.out", "w");

  fprintf(line_list, "#index nucindex Z A nucgammmaindex en_gamma_mev gammaline_probability\n");
  for (int i = 0; i < total_lines; i++) {
    const int nucindex = allnuc_gamma_line_list[i].nucindex;
    const int index = allnuc_gamma_line_list[i].nucgammaindex;
    fprintf(line_list, "%d %d %d %d %d %g %g \n", i, allnuc_gamma_line_list[i].nucindex,
            decay::get_nuc_z(allnuc_gamma_line_list[i].nucindex), decay::get_nuc_a(allnuc_gamma_line_list[i].nucindex),
            allnuc_gamma_line_list[i].nucgammaindex, gamma_spectra[nucindex].energy[index] / MEV,
            gamma_spectra[nucindex].probability[index]);
  }
  fclose(line_list);
}

void normalise_grey(int nts) {
  const double dt = globals::timesteps[nts].width;
  globals::timesteps[nts].gamma_dep_pathint = 0.;
  for (int mgi = 0; mgi < grid::get_npts_model(); mgi++) {
    if (grid::get_numassociatedcells(mgi) > 0) {
      const double dV = grid::get_modelcell_assocvolume_tmin(mgi) * pow(globals::timesteps[nts].mid / globals::tmin, 3);

      globals::timesteps[nts].gamma_dep_pathint += globals::rpkt_emiss[mgi] / globals::nprocs;

      globals::rpkt_emiss[mgi] = globals::rpkt_emiss[mgi] * ONEOVER4PI / dV / dt / globals::nprocs;

      assert_testmodeonly(globals::rpkt_emiss[mgi] >= 0.);
      assert_testmodeonly(isfinite(globals::rpkt_emiss[mgi]));
    }
  }
}

static void choose_gamma_ray(struct packet *pkt_ptr) {
  // Routine to choose which gamma ray line it'll be.

  const int nucindex = pkt_ptr->pellet_nucindex;
  double const E_gamma = decay::nucdecayenergygamma(nucindex);  // Average energy per gamma line of a decay

  const double zrand = rng_uniform();
  int nselected = -1;
  double runtot = 0.;
  for (int n = 0; n < gamma_spectra[nucindex].nlines; n++) {
    runtot += gamma_spectra[nucindex].probability[n] * gamma_spectra[nucindex].energy[n] / E_gamma;
    if (zrand <= runtot) {
      nselected = n;
      break;
    }
  }

  if (nselected < 0) {
    printout("Failure to choose line (packet type %d pellet_nucindex %d). Abort. zrand %g runtot %g\n", pkt_ptr->type,
             pkt_ptr->pellet_nucindex, zrand, runtot);
    abort();
  }

  pkt_ptr->nu_cmf = gamma_spectra[nucindex].energy[nselected] / H;
  // printout("%s PELLET %g\n", gammaspec->filename, gammaspec->energy[nselected]);
}

void pellet_gamma_decay(struct packet *pkt_ptr) {
  // Subroutine to convert a pellet to a gamma ray (or kpkt if no gamma spec loaded)

  // pkt_ptr is a pointer to the packet that is decaying.

  // Start by getting the position of the pellet at the point of decay. Pellet
  // is moving with the matter.

  // if no gamma spectra is known, then covert straight to kpkts (e.g., Fe52, Mn52)
  if (gamma_spectra[pkt_ptr->pellet_nucindex].nlines == 0) {
    pkt_ptr->type = TYPE_KPKT;
    pkt_ptr->absorptiontype = -6;
    return;
  }

  // Now let's give the gamma ray a direction.
  // Assuming isotropic emission in cmf

  double dir_cmf[3];
  get_rand_isotropic_unitvec(dir_cmf);

  // This direction is in the cmf - we want to convert it to the rest
  // frame - use aberation of angles. We want to convert from cmf to
  // rest so need -ve velocity.

  double vel_vec[3];
  get_velocity(pkt_ptr->pos, vel_vec, -1. * pkt_ptr->tdecay);
  // negative time since we want the backwards transformation here

  angle_ab(dir_cmf, vel_vec, pkt_ptr->dir);

  // Now need to assign the frequency of the packet in the co-moving frame.

  choose_gamma_ray(pkt_ptr);

  // Finally we want to put in the rest frame energy and frequency. And record
  // that it's now a gamma ray.

  pkt_ptr->prop_time = pkt_ptr->tdecay;
  const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr);
  pkt_ptr->nu_rf = pkt_ptr->nu_cmf / dopplerfactor;
  pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;

  pkt_ptr->type = TYPE_GAMMA;
  pkt_ptr->last_cross = BOUNDARY_NONE;

  // initialise polarisation information
  pkt_ptr->stokes[0] = 1.0;
  pkt_ptr->stokes[1] = pkt_ptr->stokes[2] = 0.0;
  double dummy_dir[3];

  dummy_dir[0] = dummy_dir[1] = 0.0;
  dummy_dir[2] = 1.0;
  cross_prod(pkt_ptr->dir, dummy_dir, pkt_ptr->pol_dir);
  if ((dot(pkt_ptr->pol_dir, pkt_ptr->pol_dir)) < 1.e-8) {
    dummy_dir[0] = dummy_dir[2] = 0.0;
    dummy_dir[1] = 1.0;
    cross_prod(pkt_ptr->dir, dummy_dir, pkt_ptr->pol_dir);
  }

  vec_norm(pkt_ptr->pol_dir, pkt_ptr->pol_dir);
  // printout("initialise pol state of packet %g, %g, %g, %g,
  // %g\n",pkt_ptr->stokes_qu[0],pkt_ptr->stokes_qu[1],pkt_ptr->pol_dir[0],pkt_ptr->pol_dir[1],pkt_ptr->pol_dir[2]);
  // printout("pkt direction %g, %g, %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
}

static constexpr auto sigma_compton_partial(const double x, const double f) -> double
// Routine to compute the partial cross section for Compton scattering.
//   xx is the photon energy (in units of electron mass) and f
//  is the energy loss factor up to which we wish to integrate.
{
  const double term1 = ((x * x) - (2 * x) - 2) * std::log(f) / x / x;
  const double term2 = (((f * f) - 1) / (f * f)) / 2;
  const double term3 = ((f - 1) / x) * ((1 / x) + (2 / f) + (1 / (x * f)));

  return (3 * SIGMA_T * (term1 + term2 + term3) / (8 * x));
}

static auto get_kappa_compton_rf(const struct packet *pkt_ptr) -> double {
  // calculate the absorption coefficient [cm^-1] for Compton scattering in the observer reference frame
  // Start by working out the compton x-section in the co-moving frame.

  double const xx = H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT;

  // Use this to decide whether the Thompson limit is acceptable.

  double sigma_cmf;
  if (xx < THOMSON_LIMIT) {
    sigma_cmf = SIGMA_T;
  } else {
    double const fmax = (1 + (2 * xx));
    sigma_cmf = sigma_compton_partial(xx, fmax);
  }

  // Now need to multiply by the electron number density.
  const double kappa_cmf = sigma_cmf * grid::get_nnetot(grid::get_cell_modelgridindex(pkt_ptr->where));

  // convert between frames
  const double kappa_rf = kappa_cmf * doppler_packet_nucmf_on_nurf(pkt_ptr);

  assert_testmodeonly(std::isfinite(kappa_rf));

  return kappa_rf;
}

static auto choose_f(const double xx, const double zrand) -> double
// To choose the value of f to integrate to - idea is we want
//   sigma_compton_partial(xx,f) = zrand.
{
  double f_max = 1 + (2 * xx);
  double f_min = 1;

  const double norm = zrand * sigma_compton_partial(xx, f_max);

  int count = 0;
  double err = 1e20;

  // printout("new\n");

  double ftry = (f_max + f_min) / 2;
  while ((err > 1.e-4) && (count < 1000)) {
    ftry = (f_max + f_min) / 2;
    const double sigma_try = sigma_compton_partial(xx, ftry);
    // printout("ftry %g %g %g %g %g\n",ftry, f_min, f_max, try, norm);
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

static auto thomson_angle() -> double {
  // For Thomson scattering we can get the new angle from a random number very easily.

  const double zrand = rng_uniform();

  const double B_coeff = (8. * zrand) - 4.;

  double t_coeff = sqrt((B_coeff * B_coeff) + 4);
  t_coeff = t_coeff - B_coeff;
  t_coeff = t_coeff / 2;
  t_coeff = cbrt(t_coeff);

  const double mu = (1 / t_coeff) - t_coeff;

  if (fabs(mu) > 1) {
    printout("Error in Thomson. Abort.\n");
    abort();
  }

  return mu;
}

static void compton_scatter(struct packet *pkt_ptr)
// Routine to deal with physical Compton scattering event.
{
  double f = NAN;

  //  printout("Compton scattering.\n");

  const double xx = H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT;

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

  bool stay_gamma = false;
  if (xx < THOMSON_LIMIT) {
    f = 1.0;  // no energy loss
    stay_gamma = true;
  } else {
    const double zrand = rng_uniform();
    f = choose_f(xx, zrand);

    // Check that f lies between 1.0 and (2xx  + 1)

    if ((f < 1) || (f > (2 * xx + 1))) {
      printout("Compton f out of bounds. Abort.\n");
      abort();
    }

    // Prob of keeping gamma ray is...

    const double prob_gamma = 1. / f;

    const double zrand2 = rng_uniform();
    stay_gamma = (zrand2 < prob_gamma);
  }

  if (stay_gamma) {
    // It stays as a gamma ray. Change frequency and direction in
    // co-moving frame then transfer back to rest frame.

    pkt_ptr->nu_cmf = pkt_ptr->nu_cmf / f;  // reduce frequency

    // The packet has stored the direction in the rest frame.
    // Use aberation of angles to get this into the co-moving frame.

    double vel_vec[3];
    get_velocity(pkt_ptr->pos, vel_vec, pkt_ptr->prop_time);

    double cmf_dir[3];
    angle_ab(pkt_ptr->dir, vel_vec, cmf_dir);

    // Now change the direction through the scattering angle.

    const double cos_theta = (xx < THOMSON_LIMIT) ? thomson_angle() : 1. - ((f - 1) / xx);

    double new_dir[3] = {NAN, NAN, NAN};
    scatter_dir(cmf_dir, cos_theta, new_dir);

    const double test = dot(new_dir, new_dir);
    if (fabs(1. - test) > 1.e-8) {
      printout("Not a unit vector - Compton. Abort. %g %g %g\n", f, xx, test);
      printout("new_dir %g %g %g\n", new_dir[0], new_dir[1], new_dir[2]);
      printout("cmf_dir %g %g %g\n", cmf_dir[0], cmf_dir[1], cmf_dir[2]);
      printout("cos_theta %g", cos_theta);
      abort();
    }

    const double test2 = dot(new_dir, cmf_dir);
    if (fabs(test2 - cos_theta) > 1.e-8) {
      printout("Problem with angle - Compton. Abort.\n");
      abort();
    }

    // Now convert back again.

    // get_velocity(pkt_ptr->pos, vel_vec, (-1 * pkt_ptr->prop_time));
    vec_scale(vel_vec, -1.);

    double final_dir[3];
    angle_ab(new_dir, vel_vec, final_dir);

    vec_copy(pkt_ptr->dir, final_dir);

    assert_testmodeonly(std::fabs(vec_len(pkt_ptr->dir) - 1.) < 1e-10);

    // It now has a rest frame direction and a co-moving frequency.
    //  Just need to set the rest frame energy.
    const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr);
    pkt_ptr->nu_rf = pkt_ptr->nu_cmf / dopplerfactor;
    pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;

    pkt_ptr->last_cross = BOUNDARY_NONE;  // allow it to re-cross a boundary
  } else {
    // It's converted to an e-minus packet.
    pkt_ptr->type = TYPE_NTLEPTON;
    pkt_ptr->absorptiontype = -3;
    stats::increment(stats::COUNTER_NT_STAT_FROM_GAMMA);
  }
}

static auto get_kappa_photo_electric_rf(const struct packet *pkt_ptr) -> double {
  // calculate the absorption coefficient [cm^-1] for photo electric effect scattering in the observer reference frame

  double kappa_cmf;
  // Start by working out the x-section in the co-moving frame.

  const int mgi = grid::get_cell_modelgridindex(pkt_ptr->where);
  const double rho = grid::get_rho(mgi);

  if (globals::gamma_grey < 0) {
    // double sigma_cmf_cno = 0.0448e-24 * pow(pkt_ptr->nu_cmf / 2.41326e19, -3.2);

    double sigma_cmf_si = 1.16e-24 * pow(pkt_ptr->nu_cmf / 2.41326e19, -3.13);

    double sigma_cmf_fe = 25.7e-24 * pow(pkt_ptr->nu_cmf / 2.41326e19, -3.0);

    // 2.41326e19 = 100keV in frequency.

    // Now need to multiply by the particle number density.

    // sigma_cmf_cno *= rho * (1. - f_fe) / MH / 14;
    //  Assumes Z = 7. So mass = 14.

    const double kappa_cmf_si = sigma_cmf_si * rho / MH / 28;
    // Assumes Z = 14. So mass = 28.

    const double kappa_cmf_fe = sigma_cmf_fe *= rho / MH / 56;
    // Assumes Z = 28. So mass = 56.

    const double f_fe = grid::get_ffegrp(mgi);

    kappa_cmf = (kappa_cmf_fe * f_fe) + (kappa_cmf_si * (1. - f_fe));
  } else {
    kappa_cmf = globals::gamma_grey * rho;
  }

  // Now need to convert between frames.

  const double kappa_rf = kappa_cmf * doppler_packet_nucmf_on_nurf(pkt_ptr);
  return kappa_rf;
}

static auto get_kappa_pair_prod_rf(const struct packet *pkt_ptr) -> double {
  // calculate the absorption coefficient [cm^-1] for pair production in the observer reference frame

  // Start by working out the x-section in the co-moving frame.

  const int mgi = grid::get_cell_modelgridindex(pkt_ptr->where);
  const double rho = grid::get_rho(mgi);

  if (globals::gamma_grey >= 0.) {
    return 0.;
  }

  // 2.46636e+20 = 1022 keV in frequency
  if (pkt_ptr->nu_cmf <= 2.46636e+20) {
    return 0.;
  }

  // double sigma_cmf_cno;
  double sigma_cmf_si;
  double sigma_cmf_fe;
  const double f_fe = grid::get_ffegrp(mgi);

  // 3.61990e+20 = 1500 keV in frequency
  if (pkt_ptr->nu_cmf > 3.61990e+20) {
    // sigma_cmf_cno = (0.0481 + (0.301 * ((pkt_ptr->nu_cmf/2.41326e+20) - 1.5))) * 49.e-27;

    sigma_cmf_si = (0.0481 + (0.301 * ((pkt_ptr->nu_cmf / 2.41326e+20) - 1.5))) * 196.e-27;

    sigma_cmf_fe = (0.0481 + (0.301 * ((pkt_ptr->nu_cmf / 2.41326e+20) - 1.5))) * 784.e-27;
  } else {
    // sigma_cmf_cno = 1.0063 * ((pkt_ptr->nu_cmf/2.41326e+20) - 1.022) * 49.e-27;

    sigma_cmf_si = 1.0063 * ((pkt_ptr->nu_cmf / 2.41326e+20) - 1.022) * 196.e-27;

    sigma_cmf_fe = 1.0063 * ((pkt_ptr->nu_cmf / 2.41326e+20) - 1.022) * 784.e-27;
  }

  // Now need to multiply by the particle number density.

  // sigma_cmf_cno *= rho * (1. - f_fe) / MH / 14;
  // Assumes Z = 7. So mass = 14.

  const double kappa_cmf_si = sigma_cmf_si * rho / MH / 28;
  // Assumes Z = 14. So mass = 28.

  const double kappa_cmf_fe = sigma_cmf_fe * rho / MH / 56;
  // Assumes Z = 28. So mass = 56.

  const double kappa_cmf = (kappa_cmf_fe * f_fe) + (kappa_cmf_si * (1. - f_fe));

  // Now need to convert between frames.

  const double kappa_rf = kappa_cmf * doppler_packet_nucmf_on_nurf(pkt_ptr);

  assert_testmodeonly(kappa_rf >= 0.);

  return kappa_rf;
}

constexpr auto meanf_sigma(const double x) -> double
// Routine to compute the mean energy converted to non-thermal electrons times
// the Klein-Nishina cross section.
{
  double const f = 1 + (2 * x);

  double const term0 = 2 / x;
  double const term1 = (1 - (2 / x) - (3 / (x * x))) * log(f);
  double const term2 = ((4 / x) + (3 / (x * x)) - 1) * 2 * x / f;
  double const term3 = (1 - (2 / x) - (1 / (x * x))) * 2 * x * (1 + x) / f / f;
  double const term4 = -2. * x * ((4 * x * x) + (6 * x) + 3) / 3 / f / f / f;

  double const tot = 3 * SIGMA_T * (term0 + term1 + term2 + term3 + term4) / (8 * x);

  return tot;
}

static void rlc_emiss_gamma(const struct packet *pkt_ptr, const double dist) {
  // Subroutine to record the heating rate in a cell due to gamma rays.
  // By heating rate I mean, for now, really the rate at which the code is making
  // k-packets in that cell which will then convert into r-packets. This is (going
  // to be) used for the new light_curve syn-style calculation.

  // The intention is that rpkt_emiss will contain the emissivity of r-packets
  // in the co-moving frame (which is going to be isotropic).

  // This is only done to order v/c for now.

  // Called with a packet that is about to travel a
  // distance dist in the lab frame.

  if (!(dist > 0)) {
    return;
  }

  double vel_vec[3];
  get_velocity(pkt_ptr->pos, vel_vec, pkt_ptr->prop_time);

  double const doppler_sq = doppler_squared_nucmf_on_nurf(pkt_ptr->dir, vel_vec);

  const int mgi = grid::get_cell_modelgridindex(pkt_ptr->where);
  const double xx = H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT;
  double heating_cont = ((meanf_sigma(xx) * grid::get_nnetot(mgi)) + get_kappa_photo_electric_rf(pkt_ptr) +
                         (get_kappa_pair_prod_rf(pkt_ptr) * (1. - (2.46636e+20 / pkt_ptr->nu_cmf))));
  heating_cont = heating_cont * pkt_ptr->e_rf * dist * doppler_sq;

  // The terms in the above are for Compton, photoelectric and pair production. The pair production one
  // assumes that a fraction (1. - (1.022 MeV / nu)) of the gamma's energy is thermalised.
  // The remaining 1.022 MeV is made into gamma rays

  // For normalisation this needs to be
  //  1) divided by volume
  //  2) divided by the length of the time step
  //  3) divided by 4 pi sr
  //  This will all be done later
  assert_testmodeonly(heating_cont >= 0.);
  assert_testmodeonly(isfinite(heating_cont));
  safeadd(globals::rpkt_emiss[mgi], heating_cont);
}

void pair_prod(struct packet *pkt_ptr) {
  // Routine to deal with pair production.

  //  In pair production, the original gamma makes an electron positron pair - kinetic energy equal to
  //  gamma ray energy - 1.022 MeV. We assume that the electron deposits any kinetic energy directly to
  //  the thermal pool. The positron annihilates with an electron locally making a pair of gamma rays
  //  at 0.511 MeV in the local cmf (isotropic). So all the thermal energy goes to the thermal pool
  //  immediately and the remainder goes into gamma-rays at 0.511 MeV.

  const double prob_gamma = 1.022 * MEV / (H * pkt_ptr->nu_cmf);

  if (prob_gamma < 0) {
    printout("prob_gamma < 0. pair_prod. Abort. %g\n", prob_gamma);
    abort();
  }

  const double zrand = rng_uniform();

  if (zrand > prob_gamma) {
    // Convert it to an e-minus packet - actually it could be positron EK too, but this works
    // for consistency with compton_scatter.
    pkt_ptr->type = TYPE_NTLEPTON;
    pkt_ptr->absorptiontype = -5;
    stats::increment(stats::COUNTER_NT_STAT_FROM_GAMMA);
  } else {
    // The energy goes into emission at 511 keV.
    pkt_ptr->nu_cmf = 0.511 * MEV / H;

    // Now let's give the gamma ray a direction.

    double dir_cmf[3];
    get_rand_isotropic_unitvec(dir_cmf);

    // This direction is in the cmf - we want to convert it to the rest
    // frame - use aberation of angles. We want to convert from cmf to
    // rest so need -ve velocity.

    double vel_vec[3];
    get_velocity(pkt_ptr->pos, vel_vec, -1. * pkt_ptr->prop_time);
    // negative time since we want the backwards transformation here

    angle_ab(dir_cmf, vel_vec, pkt_ptr->dir);

    const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr);
    pkt_ptr->nu_rf = pkt_ptr->nu_cmf / dopplerfactor;
    pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;

    pkt_ptr->type = TYPE_GAMMA;
    pkt_ptr->last_cross = BOUNDARY_NONE;
  }
}

void do_gamma(struct packet *pkt_ptr, double t2)
// Now routine for moving a gamma packet. Idea is that we have as input
// a gamma packet with known properties at time t1 and we want to follow it
// until time t2.
{
  // Assign optical depth to next physical event. And start counter of
  // optical depth for this path.
  double zrand = rng_uniform_pos();
  const double tau_next = -1. * log(zrand);
  const double tau_current = 0.0;

  // Start by finding the distance to the crossing of the grid cell
  // boundaries. sdist is the boundary distance and snext is the
  // grid cell into which we pass.

  int snext = 0;
  double sdist = grid::boundary_distance(pkt_ptr->dir, pkt_ptr->pos, pkt_ptr->prop_time, pkt_ptr->where, &snext,
                                         &pkt_ptr->last_cross);

  const double maxsdist = (GRID_TYPE == GRID_CARTESIAN3D)
                              ? globals::rmax * pkt_ptr->prop_time / globals::tmin
                              : 2 * globals::rmax * (pkt_ptr->prop_time + sdist / CLIGHT_PROP) / globals::tmin;
  if (sdist > maxsdist) {
    printout("Unreasonably large sdist (gamma). Abort. %g %g %g\n", globals::rmax, pkt_ptr->prop_time / globals::tmin,
             sdist);
    abort();
  }

  if (sdist < 0) {
    printout("Negative distance (sdist). Abort?\n");
    sdist = 0;
  }

  if (((snext < 0) && (snext != -99)) || (snext >= grid::ngrid)) {
    printout("Heading for inappropriate grid cell. Abort.\n");
    printout("Current cell %d, target cell %d.\n", pkt_ptr->where, snext);
    abort();
  }

  if (sdist > globals::max_path_step) {
    sdist = globals::max_path_step;
    snext = pkt_ptr->where;
  }

  // Now consider the scattering/destruction processes.
  // Compton scattering - need to determine the scattering co-efficient.
  // Routine returns the value in the rest frame.

  const double kap_compton = (globals::gamma_grey < 0) ? get_kappa_compton_rf(pkt_ptr) : 0.0;

  const double kap_photo_electric = get_kappa_photo_electric_rf(pkt_ptr);
  const double kap_pair_prod = get_kappa_pair_prod_rf(pkt_ptr);
  const double kap_tot = kap_compton + kap_photo_electric + kap_pair_prod;

  assert_testmodeonly(std::isfinite(kap_compton));
  assert_testmodeonly(std::isfinite(kap_photo_electric));
  assert_testmodeonly(std::isfinite(kap_pair_prod));

  // So distance before physical event is...

  double const edist = (tau_next - tau_current) / kap_tot;

  if (edist < 0) {
    printout("Negative distance (edist). Abort. \n");
    abort();
  }

  // Find how far it can travel during the time inverval.

  double const tdist = (t2 - pkt_ptr->prop_time) * CLIGHT_PROP;

  if (tdist < 0) {
    printout("Negative distance (tdist). Abort. \n");
    abort();
  }

  // printout("sdist, tdist, edist %g %g %g\n",sdist, tdist, edist);

  if ((sdist < tdist) && (sdist < edist)) {
    pkt_ptr->prop_time += sdist / 2. / CLIGHT_PROP;
    move_pkt(pkt_ptr, sdist / 2.);

    // Move it into the new cell.
    if (kap_tot > 0) {
      rlc_emiss_gamma(pkt_ptr, sdist);
    }

    pkt_ptr->prop_time += sdist / 2. / CLIGHT_PROP;
    move_pkt(pkt_ptr, sdist / 2.);

    if (snext != pkt_ptr->where) {
      grid::change_cell(pkt_ptr, snext);
    }
  } else if ((tdist < sdist) && (tdist < edist)) {
    // Doesn't reach boundary.
    pkt_ptr->prop_time += tdist / 2. / CLIGHT_PROP;
    move_pkt(pkt_ptr, tdist / 2.);

    if (kap_tot > 0) {
      rlc_emiss_gamma(pkt_ptr, tdist);
    }
    pkt_ptr->prop_time = t2;
    move_pkt(pkt_ptr, tdist / 2.);
  } else if ((edist < sdist) && (edist < tdist)) {
    pkt_ptr->prop_time += edist / 2. / CLIGHT_PROP;
    move_pkt(pkt_ptr, edist / 2.);
    if (kap_tot > 0) {
      rlc_emiss_gamma(pkt_ptr, edist);
    }
    pkt_ptr->prop_time += edist / 2. / CLIGHT_PROP;
    move_pkt(pkt_ptr, edist / 2.);

    // event occurs. Choose which event and call the appropriate subroutine.
    zrand = rng_uniform();
    if (kap_compton > (zrand * kap_tot)) {
      // Compton scattering.
      compton_scatter(pkt_ptr);
    } else if ((kap_compton + kap_photo_electric) > (zrand * kap_tot)) {
      // Photo electric effect - makes it a k-packet for sure.
      pkt_ptr->type = TYPE_NTLEPTON;
      pkt_ptr->absorptiontype = -4;
      // pkt_ptr->type = TYPE_PRE_KPKT;
      stats::increment(stats::COUNTER_NT_STAT_FROM_GAMMA);
    } else if ((kap_compton + kap_photo_electric + kap_pair_prod) > (zrand * kap_tot)) {
      // It's a pair production
      pair_prod(pkt_ptr);
    } else {
      printout("Failed to identify event. Gamma (1). kap_compton %g kap_photo_electric %g kap_tot %g zrand %g Abort.\n",
               kap_compton, kap_photo_electric, kap_tot, zrand);
      const int cellindex = pkt_ptr->where;
      printout(
          " globals::cell[pkt_ptr->where].rho %g pkt_ptr->nu_cmf %g pkt_ptr->dir[0] %g pkt_ptr->dir[1] %g "
          "pkt_ptr->dir[2] %g pkt_ptr->pos[0] %g pkt_ptr->pos[1] %g pkt_ptr->pos[2] %g \n",
          grid::get_rho(grid::get_cell_modelgridindex(cellindex)), pkt_ptr->nu_cmf, pkt_ptr->dir[0], pkt_ptr->dir[0],
          pkt_ptr->dir[1], pkt_ptr->dir[2], pkt_ptr->pos[1], pkt_ptr->pos[2]);

      abort();
    }
  } else {
    printout("Failed to identify event. Gamma (2). edist %g, sdist %g, tdist %g Abort.\n", edist, sdist, tdist);
    abort();
  }
}

}  // namespace gammapkt