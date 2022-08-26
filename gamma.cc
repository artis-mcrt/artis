#include "gamma.h"

#include <cstring>
#include <fstream>

#include "boundary.h"
#include "decay.h"
#include "emissivities.h"
#include "grey_emissivities.h"
#include "grid.h"
#include "nonthermal.h"
#include "packet.h"
#include "photo_electric.h"
#include "sn3d.h"
#include "stats.h"
#include "vectors.h"

namespace gamma {
// Code for handing gamma rays - creation and propagation

struct gamma_spec {
  double *energy;
  double *probability;
  int nlines;
};

static struct gamma_spec *gamma_spectra;

static const int RED_OF_LIST = -956;  // must be negative

struct gamma_ll {
  int *nucindex;  // is it a Ni56, Co56, a fake line, etc
  int *index;     // which of the lines of that element is it
  int total;      // the total number of lines in the list
};

static struct gamma_ll gam_line_list;

static void read_gamma_spectrum(const int z, const int a, const char filename[50])
// reads in gamma_spectra and returns the average energy in gamma rays per nuclear decay
{
  const int nucindex = decay::get_nuc_index(z, a);

  FILE *filein = fopen_required(filename, "r");
  int nlines = 0;
  assert_always(fscanf(filein, "%d", &nlines) == 1);

  gamma_spectra[nucindex].nlines = nlines;

  gamma_spectra[nucindex].energy = (double *)malloc(nlines * sizeof(double));
  gamma_spectra[nucindex].probability = (double *)malloc(nlines * sizeof(double));

  double E_gamma_avg = 0.;
  for (int n = 0; n < nlines; n++) {
    double en_mev;
    double prob;
    assert_always(fscanf(filein, "%lg %lg", &en_mev, &prob) == 2);
    gamma_spectra[nucindex].energy[n] = en_mev * MEV;
    gamma_spectra[nucindex].probability[n] = prob;
    E_gamma_avg += en_mev * MEV * prob;
  }
  fclose(filein);

  decay::set_nucdecayenergygamma(z, a, E_gamma_avg);

  printout("gamma spectrum for Z=%d A=%d read from %s: nlines %d avg_en_gamma %g MeV\n", z, a, filename, nlines,
           E_gamma_avg / MEV);
}

static void set_trivial_gamma_spectrum(const int z, const int a) {
  // printout("Setting trivial gamma spectrum for z %d a %d engamma %g\n", z, a, decay::nucdecayenergygamma(z, a));
  const int nucindex = decay::get_nuc_index(z, a);
  const int nlines = 1;
  gamma_spectra[nucindex].nlines = nlines;
  gamma_spectra[nucindex].energy = (double *)malloc(nlines * sizeof(double));
  gamma_spectra[nucindex].probability = (double *)malloc(nlines * sizeof(double));
  gamma_spectra[nucindex].energy[0] = decay::nucdecayenergygamma(z, a);
  gamma_spectra[nucindex].probability[0] = 1.;
}

static void read_decaydata(void) {
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

  gamma_spectra = (struct gamma_spec *)calloc(decay::get_num_nuclides(), sizeof(struct gamma_spec));

  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    gamma_spectra[nucindex].nlines = 0;
    gamma_spectra[nucindex].energy = NULL;
    gamma_spectra[nucindex].probability = NULL;
    const int z = decay::get_nuc_z(nucindex);
    const int a = decay::get_nuc_a(nucindex);
    if (z < 1) {
      continue;
    }

    const char *elname = decay::get_elname(z);
    char elnamelower[sizeof(elname)] = "";
    for (int i = 0; i < (int)sizeof(elname); i++) {
      elnamelower[i] = tolower(elname[i]);
    }

    // look in the current folder
    char filename[128];
    snprintf(filename, 128, "%s%d_lines.txt", elnamelower, a);

    // look in the 'data' subfolder
    char filename2[128];
    snprintf(filename2, 128, "data/%s%d_lines.txt", elnamelower, a);

    if (std::ifstream(filename)) {
      read_gamma_spectrum(z, a, filename);
    } else if (std::ifstream(filename2)) {
      read_gamma_spectrum(z, a, filename2);
    } else if (decay::nucdecayenergygamma(z, a) > 0.) {
      // printout("%s does not exist. Setting 100%% chance of single gamma-line with energy %g MeV\n",
      //   filename, decay::nucdecayenergygamma(z, a) / EV / 1e6);
      set_trivial_gamma_spectrum(z, a);
    } else {
      // printout("%s does not exist. No gamma decay from this nuclide.\n", filename);
    }
  }

  // read_gamma_spectrum(28, 56, "ni56_lines.txt");
  //
  // read_gamma_spectrum(27, 56, "co56_lines.txt");
  //
  // read_gamma_spectrum(23, 48, "v48_lines.txt");
  //
  // read_gamma_spectrum(24, 48, "cr48_lines.txt");
  //
  // read_gamma_spectrum(28, 57, "ni57_lines.txt");
  //
  // read_gamma_spectrum(27, 57, "co57_lines.txt");
  if (decay::nuc_exists(26, 52)) {
    decay::set_nucdecayenergygamma(26, 52, 0.86 * MEV);  // Fe52
  }
  if (decay::nuc_exists(25, 52)) {
    decay::set_nucdecayenergygamma(25, 52, 3.415 * MEV);  // Mn52
  }
}

// construct an energy ordered gamma ray line list.
void init_gamma_linelist(void) {
  read_decaydata();

  /* Start by setting up the grid of fake lines and their energies. */
  gamma_spectra[FAKE_GAM_LINE_ID].nlines = globals::nfake_gam;
  gamma_spectra[FAKE_GAM_LINE_ID].energy = (double *)malloc(globals::nfake_gam * sizeof(double));
  gamma_spectra[FAKE_GAM_LINE_ID].probability = (double *)malloc(globals::nfake_gam * sizeof(double));

  const double deltanu = (globals::nusyn_max - globals::nusyn_min) / (gamma_spectra[FAKE_GAM_LINE_ID].nlines - 3);
  for (int i = 0; i < gamma_spectra[FAKE_GAM_LINE_ID].nlines; i++) {
    gamma_spectra[FAKE_GAM_LINE_ID].energy[i] = (globals::nusyn_min + deltanu * (i - 1)) * H;
    gamma_spectra[FAKE_GAM_LINE_ID].probability[i] = 0.0;
  }

  /* Now do the sorting. */

  int total_lines = 0;
  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    total_lines += gamma_spectra[nucindex].nlines;
  }
  printout("total gamma-ray lines %d\n", total_lines);

  gam_line_list.total = total_lines;
  gam_line_list.nucindex = (int *)malloc(total_lines * sizeof(int));
  gam_line_list.index = (int *)malloc(total_lines * sizeof(int));

  double energy_last = 0.0;
  int next = -99;
  int next_type;

  for (int i = 0; i < total_lines; i++) {
    double energy_try = 1.e50;

    for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
      // printout("nucindex %d nlines %d\n", nucindex, gamma_spectra[nucindex].nlines);
      for (int j = 0; j < gamma_spectra[nucindex].nlines; j++) {
        if (gamma_spectra[nucindex].energy[j] > energy_last && gamma_spectra[nucindex].energy[j] < energy_try) {
          // next_type = spec_type[iso];
          next_type = nucindex;
          next = j;
          energy_try = gamma_spectra[nucindex].energy[j];
        }
      }
    }

    assert_always(next >= 0);
    gam_line_list.nucindex[i] = next_type;
    gam_line_list.index[i] = next;
    energy_last = energy_try;
  }

  FILE *const line_list = fopen_required("gammalinelist.out", "w+");

  fprintf(line_list, "#index nucindex Z A nucgammmaindex en_gamma_mev gammaline_probability\n");
  for (int i = 0; i < total_lines; i++) {
    const int nucindex = gam_line_list.nucindex[i];
    const int index = gam_line_list.index[i];
    fprintf(line_list, "%d %d %d %d %d %g %g \n", i, gam_line_list.nucindex[i],
            decay::get_nuc_z(gam_line_list.nucindex[i]), decay::get_nuc_a(gam_line_list.nucindex[i]),
            gam_line_list.index[i], gamma_spectra[nucindex].energy[index] / MEV,
            gamma_spectra[nucindex].probability[index]);
  }
  fclose(line_list);
}

static void choose_gamma_ray(struct packet *pkt_ptr) {
  // Routine to choose which gamma ray line it'll be.

  const int nucindex = pkt_ptr->pellet_nucindex;
  const int z = decay::get_nuc_z(nucindex);
  const int a = decay::get_nuc_a(nucindex);
  double E_gamma = decay::nucdecayenergygamma(z, a);  // Average energy per gamma line of a decay

  const double zrand = gsl_rng_uniform(rng);
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

void pellet_gamma_decay(const int nts, struct packet *pkt_ptr) {
  // Subroutine to convert a pellet to a gamma ray (or kpkt if no gamma spec loaded)

  // nts defines the time step we are in. pkt_ptr is a pointer to the packet
  // that is decaying.

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
  pkt_ptr->last_cross = NONE;

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

static double sigma_compton_partial(const double x, const double f)
// Routine to compute the partial cross section for Compton scattering.
//   xx is the photon energy (in units of electron mass) and f
//  is the energy loss factor up to which we wish to integrate.
{
  const double term1 = ((x * x) - (2 * x) - 2) * log(f) / x / x;
  const double term2 = (((f * f) - 1) / (f * f)) / 2;
  const double term3 = ((f - 1) / x) * ((1 / x) + (2 / f) + (1 / (x * f)));

  return (3 * SIGMA_T * (term1 + term2 + term3) / (8 * x));
}

static double sig_comp(const struct packet *pkt_ptr) {
  // Start by working out the compton x-section in the co-moving frame.

  double xx = H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT;

  // Use this to decide whether the Thompson limit is acceptable.

  double sigma_cmf;
  if (xx < THOMSON_LIMIT) {
    sigma_cmf = SIGMA_T;
  } else {
    double fmax = (1 + (2 * xx));
    sigma_cmf = sigma_compton_partial(xx, fmax);
  }

  // Now need to multiply by the electron number density.
  const int cellindex = pkt_ptr->where;
  sigma_cmf *= grid::get_nnetot(grid::get_cell_modelgridindex(cellindex));

  // Now need to convert between frames.

  const double sigma_rf = sigma_cmf * doppler_packet_nucmf_on_nurf(pkt_ptr);

  assert_testmodeonly(std::isfinite(sigma_rf));

  return sigma_rf;
}

static double choose_f(double xx, double zrand)
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

static double thomson_angle(void) {
  // For Thomson scattering we can get the new angle from a random number very easily.

  const double zrand = gsl_rng_uniform(rng);

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
  double f;

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

  bool stay_gamma;
  if (xx < THOMSON_LIMIT) {
    f = 1.0;  // no energy loss
    stay_gamma = true;
  } else {
    const double zrand = gsl_rng_uniform(rng);
    f = choose_f(xx, zrand);

    // Check that f lies between 1.0 and (2xx  + 1)

    if ((f < 1) || (f > (2 * xx + 1))) {
      printout("Compton f out of bounds. Abort.\n");
      abort();
    }

    // Prob of keeping gamma ray is...

    const double prob_gamma = 1. / f;

    const double zrand2 = gsl_rng_uniform(rng);
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

    double cos_theta;
    if (xx < THOMSON_LIMIT) {
      cos_theta = thomson_angle();
    } else {
      cos_theta = 1. - ((f - 1) / xx);
    }

    double new_dir[3];
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

    // It now has a rest frame direction and a co-moving frequency.
    //  Just need to set the rest frame energy.
    const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr);
    pkt_ptr->nu_rf = pkt_ptr->nu_cmf / dopplerfactor;
    pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;

    pkt_ptr->last_cross = NONE;  // allow it to re-cross a boundary
  } else {
    // It's converted to an e-minus packet.
    pkt_ptr->type = TYPE_NTLEPTON;
    pkt_ptr->absorptiontype = -3;
    stats::increment(stats::COUNTER_NT_STAT_FROM_GAMMA);
  }
}

void do_gamma(struct packet *pkt_ptr, double t2)
// Now routine for moving a gamma packet. Idea is that we have as input
// a gamma packet with known properties at time t1 and we want to follow it
// until time t2.
{
  // Assign optical depth to next physical event. And start counter of
  // optical depth for this path.
  double zrand = gsl_rng_uniform_pos(rng);
  const double tau_next = -1. * log(zrand);
  const double tau_current = 0.0;

  // Start by finding the distance to the crossing of the grid cell
  // boundaries. sdist is the boundary distance and snext is the
  // grid cell into which we pass.

  int snext;
  double sdist = boundary_cross(pkt_ptr, pkt_ptr->prop_time, &snext);

  const double maxsdist = (grid::grid_type == GRID_SPHERICAL1D)
                              ? 2 * globals::rmax * (pkt_ptr->prop_time + sdist / globals::CLIGHT_PROP) / globals::tmin
                              : globals::rmax * pkt_ptr->prop_time / globals::tmin;
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

  /* Now consider the scattering/destruction processes. */
  /* Compton scattering - need to determine the scattering co-efficient.*/
  /* Routine returns the value in the rest frame. */

  double kap_compton = 0.0;
  if (globals::gamma_grey < 0) {
    kap_compton = sig_comp(pkt_ptr);
  }

  const double kap_photo_electric = sig_photo_electric(pkt_ptr);
  const double kap_pair_prod = sig_pair_prod(pkt_ptr);
  const double kap_tot = kap_compton + kap_photo_electric + kap_pair_prod;

  assert_testmodeonly(std::isfinite(kap_compton));
  assert_testmodeonly(std::isfinite(kap_photo_electric));
  assert_testmodeonly(std::isfinite(kap_pair_prod));

  // So distance before physical event is...

  double edist = (tau_next - tau_current) / kap_tot;

  if (edist < 0) {
    printout("Negative distance (edist). Abort. \n");
    abort();
  }

  // Find how far it can travel during the time inverval.

  double tdist = (t2 - pkt_ptr->prop_time) * globals::CLIGHT_PROP;

  if (tdist < 0) {
    printout("Negative distance (tdist). Abort. \n");
    abort();
  }

  // printout("sdist, tdist, edist %g %g %g\n",sdist, tdist, edist);

  if ((sdist < tdist) && (sdist < edist)) {
    pkt_ptr->prop_time += sdist / 2. / globals::CLIGHT_PROP;
    move_pkt(pkt_ptr, sdist / 2., pkt_ptr->prop_time);

    // Move it into the new cell.
    if (kap_tot > 0) {
      if (globals::do_comp_est) {
        compton_emiss_cont(pkt_ptr, sdist);
        pp_emiss_cont(pkt_ptr, sdist);
      }
      if (globals::do_rlc_est != 0) {
        rlc_emiss_gamma(pkt_ptr, sdist);
      }
    }

    pkt_ptr->prop_time += sdist / 2. / globals::CLIGHT_PROP;
    move_pkt(pkt_ptr, sdist / 2., pkt_ptr->prop_time);

    if (snext != pkt_ptr->where) {
      change_cell(pkt_ptr, snext, pkt_ptr->prop_time);
    }
  } else if ((tdist < sdist) && (tdist < edist)) {
    // Doesn't reach boundary.
    pkt_ptr->prop_time += tdist / 2. / globals::CLIGHT_PROP;
    move_pkt(pkt_ptr, tdist / 2., pkt_ptr->prop_time);

    if (kap_tot > 0) {
      if (globals::do_comp_est) {
        compton_emiss_cont(pkt_ptr, tdist);
        pp_emiss_cont(pkt_ptr, tdist);
      }
      if (globals::do_rlc_est != 0) {
        rlc_emiss_gamma(pkt_ptr, tdist);
      }
    }
    pkt_ptr->prop_time = t2;
    move_pkt(pkt_ptr, tdist / 2., pkt_ptr->prop_time);
  } else if ((edist < sdist) && (edist < tdist)) {
    pkt_ptr->prop_time += edist / 2. / globals::CLIGHT_PROP;
    move_pkt(pkt_ptr, edist / 2., pkt_ptr->prop_time);
    if (kap_tot > 0) {
      if (globals::do_comp_est) {
        compton_emiss_cont(pkt_ptr, edist);
        pp_emiss_cont(pkt_ptr, edist);
      }
      if (globals::do_rlc_est != 0) {
        rlc_emiss_gamma(pkt_ptr, edist);
      }
    }
    pkt_ptr->prop_time += edist / 2. / globals::CLIGHT_PROP;
    move_pkt(pkt_ptr, edist / 2., pkt_ptr->prop_time);

    // event occurs. Choose which event and call the appropriate subroutine.
    zrand = gsl_rng_uniform(rng);
    if (kap_compton > (zrand * kap_tot)) {
      // Compton scattering.
      compton_scatter(pkt_ptr);
    } else if ((kap_compton + kap_photo_electric) > (zrand * kap_tot)) {
      // Photo electric effect - makes it a k-packet for sure.
      pkt_ptr->type = TYPE_NTLEPTON;
      pkt_ptr->absorptiontype = -4;
#ifndef FORCE_LTE
      // kgammadep[pkt_ptr->where] += pkt_ptr->e_cmf;
#endif
      // pkt_ptr->type = TYPE_PRE_KPKT;
      // pkt_ptr->type = TYPE_GAMMA_KPKT;
      // if (tid == 0) nt_stat_from_gamma++;
      stats::increment(stats::COUNTER_NT_STAT_FROM_GAMMA);
    } else if ((kap_compton + kap_photo_electric + kap_pair_prod) > (zrand * kap_tot)) {
      // It's a pair production
      pair_prod(pkt_ptr);
    } else {
      printout("Failed to identify event. Gamma (1). kap_compton %g kap_photo_electric %g kap_tot %g zrand %g Abort.\n",
               kap_compton, kap_photo_electric, kap_tot, zrand);
      const int cellindex = pkt_ptr->where;
      printout(
          " /*globals::cell[*/pkt_ptr->where].rho %g pkt_ptr->nu_cmf %g pkt_ptr->dir[0] %g pkt_ptr->dir[1] %g "
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

double get_gam_freq(const int n) {
  if (n == RED_OF_LIST) {
    return 0.0;
  }

  // returns the frequency of line n
  const int nucindex = gam_line_list.nucindex[n];
  const int lineid = gam_line_list.index[n];

  if (nucindex >= decay::get_num_nuclides() || lineid >= gamma_spectra[nucindex].nlines) {
    printout("Unknown line. %d Abort.\n", n);
    printout("line_list->nucindex[n] %d line_list->index[n] %d\n", gam_line_list.nucindex[n], gam_line_list.index[n]);
    abort();
  }

  return gamma_spectra[nucindex].energy[lineid] / H;
}

int get_nul(double freq) {
  const double freq_max = get_gam_freq(gam_line_list.total - 1);
  const double freq_min = get_gam_freq(0);

  if (freq > freq_max) {
    return (gam_line_list.total - 1);
  } else if (freq < freq_min) {
    return RED_OF_LIST;
  } else {
    int too_high = gam_line_list.total - 1;
    int too_low = 0;

    while (too_high != too_low + 1) {
      const int tryindex = (too_high + too_low) / 2;
      const double freq_try = get_gam_freq(tryindex);
      if (freq_try >= freq) {
        too_high = tryindex;
      } else {
        too_low = tryindex;
      }
    }

    return too_low;
  }
}

}  // namespace gamma