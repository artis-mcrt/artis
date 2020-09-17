#include <assert.h>
#include "sn3d.h"
#include "boundary.h"
#include "emissivities.h"
#include "gamma.h"
#include "grey_emissivities.h"
#include "grid_init.h"
#include "nonthermal.h"
#include "photo_electric.h"
#include "vectors.h"

// Material for handing gamma rays - creation and propagation.

struct gamma_spec
{
  double *energy;
  double *probability;
  int nlines;
};

static struct gamma_spec gamma_spectra[RADIONUCLIDE_COUNT];

static const int RED_OF_LIST = -956;  // must be negative

struct gamma_ll
{
  enum radionuclides *nuclidetype; // is it a Ni56, Co56, a fake line, etc
  int *index;               // which of the lines of that element is it
  int total;                // the total number of lines in the list
};

static struct gamma_ll gam_line_list;


static void read_gamma_spectrum(enum radionuclides isotope, const char filename[50])
// reads in gamma_spectra and returns the average energy in gamma rays per nuclear decay
{
  assert(isotope < RADIONUCLIDE_COUNT);

  FILE *filein = fopen_required(filename, "r");
  int nlines = 0;
  fscanf(filein, "%d", &nlines);

  gamma_spectra[isotope].nlines = nlines;

  gamma_spectra[isotope].energy = (double *) calloc(nlines, sizeof(double));
  gamma_spectra[isotope].probability = (double *) calloc(nlines, sizeof(double));

  double E_gamma_avg = 0.0;
  for (int n = 0; n < nlines; n++)
  {
    double en_mev;
    double prob;
    fscanf(filein, "%lg %lg", &en_mev, &prob);
    gamma_spectra[isotope].energy[n] = en_mev * MEV;
    gamma_spectra[isotope].probability[n] = prob;
    E_gamma_avg += en_mev * MEV * prob;
  }
  fclose(filein);

  set_nucdecayenergygamma(isotope, E_gamma_avg);
}


static void read_decaydata(void)
{
  for (int iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
  {
    gamma_spectra[iso].nlines = 0;
    gamma_spectra[iso].energy = NULL;
    gamma_spectra[iso].probability = NULL;
    set_nucdecayenergygamma((enum radionuclides) iso, 0.);
  }

  read_gamma_spectrum(NUCLIDE_NI56, "ni_lines.txt");

  read_gamma_spectrum(NUCLIDE_CO56, "co_lines.txt");

  read_gamma_spectrum(NUCLIDE_V48, "v48_lines.txt");

  read_gamma_spectrum(NUCLIDE_CR48, "cr48_lines.txt");

  read_gamma_spectrum(NUCLIDE_NI57, "ni57_lines.txt");

  read_gamma_spectrum(NUCLIDE_CO57, "co57_lines.txt");

  set_nucdecayenergygamma(NUCLIDE_FE52, 0.86 * MEV);
  set_nucdecayenergygamma(NUCLIDE_MN52, 3.415 * MEV);
}


// construct an energy ordered gamma ray line list.
void init_gamma_linelist(void)
{
  read_decaydata();

  /* Start by setting up the grid of fake lines and their energies. */
  gamma_spectra[FAKE_GAM_LINE_ID].nlines = nfake_gam;
  gamma_spectra[FAKE_GAM_LINE_ID].energy = (double *) malloc(nfake_gam * sizeof(double));
  gamma_spectra[FAKE_GAM_LINE_ID].probability = (double *) malloc(nfake_gam * sizeof(double));

  const double deltanu = (nusyn_max - nusyn_min) / (gamma_spectra[FAKE_GAM_LINE_ID].nlines - 3);
  for (int i = 0; i < gamma_spectra[FAKE_GAM_LINE_ID].nlines; i++)
  {
    gamma_spectra[FAKE_GAM_LINE_ID].energy[i] = (nusyn_min + deltanu * (i - 1)) * H;
    gamma_spectra[FAKE_GAM_LINE_ID].probability[i] = 0.0;
  }

  /* Now do the sorting. */

  int total_lines = 0;
  for (int iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
  {
    total_lines += gamma_spectra[iso].nlines;
  }
  printout("total gamma-ray lines %d\n", total_lines);

  gam_line_list.total = total_lines;
  gam_line_list.nuclidetype = (enum radionuclides *) malloc(total_lines * sizeof(enum radionuclides));
  gam_line_list.index = (int *) malloc(total_lines * sizeof(int));

  double energy_last = 0.0;
  int next = -99;
  enum radionuclides next_type;

  for (int i = 0; i < total_lines; i++)
  {
    double energy_try = 1.e50;

    for (int isoint = 0; isoint < RADIONUCLIDE_COUNT; isoint++)
    {
      enum radionuclides iso = (enum radionuclides) isoint;
      // printout("iso %d nlines %d\n", iso, gamma_spectra[iso].nlines);
      for (int j = 0; j < gamma_spectra[iso].nlines; j++)
      {
        if (gamma_spectra[iso].energy[j] > energy_last && gamma_spectra[iso].energy[j] < energy_try)
        {
          // next_type = spec_type[iso];
          next_type = iso;
          next = j;
          energy_try = gamma_spectra[iso].energy[j];
        }
      }
    }

    assert(next >= 0);
    gam_line_list.nuclidetype[i] = next_type;
    gam_line_list.index[i] = next;
    energy_last = energy_try;
  }

  FILE *const line_list = fopen_required("gammalinelist.out", "w+");

  for (int i = 0; i < total_lines; i++)
  {
    const enum radionuclides iso = gam_line_list.nuclidetype[i];
    const int index = gam_line_list.index[i];
    fprintf(line_list, "%d %d %d %g %g \n",
            i, gam_line_list.nuclidetype[i], gam_line_list.index[i],
            gamma_spectra[iso].energy[index] / MEV, gamma_spectra[iso].probability[index]);
  }
  fclose(line_list);
}


static void choose_gamma_ray(PKT *pkt_ptr)
{
  // Routine to choose which gamma ray line it'll be.

  enum radionuclides iso;
  switch (pkt_ptr->type)
  {
    case TYPE_56NI_PELLET:
      iso = NUCLIDE_NI56;
      break;

    case TYPE_56CO_PELLET:
      iso = NUCLIDE_CO56;
      break;

    case TYPE_57NI_PELLET:
      iso = NUCLIDE_NI57;
      break;

    case TYPE_57CO_PELLET:
      iso = NUCLIDE_CO57;
      break;

    case TYPE_48CR_PELLET:
      iso = NUCLIDE_CR48;
      break;

    case TYPE_48V_PELLET:
      iso = NUCLIDE_V48;
      break;

    default:
      printout("Unrecognised pellet. Abort.\n");
      abort();
  }

  double E_gamma = nucdecayenergygamma(iso); // Average energy per gamma line of a decay

  const double zrand = gsl_rng_uniform(rng);
  int nselected = -1;
  double runtot = 0.;
  for (int n = 0; n < gamma_spectra[iso].nlines; n++)
  {
    runtot += gamma_spectra[iso].probability[n] * gamma_spectra[iso].energy[n] / E_gamma;
    if (zrand <= runtot)
    {
      nselected = n;
      break;
    }
  }

  if (nselected < 0)
  {
    printout("Failure to choose line (packet type %d). Abort. zrand %g runtot %g\n", pkt_ptr->type, zrand, runtot);
    abort();
  }

  pkt_ptr->nu_cmf = gamma_spectra[iso].energy[nselected] / H;
  // printout("%s PELLET %g\n", gammaspec->filename, gammaspec->energy[nselected]);
}


void pellet_decay(const int nts, PKT *pkt_ptr)
{
  // Subroutine to convert a pellet to a gamma ray.
  // nts defines the time step we are in. pkt_ptr is a pointer to the packet
  // that is decaying.
  // Record decay.

  // Start by getting the position of the pellet at the point of decay. Pellet
  // is moving with the matter.

  // Now let's give the gamma ray a direction.

  // Assuming isotropic emission in cmf

  double dir_cmf[3];
  get_rand_isotropic_unitvec(dir_cmf);

  // This direction is in the cmf - we want to convert it to the rest
  // frame - use aberation of angles. We want to convert from cmf to
  // rest so need -ve velocity.

  double vel_vec[3];
  get_velocity(pkt_ptr->pos, vel_vec, -1. * pkt_ptr->tdecay);
  //negative time since we want the backwards transformation here

  angle_ab(dir_cmf, vel_vec, pkt_ptr->dir);

  // Now need to assign the frequency of the packet in the co-moving frame.

  choose_gamma_ray(pkt_ptr);

  // Finally we want to put in the rest frame energy and frequency. And record
  // that it's now a gamma ray.

  pkt_ptr->prop_time = pkt_ptr->tdecay;
  const double dopplerfactor = doppler_packetpos(pkt_ptr);
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
  if ((dot(pkt_ptr->pol_dir, pkt_ptr->pol_dir)) < 1.e-8)
  {
    dummy_dir[0] = dummy_dir[2] = 0.0;
    dummy_dir[1] = 1.0;
    cross_prod(pkt_ptr->dir, dummy_dir, pkt_ptr->pol_dir);
  }

  vec_norm(pkt_ptr->pol_dir, pkt_ptr->pol_dir);
  //printout("initialise pol state of packet %g, %g, %g, %g, %g\n",pkt_ptr->stokes_qu[0],pkt_ptr->stokes_qu[1],pkt_ptr->pol_dir[0],pkt_ptr->pol_dir[1],pkt_ptr->pol_dir[2]);
  //printout("pkt direction %g, %g, %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
}


static double sigma_compton_partial(const double x, const double f)
// Routine to compute the partial cross section for Compton scattering.
//   xx is the photon energy (in units of electron mass) and f
//  is the energy loss factor up to which we wish to integrate.
{
  const double term1 = ( (x * x) - (2 * x) - 2 ) * log(f) / x / x;
  const double term2 = ( ((f * f) -1) / (f * f)) / 2;
  const double term3 = ( (f - 1) / x) * ( (1 / x) + (2 / f) + (1 / (x * f)));

  return (3 * SIGMA_T * (term1 + term2 + term3) / (8 * x));
}


static double sig_comp(const PKT *pkt_ptr)
{
  // Start by working out the compton x-section in the co-moving frame.

  double xx = H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT;

  // Use this to decide whether the Thompson limit is acceptable.

  double sigma_cmf;
  if (xx < THOMSON_LIMIT)
  {
    sigma_cmf = SIGMA_T;
  }
  else
  {
    double fmax = (1 + (2*xx));
    sigma_cmf = sigma_compton_partial(xx, fmax);
  }

  // Now need to multiply by the electron number density.
  const int cellindex = pkt_ptr->where;
  sigma_cmf *= get_nnetot(cell[cellindex].modelgridindex);

  // Now need to convert between frames.

  const double sigma_rf = sigma_cmf * doppler_packetpos(pkt_ptr);

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

  //printout("new\n");

  double ftry;
  while ((err > 1.e-4) && (count < 1000))
  {
    ftry = (f_max + f_min) / 2;
    const double sigma_try = sigma_compton_partial(xx, ftry);
    //printout("ftry %g %g %g %g %g\n",ftry, f_min, f_max, try, norm);
    if (sigma_try > norm)
    {
      f_max = ftry;
      err = (sigma_try - norm) / norm;
    }
    else
    {
      f_min = ftry;
      err = (norm - sigma_try) / norm;
    }
    //      printout("error %g\n",err);
    count++;
    if (count == 1000)
    {
      printout("Compton hit 1000 tries. %g %g %g %g %g\n", f_max, f_min, ftry, sigma_try, norm);
    }
  }

  return ftry;
}


static double thomson_angle(void)
{
  // For Thomson scattering we can get the new angle from a random number very easily.

  const double zrand = gsl_rng_uniform(rng);

  const double B_coeff = (8. * zrand) - 4.;

  double t_coeff = sqrt((B_coeff * B_coeff) + 4);
  t_coeff = t_coeff - B_coeff;
  t_coeff = t_coeff / 2;
  t_coeff = cbrt(t_coeff);

  const double mu = (1 / t_coeff) - t_coeff;

  if (fabs(mu) > 1)
  {
    printout("Error in Thomson. Abort.\n");
    abort();
  }

  return mu;
}


static void compton_scatter(PKT *pkt_ptr)
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
  if (xx < THOMSON_LIMIT)
  {
    f = 1.0; // no energy loss
    stay_gamma = true;
  }
  else
  {
    const double zrand = gsl_rng_uniform(rng);
    f = choose_f(xx, zrand);

    // Check that f lies between 1.0 and (2xx  + 1)

    if ((f < 1) || (f > (2 * xx + 1)))
    {
      printout("Compton f out of bounds. Abort.\n");
      abort();
    }

    // Prob of keeping gamma ray is...

    const double prob_gamma = 1. / f;

    const double zrand2 = gsl_rng_uniform(rng);
    stay_gamma = (zrand2 < prob_gamma);
  }

  if (stay_gamma)
  {
    // It stays as a gamma ray. Change frequency and direction in
    // co-moving frame then transfer back to rest frame.

    pkt_ptr->nu_cmf = pkt_ptr->nu_cmf / f; // reduce frequency

    // The packet has stored the direction in the rest frame.
    // Use aberation of angles to get this into the co-moving frame.

    double vel_vec[3];
    get_velocity(pkt_ptr->pos, vel_vec, pkt_ptr->prop_time);

    double cmf_dir[3];
    angle_ab(pkt_ptr->dir, vel_vec, cmf_dir);

    // Now change the direction through the scattering angle.

    double cos_theta;
    if (xx <  THOMSON_LIMIT)
    {
      cos_theta = thomson_angle();
    }
    else
    {
      cos_theta = 1. - ((f - 1) / xx);
    }

    double new_dir[3];
    scatter_dir(cmf_dir, cos_theta, new_dir);

    const double test = dot(new_dir, new_dir);
    if (fabs(1. - test) > 1.e-8)
    {
      printout("Not a unit vector - Compton. Abort. %g %g %g\n", f, xx, test);
      printout("new_dir %g %g %g\n", new_dir[0], new_dir[1], new_dir[2]);
      printout("cmf_dir %g %g %g\n", cmf_dir[0], cmf_dir[1], cmf_dir[2]);
      printout("cos_theta %g", cos_theta);
      abort();
    }

    const double test2 = dot(new_dir, cmf_dir);
    if (fabs(test2 - cos_theta) > 1.e-8)
    {
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
    const double dopplerfactor = doppler_packetpos(pkt_ptr);
    pkt_ptr->nu_rf = pkt_ptr->nu_cmf / dopplerfactor;
    pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;

    pkt_ptr->last_cross = NONE; // allow it to re-cross a boundary
  }
  else
  {
    // It's converted to an e-minus packet.
    pkt_ptr->type = TYPE_NTLEPTON;
    pkt_ptr->absorptiontype = -3;
    nt_stat_from_gamma++;
  }
}


void do_gamma(PKT *pkt_ptr, double t2)
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

  const double maxsdist = (grid_type == GRID_SPHERICAL1D) ? 2 * rmax * (pkt_ptr->prop_time + sdist / CLIGHT_PROP) / tmin : rmax * pkt_ptr->prop_time / tmin;
  if (sdist > maxsdist)
  {
    printout("Unreasonably large sdist (gamma). Abort. %g %g %g\n", rmax, pkt_ptr->prop_time/tmin, sdist);
    abort();
  }

  if (sdist < 1)
  {
    printout("Negative distance (sdist). Abort?\n");
    sdist = 0;
  }

  if (((snext < 0) && (snext != -99)) || (snext >= ngrid))
  {
    printout("Heading for inappropriate grid cell. Abort.\n");
    printout("Current cell %d, target cell %d.\n", pkt_ptr->where, snext);
    abort();
  }

  if (sdist > max_path_step)
  {
    sdist = max_path_step;
    snext = pkt_ptr->where;
  }

  /* Now consider the scattering/destruction processes. */
  /* Compton scattering - need to determine the scattering co-efficient.*/
  /* Routine returns the value in the rest frame. */

  double kap_compton = 0.0;
  if (gamma_grey < 0)
  {
    kap_compton = sig_comp(pkt_ptr);
  }

  const double kap_photo_electric = sig_photo_electric(pkt_ptr);
  const double kap_pair_prod = sig_pair_prod(pkt_ptr);
  const double kap_tot = kap_compton + kap_photo_electric + kap_pair_prod;

  // So distance before physical event is...

  double edist = (tau_next - tau_current) / kap_tot;

  if (edist < 0)
  {
    printout("Negative distance (edist). Abort. \n");
    abort();
  }

  // Find how far it can travel during the time inverval.

  double tdist = (t2 - pkt_ptr->prop_time) * CLIGHT_PROP;

  if (tdist < 0)
  {
    printout("Negative distance (tdist). Abort. \n");
    abort();
  }

  //printout("sdist, tdist, edist %g %g %g\n",sdist, tdist, edist);

  if ((sdist < tdist) && (sdist < edist))
  {
    sdist = sdist / 2.;
    pkt_ptr->prop_time += sdist / CLIGHT_PROP;
    move_pkt(pkt_ptr, sdist, pkt_ptr->prop_time);

    // Move it into the new cell.
    if (kap_tot > 0)
    {
      if (do_comp_est)
      {
        sdist = sdist * 2.;
        compton_emiss_cont(pkt_ptr, sdist);
        pp_emiss_cont(pkt_ptr, sdist);
        sdist = sdist / 2.;
      }
      if (do_rlc_est != 0)
      {
        sdist = sdist * 2.;
        rlc_emiss_gamma(pkt_ptr, sdist);
        sdist = sdist / 2.;
      }
    }

    pkt_ptr->prop_time += sdist / CLIGHT_PROP;
    move_pkt(pkt_ptr, sdist, pkt_ptr->prop_time);
    sdist = sdist * 2.;
    if (snext != pkt_ptr->where)
    {
      bool end_packet = false; //tells us when to stop working on this packet
      change_cell(pkt_ptr, snext, &end_packet, pkt_ptr->prop_time);
    }
  }
  else if ((tdist < sdist) && (tdist < edist))
  {
    // Doesn't reach boundary.
    tdist = tdist / 2.;
    pkt_ptr->prop_time += tdist / CLIGHT_PROP;
    move_pkt(pkt_ptr, tdist, pkt_ptr->prop_time);

    if (kap_tot > 0)
    {
      if (do_comp_est)
      {
        tdist = tdist * 2.;
        compton_emiss_cont(pkt_ptr, tdist);
        pp_emiss_cont(pkt_ptr, tdist);
        tdist = tdist / 2.;
      }
      if (do_rlc_est != 0)
      {
        tdist = tdist * 2.;
        rlc_emiss_gamma(pkt_ptr, tdist);
        tdist = tdist / 2.;
      }
    }
    pkt_ptr->prop_time = t2;
    move_pkt(pkt_ptr, tdist, pkt_ptr->prop_time);
    tdist = tdist * 2.;
  }
  else if ((edist < sdist) && (edist < tdist))
  {
    edist = edist / 2.;
    pkt_ptr->prop_time += edist / CLIGHT_PROP;
    move_pkt(pkt_ptr, edist, pkt_ptr->prop_time);
    if (kap_tot > 0)
    {
      if (do_comp_est)
      {
        edist = edist * 2.;
        compton_emiss_cont(pkt_ptr, edist);
        pp_emiss_cont(pkt_ptr, edist);
        edist = edist / 2.;
      }
      if (do_rlc_est != 0)
      {
        edist = edist * 2.;
        rlc_emiss_gamma(pkt_ptr, edist);
        edist = edist / 2.;
      }
    }
    pkt_ptr->prop_time += edist / CLIGHT_PROP;
    move_pkt(pkt_ptr, edist, pkt_ptr->prop_time);
    edist = edist * 2.;

    // event occurs. Choose which event and call the appropriate subroutine.
    zrand = gsl_rng_uniform(rng);
    if (kap_compton > (zrand * kap_tot))
    {
      // Compton scattering.
      compton_scatter(pkt_ptr);
    }
    else if ((kap_compton + kap_photo_electric) > (zrand * kap_tot))
    {
      // Photo electric effect - makes it a k-packet for sure.
      pkt_ptr->type = TYPE_NTLEPTON;
      pkt_ptr->absorptiontype = -4;
      #ifndef FORCE_LTE
        //kgammadep[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif
      //pkt_ptr->type = TYPE_PRE_KPKT;
      //pkt_ptr->type = TYPE_GAMMA_KPKT;
      //if (tid == 0) nt_stat_from_gamma++;
      nt_stat_from_gamma++;
    }
    else if ((kap_compton + kap_photo_electric + kap_pair_prod) > (zrand * kap_tot))
    {
      // It's a pair production
      pair_prod(pkt_ptr);
    }
    else
    {
      printout("Failed to identify event. Gamma (1). kap_compton %g kap_photo_electric %g kap_tot %g zrand %g Abort.\n", kap_compton, kap_photo_electric, kap_tot, zrand);
      const int cellindex = pkt_ptr->where;
      printout(" /*cell[*/pkt_ptr->where].rho %g pkt_ptr->nu_cmf %g pkt_ptr->dir[0] %g pkt_ptr->dir[1] %g pkt_ptr->dir[2] %g pkt_ptr->pos[0] %g pkt_ptr->pos[1] %g pkt_ptr->pos[2] %g \n",get_rho(cell[cellindex].modelgridindex), pkt_ptr->nu_cmf,pkt_ptr->dir[0],pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2],pkt_ptr->pos[1],pkt_ptr->pos[2]);

      abort();
    }
  }
  else
  {
    printout("Failed to identify event. Gamma (2). edist %g, sdist %g, tdist %g Abort.\n", edist, sdist, tdist);
    abort();
  }
}


double get_gam_freq(const int n)
{
  if (n == RED_OF_LIST)
  {
    return 0.0;
  }

  // returns the frequency of line n
  enum radionuclides iso = gam_line_list.nuclidetype[n];
  const int lineid = gam_line_list.index[n];

  if (iso >= RADIONUCLIDE_COUNT || lineid >= gamma_spectra[iso].nlines)
  {
    printout("Unknown line. %d Abort.\n", n);
    printout("line_list->nuclidetype[n] %d line_list->index[n] %d\n", gam_line_list.nuclidetype[n], gam_line_list.index[n]);
    // printout(" %d %d \n", gam_line_list.nuclidetype[n], gam_line_list.index[n]);
    abort();
  }

  return gamma_spectra[iso].energy[lineid] / H;
}


int get_nul(double freq)
{
  const double freq_max = get_gam_freq(gam_line_list.total - 1);
  const double freq_min = get_gam_freq(0);

  if (freq > freq_max)
  {
    return (gam_line_list.total-1);
  }
  else if (freq < freq_min)
  {
    return RED_OF_LIST;
  }
  else
  {
    int too_high = gam_line_list.total - 1;
    int too_low = 0;

    while (too_high != too_low + 1)
  	{
  	  const int tryindex = (too_high + too_low) / 2;
  	  const double freq_try = get_gam_freq(tryindex);
  	  if (freq_try >= freq)
	    {
	      too_high = tryindex;
	    }
  	  else
	    {
	      too_low = tryindex;
	    }
  	}

    return too_low;
  }
}
