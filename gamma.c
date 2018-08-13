#include "sn3d.h"
#include "boundary.h"
#include "compton.h"
#include "emissivities.h"
#include "gamma.h"
#include "grey_emissivities.h"
#include "grid_init.h"
#include "photo_electric.h"
#include "move.h"
#include "vectors.h"

/* Material for handing gamma rays - creation and propagation. */

struct gamma_spec
{
  double *energy;
  double *probability;
  int nlines;
};

static struct gamma_spec gamma_spectra[RADIONUCLIDE_COUNT];

static LIST gam_line_list;

#define RED_OF_LIST -956  //must be negative


static double read_gamma_spectrum(enum radionuclides isotope, const char filename[50])
// reads in gamma_spectra and returns the average energy in gamma rays per nuclear decay
{
  assert(isotope < RADIONUCLIDE_COUNT);

  FILE *filein = fopen_required(filename, "r");
  int nlines = 0;
  fscanf(filein, "%d", &nlines);

  gamma_spectra[isotope].nlines = nlines;

  gamma_spectra[isotope].energy = (double *) malloc(nlines * sizeof(double));
  gamma_spectra[isotope].probability = (double *) malloc(nlines * sizeof(double));

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

  return E_gamma_avg;
}


static void read_decaydata(void)
{
  for (enum radionuclides iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
  {
    gamma_spectra[iso].nlines = 0;
    gamma_spectra[iso].energy = NULL;
    gamma_spectra[iso].probability = NULL;
  }

  E56NI = read_gamma_spectrum(NUCLIDE_NI56, "ni_lines.txt");

  E56CO_GAMMA = read_gamma_spectrum(NUCLIDE_CO56, "co_lines.txt");
  /// Average energy per gamma line of Co56 decay and positron annihilation
  /// For total deposited energy we need to add the kinetic energy per emitted positron
  E56CO = E56CO_GAMMA + 0.63 * MEV * 0.19;

  E48V = read_gamma_spectrum(NUCLIDE_V48, "v48_lines.txt");

  E48CR = read_gamma_spectrum(NUCLIDE_CR48, "cr48_lines.txt");

  E57NI_GAMMA = read_gamma_spectrum(NUCLIDE_NI57, "ni57_lines.txt");
  E57NI = E57NI_GAMMA + 0.354 * MEV * 0.436;

  E57CO = read_gamma_spectrum(NUCLIDE_CO57, "co57_lines.txt");
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
  for (enum radionuclides iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
  {
    total_lines += gamma_spectra[iso].nlines;
  }
  printout("total gamma-ray lines %d\n", total_lines);

  gam_line_list.total = total_lines;
  gam_line_list.nuclidetype = (enum radionuclides *) malloc(total_lines * sizeof(enum radionuclides));
  gam_line_list.index = (int *) malloc(total_lines * sizeof(int));

  double energy_last = 0.0;
  int next = -99;
  enum radionuclides next_type = -99;

  for (int i = 0; i < total_lines; i++)
  {
    double energy_try = 1.e50;

    for (enum radionuclides iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
    {
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
  /* Routine to choose which gamma ray line it'll be. */
  enum radionuclides iso;
  double E_gamma;  // Average energy per gamma line of a decay
  switch (pkt_ptr->type)
  {
    case TYPE_56NI_PELLET:
      iso = NUCLIDE_NI56;
      E_gamma = E56NI;
      break;

    case TYPE_56CO_PELLET:
      iso = NUCLIDE_CO56;
      E_gamma = E56CO_GAMMA;
      break;

    case TYPE_57NI_PELLET:
      iso = NUCLIDE_NI57;
      E_gamma = E57NI_GAMMA;
      break;

    case TYPE_57CO_PELLET:
      iso = NUCLIDE_CO57;
      E_gamma = E57CO;
      break;

    case TYPE_48CR_PELLET:
      iso = NUCLIDE_CR48;
      E_gamma = E48CR;
      break;

    case TYPE_48V_PELLET:
      iso = NUCLIDE_V48;
      E_gamma = E48V;
      break;

    default:
      printout("Unrecognised pellet. Abort.\n");
      abort();
  }

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
  #ifdef _OPENMP
    #pragma omp atomic
  #endif
  time_step[nts].pellet_decays++;

  // Start by getting the position of the pellet at the point of decay. Pellet
  // is moving with the matter.

  vec_scale(pkt_ptr->pos, pkt_ptr->tdecay / time_step[nts].start);

  // Now let's give the gamma ray a direction.

  const double zrand = gsl_rng_uniform(rng);
  const double zrand2 = gsl_rng_uniform(rng);

  // Assuming isotropic emission in cmf, use these two random numbers to set
  // up a cmf direction in cos(theta) and phi.

  const double mu = -1 + (2.*zrand);
  const double phi = zrand2 * 2 * PI;
  const double sintheta = sqrt(1. - (mu * mu));

  pkt_ptr->dir[0] = sintheta * cos(phi);
  pkt_ptr->dir[1] = sintheta * sin(phi);
  pkt_ptr->dir[2] = mu;

  /* This direction is in the cmf - we want to convert it to the rest
  frame - use aberation of angles. We want to convert from cmf to
  rest so need -ve velocity. */

  double vel_vec[3];
  get_velocity(pkt_ptr->pos, vel_vec, -1. * pkt_ptr->tdecay);
  //negative time since we want the backwards transformation here

  double dummy_dir[3];
  angle_ab(pkt_ptr->dir, vel_vec, dummy_dir);

  vec_copy(pkt_ptr->dir, dummy_dir);

  /*Check unit vector.*/
  if (fabs(vec_len(pkt_ptr->dir) - 1) > 1.e-8)
  {
    printout("Not a unit vector. Abort.\n");
    abort();
  }

  /* Now need to assign the frequency of the packet in the co-moving frame.*/

  choose_gamma_ray(pkt_ptr);

  /* Finally we want to put in the rest frame energy and frequency. And record
  that it's now a gamma ray.*/

  const double dopplerfactor = doppler_packetpos(pkt_ptr, pkt_ptr->tdecay);
  pkt_ptr->nu_rf = pkt_ptr->nu_cmf / dopplerfactor;
  pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;

  pkt_ptr->type = TYPE_GAMMA;
  pkt_ptr->last_cross = NONE;

  /* initialise polarisation information */
  pkt_ptr->stokes[0] = 1.0;
  pkt_ptr->stokes[1] = pkt_ptr->stokes[2] = 0.0;
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


double do_gamma(PKT *restrict pkt_ptr, double t1, double t2)
// Now routine for moving a gamma packet. Idea is that we have as input
// a gamma packet with known properties at time t1 and we want to follow it
// until time t2.
{
  double t_current = t1; //this will keep track of time in the calculation

  bool end_packet = false; //tells us when to stop working on this packet
  while (!end_packet)
  {
    // Assign optical depth to next physical event. And start counter of
    // optical depth for this path.
    double zrand = gsl_rng_uniform(rng);
    const double tau_next = -1. * log(zrand);
    const double tau_current = 0.0;

    // Start by finding the distance to the crossing of the grid cell
    // boundaries. sdist is the boundary distance and snext is the
    // grid cell into which we pass.

    int snext;
    double sdist = boundary_cross(pkt_ptr, t_current, &snext);

    const double maxsdist = (grid_type == GRID_SPHERICAL1D) ? 2 * rmax * t_current / tmin : rmax * t_current / tmin;
    if (sdist > maxsdist)
    {
      printout("Unreasonably large sdist (gamma). Abort. %g %g %g\n", rmax, t_current/tmin, sdist);
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
      kap_compton = sig_comp(pkt_ptr, t_current);
    }

    const double kap_photo_electric = sig_photo_electric(pkt_ptr, t_current);
    const double kap_pair_prod = sig_pair_prod(pkt_ptr, t_current);
    const double kap_tot = kap_compton + kap_photo_electric + kap_pair_prod;

    // So distance before physical event is...

    double edist = (tau_next - tau_current) / kap_tot;

    if (edist < 0)
    {
      printout("Negative distance (edist). Abort. \n");
      abort();
    }

    // Find how far it can travel during the time inverval.

    double tdist = (t2 - t_current) * CLIGHT_PROP;

    if (tdist < 0)
    {
      printout("Negative distance (tdist). Abort. \n");
      abort();
    }

    //printout("sdist, tdist, edist %g %g %g\n",sdist, tdist, edist);

    if ((sdist < tdist) && (sdist < edist))
    {
      sdist = sdist / 2.;
      t_current += sdist / CLIGHT_PROP;
      move_pkt(pkt_ptr, sdist, t_current);

      // Move it into the new cell.
      if (kap_tot > 0)
      {
        if (do_comp_est)
        {
          sdist = sdist * 2.;
          compton_emiss_cont(pkt_ptr, sdist, t_current);
          pp_emiss_cont(pkt_ptr, sdist, t_current);
          sdist = sdist / 2.;
        }
        if (do_rlc_est != 0)
        {
          sdist = sdist * 2.;
          rlc_emiss_gamma(pkt_ptr, sdist, t_current);
          sdist = sdist / 2.;
        }
      }

      t_current += sdist / CLIGHT_PROP;
      move_pkt(pkt_ptr, sdist, t_current);
      sdist = sdist * 2.;
      if (snext != pkt_ptr->where)
      {
        change_cell(pkt_ptr, snext, &end_packet, t_current);
      }
    }
    else if ((tdist < sdist) && (tdist < edist))
    {
      // Doesn't reach boundary.
      tdist = tdist / 2.;
      t_current += tdist  / CLIGHT_PROP;
      move_pkt(pkt_ptr, tdist, t_current);

      if (kap_tot > 0)
      {
        if (do_comp_est)
        {
          tdist = tdist * 2.;
          compton_emiss_cont(pkt_ptr, tdist, t_current);
          pp_emiss_cont(pkt_ptr, tdist, t_current);
          tdist = tdist / 2.;
        }
        if (do_rlc_est != 0)
        {
          tdist = tdist * 2.;
          rlc_emiss_gamma(pkt_ptr, tdist, t_current);
          tdist = tdist / 2.;
        }
      }
      t_current = t2;
      move_pkt(pkt_ptr, tdist, t_current);
      tdist = tdist * 2.;
      end_packet = true;
    }
    else if ((edist < sdist) && (edist < tdist))
    {
      edist = edist / 2.;
      t_current += edist / CLIGHT_PROP;
      move_pkt(pkt_ptr, edist, t_current);
      if (kap_tot > 0)
      {
        if (do_comp_est)
        {
          edist = edist * 2.;
          compton_emiss_cont(pkt_ptr, edist, t_current);
          pp_emiss_cont(pkt_ptr, edist, t_current);
          edist = edist / 2.;
        }
        if (do_rlc_est != 0)
        {
          edist = edist * 2.;
          rlc_emiss_gamma(pkt_ptr, edist, t_current);
          edist = edist / 2.;
        }
      }
      t_current += edist / CLIGHT_PROP;
      move_pkt(pkt_ptr, edist, t_current);
      edist = edist * 2.;

      // event occurs. Choose which event and call the appropriate subroutine.
      zrand = gsl_rng_uniform(rng);
      if (kap_compton > (zrand * kap_tot))
      {
        // Compton scattering.
        compton_scatter(pkt_ptr, t_current);
        if (pkt_ptr->type != TYPE_GAMMA)
        {
          // It's not a gamma ray any more - return.
          return t_current;
        }
      }
      else if ((kap_compton + kap_photo_electric) > (zrand*kap_tot))
      {
        // Photo electric effect - makes it a k-packet for sure.
        pkt_ptr->type = TYPE_KPKT;
        pkt_ptr->absorptiontype = -4;
        #ifndef FORCE_LTE
          //kgammadep[pkt_ptr->where] += pkt_ptr->e_cmf;
        #endif
        //pkt_ptr->type = TYPE_PRE_KPKT;
        //pkt_ptr->type = TYPE_GAMMA_KPKT;
        //if (tid == 0) k_stat_from_gamma++;
        k_stat_from_gamma++;
        return t_current;
      }
      else if ((kap_compton + kap_photo_electric + kap_pair_prod) > (zrand*kap_tot))
      {
        // It's a pair production
        pair_prod(pkt_ptr, t_current);
        if (pkt_ptr->type != TYPE_GAMMA)
        {
          // It's not a gamma ray any more - return.
          return t_current;
        }
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
  return PACKET_SAME;
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
  	  const int try = (too_high + too_low)/2;
  	  const double freq_try = get_gam_freq(try);
  	  if (freq_try >= freq)
	    {
	      too_high = try;
	    }
  	  else
	    {
	      too_low = try;
	    }
  	}

    return too_low;
  }
}
