#include "sn3d.h"
#include "grid_init.h"
#include "compton.h"
#include "vectors.h"


// Stuff for compton scattering.

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


double sig_comp(const PKT *pkt_ptr, double t_current)
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

  double vel_vec[3];
  get_velocity(pkt_ptr->pos, vel_vec, t_current);
  double sigma_rf = sigma_cmf * doppler(pkt_ptr->dir, vel_vec);

  return sigma_rf;
}


static double choose_f(double xx, double zrand)
// To choose the value of f to integrate to - idea is we want
//   sigma_compton_partial(xx,f) = zrand. */
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
    const double try = sigma_compton_partial(xx, ftry);
    //printout("ftry %g %g %g %g %g\n",ftry, f_min, f_max, try, norm);
    if (try > norm)
    {
      f_max = ftry;
      err = (try - norm) / norm;
    }
    else
    {
      f_min = ftry;
      err = (norm - try) / norm;
    }
    //      printout("error %g\n",err);
    count++;
    if (count == 1000)
    {
      printout("Compton hit 1000 tries. %g %g %g %g %g\n", f_max, f_min, ftry, try, norm);
    }
  }

  return ftry;
}


static double thomson_angle(void)
{
  // For Thomson scattering we can get the new angle from a random number very easily.

  const double zrand = gsl_rng_uniform(rng);

  const double B_coeff = (8. * zrand) - 4.;

  double t_coeff = sqrt( (B_coeff * B_coeff) + 4);
  t_coeff = t_coeff - B_coeff;
  t_coeff = t_coeff / 2;
  t_coeff = pow(t_coeff, (1 / 3));

  const double mu = (1 / t_coeff) - t_coeff;

  if (fabs(mu) > 1)
  {
    printout("Error in Thomson. Abort.\n");
    abort();
  }

  return mu;
}


void com_sca(PKT *pkt_ptr, double t_current)
// Routine to deal with physical Compton scattering event.
{
  double f;
  double final_dir[3];
  double prob_gamma;

  //  printout("Compton scattering.\n");

  const double xx = H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT;

  /* It is known that a Compton scattering event is going to take place.
     We need to do two things - (1) decide whether to convert energy
     to electron or leave as gamma (2) decide properties of new packet.*/

  /* The probability of giving energy to electron is related to the
     energy change of the gamma ray. This is equivalent to the choice of
     scattering angle. Probability of scattering into particular angle
     (i.e. final energy) is related to the partial cross-section.*/

  /* Choose a random number to get the energy. Want to find the
   factor by which the energy changes "f" such that
   sigma_partial/sigma_tot = zrand */

  double zrand = gsl_rng_uniform(rng);

  if (xx <  THOMSON_LIMIT)
  {
    f = 1.0; //no energy loss
    prob_gamma = 1.0;
  }
  else
  {
    f = choose_f(xx, zrand);

    /* Check that f lies between 1.0 and (2xx  + 1) */

    if ((f < 1) || (f > (2 * xx + 1)))
    {
      printout("Compton f out of bounds. Abort.\n");
      abort();
    }

    /* Prob of keeping gamma ray is...*/

    prob_gamma = 1. / f;
  }

  zrand = gsl_rng_uniform(rng);
  if (zrand < prob_gamma)
  {
    // It stays as a gamma ray. Change frequency and direction in
    // co-moving frame then transfer back to rest frame.

    pkt_ptr->nu_cmf = pkt_ptr->nu_cmf / f; // reduce frequency

    // The packet has stored the direction in the rest frame.
    // Use aberation of angles to get this into the co-moving frame.

    double vel_vec[3];
    get_velocity(pkt_ptr->pos, vel_vec, t_current);

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

    double test = dot(new_dir, new_dir);
    if (fabs(1. - test) > 1.e-8)
    {
      printout("Not a unit vector - Compton. Abort. %g %g %g\n", f, xx, test);
      printout("new_dir %g %g %g\n", new_dir[0], new_dir[1], new_dir[2]);
      printout("cmf_dir %g %g %g\n", cmf_dir[0], cmf_dir[1], cmf_dir[2]);
      printout("cos_theta %g", cos_theta);
      abort();
    }

    test = dot(new_dir, cmf_dir);
    if (fabs(test - cos_theta) > 1.e-8)
    {
      printout("Problem with angle - Compton. Abort.\n");
      abort();
    }

    // Now convert back again.

    get_velocity(pkt_ptr->pos, vel_vec, (-1 * t_current));
    angle_ab(new_dir, vel_vec, final_dir);

    pkt_ptr->dir[0] = final_dir[0];
    pkt_ptr->dir[1] = final_dir[1];
    pkt_ptr->dir[2] = final_dir[2];

    // It now has a rest frame direction and a co-moving frequency.
    //  Just need to set the rest frame energy.

    get_velocity(pkt_ptr->pos, vel_vec, t_current);

    const double dopplerfactor = 1 / doppler(pkt_ptr->dir, vel_vec);

    pkt_ptr->nu_rf = pkt_ptr->nu_cmf * dopplerfactor;
    pkt_ptr->e_rf = pkt_ptr->e_cmf * dopplerfactor;

    pkt_ptr->last_cross = NONE; // allow it to re-cross a boundary
  }
  else
  {
    // It's converted to an e-minus packet.
    pkt_ptr->type = TYPE_NTLEPTON;
    pkt_ptr->absorptiontype = -3;
  }
}
