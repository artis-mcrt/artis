#include "sn3d.h"

/* Stuff for photo electric effect scattering. */

double sig_photo_electric(pkt_ptr,t_current)
     PKT *pkt_ptr;
     double t_current;
{
  double vel_vec[3];
  double doppler();
  int get_velocity();
  double sigma_cmf,sigma_rf;
  double sigma_cmf_si, sigma_cmf_cno, sigma_cmf_fe;
  /* Start by working out the x-section in the co-moving frame.*/

  int mgi = cell[pkt_ptr->where].modelgridindex;
  double rho = get_rho(mgi);
  double f_fe = get_ffe(mgi);

  if (gamma_grey < 0)
    {
      sigma_cmf_cno = 0.0448e-24 * pow(pkt_ptr->nu_cmf / 2.41326e19, -3.2);

      sigma_cmf_si = 1.16e-24 * pow(pkt_ptr->nu_cmf / 2.41326e19, -3.13);

      sigma_cmf_fe = 25.7e-24 * pow(pkt_ptr->nu_cmf / 2.41326e19, -3.0);

      /* 2.41326e19 = 100keV in frequency. */

      /* Now need to multiply by the particle number density. */

      sigma_cmf_cno *= rho * (1. - f_fe) / MH / 14;
       /* Assumes Z = 7. So mass = 14. */

      sigma_cmf_si *= rho / MH / 28;
      /* Assumes Z = 14. So mass = 28. */

      sigma_cmf_fe *= rho / MH / 56;
      /* Assumes Z = 28. So mass = 56. */

      sigma_cmf = (sigma_cmf_fe * f_fe) + (sigma_cmf_si * (1. - f_fe));

    }
  else
    {
      sigma_cmf = gamma_grey * rho;
    }

  /* Now need to convert between frames. */

  get_velocity(pkt_ptr->pos, vel_vec, t_current);
  sigma_rf = sigma_cmf * doppler(pkt_ptr->dir, vel_vec);

  return(sigma_rf);
}

/* Cross section for pair production. */

double
sig_pair_prod(pkt_ptr,t_current)
     PKT *pkt_ptr;
     double t_current;
{
  double vel_vec[3];
  double doppler();
  int get_velocity();
  double sigma_cmf,sigma_rf;
  double sigma_cmf_si, sigma_cmf_cno, sigma_cmf_fe;

  /* Start by working out the x-section in the co-moving frame.*/

  int mgi = cell[pkt_ptr->where].modelgridindex;
  double rho = get_rho(mgi);
  double f_fe = get_ffe(mgi);

  if (gamma_grey < 0)
    {
      /* 2.46636e+20 = 1022 keV in frequency */
      /* 3.61990e+20 = 1500 keV in frequency */

      if (pkt_ptr->nu_cmf > 2.46636e+20)
	{
	  if (pkt_ptr->nu_cmf > 3.61990e+20)
	    {

	      sigma_cmf_cno = (0.0481 + (0.301 * ((pkt_ptr->nu_cmf/2.41326e+20) - 1.5))) * 49.e-27;

	      sigma_cmf_si  = (0.0481 + (0.301 * ((pkt_ptr->nu_cmf/2.41326e+20) - 1.5))) * 196.e-27;

	      sigma_cmf_fe  = (0.0481 + (0.301 * ((pkt_ptr->nu_cmf/2.41326e+20) - 1.5))) * 784.e-27;
	    }
	  else
	    {
	      sigma_cmf_cno = 1.0063 * ((pkt_ptr->nu_cmf/2.41326e+20) - 1.022) * 49.e-27;

	      sigma_cmf_si  = 1.0063 * ((pkt_ptr->nu_cmf/2.41326e+20) - 1.022) * 196.e-27;

	      sigma_cmf_fe  = 1.0063 * ((pkt_ptr->nu_cmf/2.41326e+20) - 1.022) * 784.e-27;
	    }

	  /* Now need to multiply by the particle number density. */

	  sigma_cmf_cno *= rho * (1. - f_fe) / MH / 14;
	  /* Assumes Z = 7. So mass = 14. */

	  sigma_cmf_si *= rho / MH / 28;
	  /* Assumes Z = 14. So mass = 28. */

	  //sigma_cmf_fe *= cell[pkt_ptr->where].rho * cell[pkt_ptr->where].f_fe / MH / 56;
	  /* Assumes Z = 28. So mass = 56. */

	  sigma_cmf_fe *= rho / MH / 56;
	  sigma_cmf = (sigma_cmf_fe * f_fe) + (sigma_cmf_si * (1. - f_fe));

	}
      else
	{
	  sigma_cmf = 0.0;
	}
    }
  else
    {
      sigma_cmf = 0.0;
    }

  /* Now need to convert between frames. */

  get_velocity(pkt_ptr->pos, vel_vec, t_current);
  sigma_rf = sigma_cmf * doppler(pkt_ptr->dir, vel_vec);

  if (sigma_rf < 0)
    {
      printout("Negative pair production sigma. Setting to zero. Abort? %g\n", sigma_rf);
      sigma_rf = 0.0;
    }

  return(sigma_rf);
}

/* Routine to deal with pair production. */
int
pair_prod(pkt_ptr,t_current)
     PKT *pkt_ptr;
     double t_current;
{
  double zrand, zrand2, zrand3;
  double prob_gamma;
  double mu, phi, sintheta;
  int get_velocity();
  int angle_ab();
  double doppler();
  double dummy_dir[3], vel_vec[3];
  double vec_len();


  /* In pair prodiction, the original gamma makes an electron positron pair - kinetic energy equal to
     gamma ray energy - 1.022 MeV. We assume that the electron deposits any kinetic energy directly to
     the thermal pool. The positron annihilates with an electron locally making a pair of gamma rays
     at 0.511 MeV in the local cmf (isotropic). So all the thermal energy goes to the thermal pool
     immediately and the remainder goes into gamma-rays at 0.511 MeV. */

  prob_gamma = 1.022 * MEV / (H * pkt_ptr->nu_cmf);

  if (prob_gamma < 0)
    {
      printout("prob_gamma < 0. pair_prod. Abort. %g\n", prob_gamma);
      exit(0);
    }

  zrand = gsl_rng_uniform(rng);

  if (zrand > prob_gamma)
    {
      /* Convert it to an e-minus packet - actually it could be positron EK too, but this works
         for consistency with com_sca.*/
      pkt_ptr->type = TYPE_EMINUS;
      pkt_ptr->absorptiontype = -5;
    }
  else
    {
      /* The energy goes into emission at 511 keV. */
      pkt_ptr->nu_cmf = 0.511 * MEV / H;

      /*Now let's give the gamma ray a direction. */

      zrand2 = gsl_rng_uniform(rng);
      zrand3 = gsl_rng_uniform(rng);

      /*Assuming isotropic emission in cmf, use these two random numbers to set
	up a cmf direction in cos(theta) and phi. */

      mu = -1 + (2.*zrand2);
      phi = zrand3 * 2 * PI;
      sintheta = sqrt(1. - (mu * mu));

      pkt_ptr->dir[0] = sintheta * cos(phi);
      pkt_ptr->dir[1] = sintheta * sin(phi);
      pkt_ptr->dir[2] = mu;

      /* This direction is in the cmf - we want to convert it to the rest
	 frame - use aberation of angles. We want to convert from cmf to
	 rest so need -ve velocity. */

      get_velocity(pkt_ptr->pos, vel_vec, (-1*(t_current)));
      //negative time since we want the backwards transformation here

      angle_ab(pkt_ptr->dir, vel_vec, dummy_dir);

      pkt_ptr->dir[0] = dummy_dir[0];
      pkt_ptr->dir[1] = dummy_dir[1];
      pkt_ptr->dir[2] = dummy_dir[2];

      /*Check unit vector.*/
      if (fabs(vec_len(pkt_ptr->dir) - 1) > 1.e-8)
	{
	  printout("Not a unit vector. Abort.\n");
	  exit(0);
	}

      get_velocity(pkt_ptr->pos, vel_vec, pkt_ptr->tdecay);
      pkt_ptr->nu_rf = pkt_ptr->nu_cmf / doppler(pkt_ptr->dir, vel_vec);
      pkt_ptr->e_rf = pkt_ptr->e_cmf * pkt_ptr->nu_rf /pkt_ptr->nu_cmf;

      pkt_ptr->type = TYPE_GAMMA;
      pkt_ptr->last_cross = NONE;

    }

  return(0);
}
