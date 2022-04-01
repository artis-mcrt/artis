#include "sn3d.h"
#include "polarization.h"
#include "vpkt.h"


void escat_rpkt(PKT *pkt_ptr)
{
  /// now make the packet a r-pkt and set further flags
  pkt_ptr->type = TYPE_RPKT;
  pkt_ptr->last_cross = NONE;  /// allow all further cell crossings

  double vel_vec[3];
  get_velocity(pkt_ptr->pos, vel_vec, pkt_ptr->prop_time);

  // Transform Stokes Parameters from the RF to the CMF

  double Qi = pkt_ptr->stokes[1];
  double Ui = pkt_ptr->stokes[2];

  double old_dir_cmf[3];
  frame_transform(pkt_ptr->dir,&Qi,&Ui,vel_vec,old_dir_cmf);

  // Outcoming direction. Compute the new cmf direction from the old direction and the scattering angles (see Kalos & Whitlock 2008)
  double M = 0.;
  double mu = 0.;
  double phisc = 0.;

#ifdef DIPOLE
  // Assume dipole function (rejecton method, see Code & Whitney 1995)
  double p = 0.;
  double x = 0.;
  do
  {
    const double zrand = gsl_rng_uniform(rng);
    const double zrand2 = gsl_rng_uniform(rng);
    const double zrand3 = gsl_rng_uniform(rng);

    M = 2 * zrand - 1;
    mu = pow(M, 2.) ;
    phisc = 2 * PI * zrand2;

    // NB: the rotational matrix R here is chosen in the clockwise direction ("+").
    // In Bulla+2015 equation (10) and (12) refer to the specific case shown in Fig.2 where the angle i2
    // is measured in the counter-clockwise direction. Therefore we use the clockwise rotation matrix but
    // with -i1. Here, instead, we calculate the angle in the clockwise direction from 0 to 2PI.
    // For instance, the i1 angle in Fig.2 of Bulla+2015 corresponds to 2PI-i1 here.
    // NB2: the i1 and i2 angles computed in the code (before and after scattering) are instead as in Bulla+2015
    p = (mu + 1) + (mu - 1) * (cos(2 * phisc) * Qi + sin(2 * phisc) * Ui);

    // generate a number between 0 and the maximum of the previous function (2)
    x = 2 * zrand3;
  }
  while (x > p);

#else

  // Assume isotropic scattering
  const double zrand = gsl_rng_uniform(rng);
  const double zrand2 = gsl_rng_uniform(rng);

  M = 2. * zrand - 1;
  mu = pow(M, 2.);
  phisc = 2 * PI * zrand2;

#endif

  const double tsc = acos(M);
  double new_dir_cmf[3];

  if (fabs(old_dir_cmf[2]) < 0.99999)
  {
    new_dir_cmf[0] = sin(tsc) / sqrt(1. - pow(old_dir_cmf[2], 2.)) * (old_dir_cmf[1] * sin(phisc) - old_dir_cmf[0] * old_dir_cmf[2] * cos(phisc) ) + old_dir_cmf[0] * cos(tsc);
    new_dir_cmf[1] = sin(tsc) / sqrt(1 - pow(old_dir_cmf[2], 2.)) * (- old_dir_cmf[0] * sin(phisc) - old_dir_cmf[1] * old_dir_cmf[2] * cos(phisc) ) + old_dir_cmf[1] * cos(tsc);
    new_dir_cmf[2] = sin(tsc) * cos(phisc) * sqrt(1 - pow(old_dir_cmf[2],2.)) + old_dir_cmf[2] * cos(tsc);
  }
  else
  {
    new_dir_cmf[0] = sin(tsc) * cos(phisc);
    new_dir_cmf[1] = sin(tsc) * sin(phisc);
    if (old_dir_cmf[2] > 0)
    {
      new_dir_cmf[2] = cos(tsc);
    }
    else
    {
      new_dir_cmf[2] = - cos(tsc);
    }
  }

  // Need to rotate Stokes Parameters in the scattering plane

  double ref1[3];
  double ref2[3];
  meridian(old_dir_cmf,ref1,ref2);

  /* This is the i1 angle of Bulla+2015, obtained by computing the angle between the
     reference axes ref1 and ref2 in the meridian frame and the corresponding axes
     ref1_sc and ref2_sc in the scattering plane. It is the supplementary angle of the
     scatt angle phisc chosen in the rejection technique above (phisc+i1=180 or phisc+i1=540) */
  const double i1 = rot_angle(old_dir_cmf,new_dir_cmf,ref1,ref2);
  const double cos2i1 = cos(2 * i1);
  const double sin2i1 = sin(2 * i1);

  const double Qold = Qi * cos2i1 - Ui * sin2i1;
  const double Uold = Qi * sin2i1 + Ui * cos2i1;

  // Scattering

  mu = dot(old_dir_cmf,new_dir_cmf);

  const double Inew = 0.75 * ( (mu * mu + 1.0) + Qold * (mu * mu - 1.0) );
  double Qnew = 0.75 * ( (mu * mu - 1.0) + Qold * (mu * mu + 1.0) );
  double Unew = 1.5 * mu * Uold ;

  Qnew = Qnew / Inew;
  Unew = Unew / Inew;
  const double I = 1.0; // Inew / Inew

  // Need to rotate Stokes Parameters out of the scattering plane to the meridian frame (Clockwise rotation of PI-i2)

  meridian(new_dir_cmf, ref1, ref2);

  // This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
  //   reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
  //   meridian frame. NB: we need to add PI to transform THETA to i2
  const double i2 = PI + rot_angle(new_dir_cmf, old_dir_cmf, ref1, ref2);
  const double cos2i2 = cos(2 * i2);
  const double sin2i2 = sin(2 * i2);

  double Q = Qnew * cos2i2 + Unew * sin2i2;
  double U = - Qnew * sin2i2 + Unew * cos2i2;

  // Transform Stokes Parameters from the CMF to the RF
  double vel_rev[3];
  vel_rev[0] = - vel_vec[0];
  vel_rev[1] = - vel_vec[1];
  vel_rev[2] = - vel_vec[2];

  double dummy_dir[3];
  frame_transform(new_dir_cmf, &Q, &U, vel_rev, dummy_dir);

  pkt_ptr->stokes[0] = I;
  pkt_ptr->stokes[1] = Q;
  pkt_ptr->stokes[2] = U;


  // ---------------------- Update rest frame direction, frequency and energy --------------------

  pkt_ptr->dir[0] = dummy_dir[0];
  pkt_ptr->dir[1] = dummy_dir[1];
  pkt_ptr->dir[2] = dummy_dir[2];

  // Check unit vector.
  #ifdef DEBUG_ON
  if (fabs(vec_len(pkt_ptr->dir) - 1) > 1.e-6)
  {
    printout("WARNING: escat_rpkt: pkt_ptr->dir is not a unit vector. x %g y %g z %g length %.10f. Normalising to unit length...\n",
             pkt_ptr->dir[0], pkt_ptr->dir[1], pkt_ptr->dir[2], vec_len(pkt_ptr->dir));
    vec_norm(pkt_ptr->dir, pkt_ptr->dir);
    if (fabs(vec_len(pkt_ptr->dir) - 1) > 1.e-6)
    {
      printout("[fatal] escat_rpkt: After normalising: pkt_ptr->dir is still not a unit vector. x %g y %g z %g length %.10f\n",
               pkt_ptr->dir[0], pkt_ptr->dir[1], pkt_ptr->dir[2], vec_len(pkt_ptr->dir));
      abort();
    }
  }
  #endif

  // Finally we want to put in the rest frame energy and frequency.
  // And record that it's now a r-pkt.

  #ifdef DEBUG_ON
    if (pkt_ptr->e_cmf >1e52)
    {
      printout("[fatal] emitt_rpkt: here %g\n", pkt_ptr->e_cmf);
      abort();
    }
  #endif

  const double dopplerfactor = doppler_packetpos(pkt_ptr);
  pkt_ptr->nu_rf = pkt_ptr->nu_cmf / dopplerfactor;
  pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;
 }

