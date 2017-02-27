#include "sn3d.h"
#include "polarization.h"
#include "vectors.h"
#include "vpkt.h"


void escat_rpkt(PKT *pkt_ptr, double t_current)
{
  double dummy_dir[3], vel_vec[3], vel_rev[3];
  double old_dir_cmf[3],new_dir_cmf[3];
  double Qi, Ui, Inew, Unew, Qnew, Uold, Qold, I, Q, U;
  double mu,M,tsc,phisc;
  double i1,i2,cos2i1,sin2i1,cos2i2,sin2i2;
  double ref1[3],ref2[3];


  /// now make the packet a r-pkt and set further flags
  pkt_ptr->type = TYPE_RPKT;
  pkt_ptr->last_cross = NONE;  /// allow all further cell crossings


  get_velocity(pkt_ptr->pos, vel_vec, t_current);


  // Transform Stokes Parameters from the RF to the CMF

  Qi = pkt_ptr->stokes[1];
  Ui = pkt_ptr->stokes[2];

  frame_transform(pkt_ptr->dir,&Qi,&Ui,vel_vec,old_dir_cmf);


  // Outcoming direction. Compute the new cmf direction from the old direction and the scattering angles (see Kalos & Whitlock 2008)

  /* Assume dipole function (rejecton method, see Code & Whitney 1995) */
  #ifdef DIPOLE

    do {

        const double zrand = gsl_rng_uniform(rng);
        const double zrand2 = gsl_rng_uniform(rng);
        zrand3 = gsl_rng_uniform(rng);

        phisc = 2 * PI * zrand ;
        M = 2 * zrand2 - 1;
        mu = pow(M,2.) ;

        // NB: the rotational matrix R here is chosen in the clockwise direction ("+").
        // In Bulla+2015 equation (10) and (12) refer to the specific case shown in Fig.2 where the angle i2
        // is measured in the counter-clockwise direction. Therefore we use the clockwise rotation matrix but
        // with -i1. Here, instead, we calculate the angle in the clockwise direction from 0 to 2PI.
        // For instance, the i1 angle in Fig.2 of Bulla+2015 corresponds to 2PI-i1 here.
        // NB2: the i1 and i2 angles computed in the code (before and after scattering) are instead as in Bulla+2015
        p = (mu+1) + (mu-1) * ( cos(2*phisc) * Qi + sin(2*phisc) * Ui  );

        // generate a number between 0 and the maximum of the previous function (2)
        x = 2 * zrand3 ;
    }

    while (x>p);

  /* Assume isotropic scattering */
  #else

    const double zrand = gsl_rng_uniform(rng);
    const double zrand2 = gsl_rng_uniform(rng);

    M = 2. * zrand - 1 ;
    mu = pow(M,2.) ;
    phisc = 2 * PI * zrand2 ;

  #endif

  tsc = acos(M);

  if( fabs(old_dir_cmf[2]) < 0.99999 ) {

      new_dir_cmf[0] = sin(tsc)/sqrt(1.-pow(old_dir_cmf[2],2.)) * ( old_dir_cmf[1] * sin(phisc) - old_dir_cmf[0] * old_dir_cmf[2] * cos(phisc) ) + old_dir_cmf[0] * cos(tsc) ;
      new_dir_cmf[1] = sin(tsc)/sqrt(1-pow(old_dir_cmf[2],2.)) * ( - old_dir_cmf[0] * sin(phisc) - old_dir_cmf[1] * old_dir_cmf[2] * cos(phisc) ) + old_dir_cmf[1] * cos(tsc) ;
      new_dir_cmf[2] = sin(tsc) * cos(phisc) * sqrt(1-pow(old_dir_cmf[2],2.))  +  old_dir_cmf[2] * cos(tsc) ;

  }

  else {

      new_dir_cmf[0] = sin(tsc) * cos(phisc) ;
      new_dir_cmf[1] = sin(tsc) * sin(phisc) ;
      if(old_dir_cmf[2]>0 ) new_dir_cmf[2] = cos(tsc) ;
      else new_dir_cmf[2] = - cos(tsc) ;

  }


  // Need to rotate Stokes Parameters in the scattering plane

  meridian(old_dir_cmf,ref1,ref2);

  /* This is the i1 angle of Bulla+2015, obtained by computing the angle between the
     reference axes ref1 and ref2 in the meridian frame and the corresponding axes
     ref1_sc and ref2_sc in the scattering plane. It is the supplementary angle of the
     scatt angle phisc chosen in the rejection technique above (phisc+i1=180 or phisc+i1=540) */
  i1 = rot_angle(old_dir_cmf,new_dir_cmf,ref1,ref2);
  cos2i1 = cos(2 * i1) ;
  sin2i1 = sin(2 * i1) ;

  Qold = Qi * cos2i1 - Ui * sin2i1;
  Uold = Qi * sin2i1 + Ui * cos2i1;


  // Scattering

  mu = dot(old_dir_cmf,new_dir_cmf);

  Inew = 0.75 * ( (mu * mu + 1.0) + Qold * (mu * mu - 1.0) ) ;
  Qnew = 0.75 * ( (mu * mu - 1.0) + Qold * (mu * mu + 1.0) ) ;
  Unew = 1.5 * mu * Uold ;

  Qnew = Qnew / Inew;
  Unew = Unew / Inew;
  I = 1.0; // Inew / Inew


  // Need to rotate Stokes Parameters out of the scattering plane to the meridian frame (Clockwise rotation of PI-i2)

  meridian(new_dir_cmf,ref1,ref2);

  /* This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
     reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
     meridian frame. NB: we need to add PI to transform THETA to i2 */
  i2 = PI + rot_angle(new_dir_cmf,old_dir_cmf,ref1,ref2);
  cos2i2 = cos(2 * i2) ;
  sin2i2 = sin(2 * i2) ;

  Q = Qnew * cos2i2 + Unew * sin2i2;
  U = - Qnew * sin2i2 + Unew * cos2i2;


  // Transform Stokes Parameters from the CMF to the RF

  vel_rev[0] = - vel_vec[0] ;
  vel_rev[1] = - vel_vec[1] ;
  vel_rev[2] = - vel_vec[2] ;

  frame_transform(new_dir_cmf,&Q,&U,vel_rev,dummy_dir);


  pkt_ptr->stokes[0]=I;
  pkt_ptr->stokes[1]=Q;
  pkt_ptr->stokes[2]=U;



  // ---------------------- Update rest frame direction, frequency and energy --------------------

  pkt_ptr->dir[0] = dummy_dir[0];
  pkt_ptr->dir[1] = dummy_dir[1];
  pkt_ptr->dir[2] = dummy_dir[2];

  // Check unit vector.
  #ifdef DEBUG_ON
    if (fabs(vec_len(pkt_ptr->dir) - 1) > 1.e-6)
    {
      printout("[fatal] do_ma: Not a unit vector. Abort.\n");
      abort();
    }
  #endif

  // Finally we want to put in the rest frame energy and frequency.
  // And record that it's now a r-pkt.

  get_velocity(pkt_ptr->pos, vel_vec, t_current);
  pkt_ptr->nu_rf = pkt_ptr->nu_cmf / doppler(pkt_ptr->dir, vel_vec);

  #ifdef DEBUG_ON
    if (pkt_ptr->e_cmf >1e50)
    {
      printout("[fatal] emitt_rpkt: here %g\n", pkt_ptr->e_cmf);
      abort();
    }
  #endif

  pkt_ptr->e_rf = pkt_ptr->e_cmf * pkt_ptr->nu_rf /pkt_ptr->nu_cmf;

 }

