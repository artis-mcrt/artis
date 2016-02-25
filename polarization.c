#include "sn3d.h"


///****************************************************************************
void escat_rpkt(PKT *pkt_ptr, double t_current)
{
  double dummy_dir[3], vel_vec[3];
  double old_dir_cmf[3],new_dir_cmf[3],dir_perp[3];
  double old_dir_rf[3];

  /// now make the packet a r-pkt and set further flags
  pkt_ptr->type = TYPE_RPKT;
  pkt_ptr->last_cross = NONE;  /// allow all further cell crossings

  /// First get the incoming direction in the cmf and store it

  /// We have incoming dir in rf - we want to convert it to the cmf
  /// - use aberation of angles.


  ///get_velocity(pkt_ptr->pos, vel_vec, t_current);
  ///angle_ab(pkt_ptr->dir, vel_vec, old_dir_cmf);

  ///trying for now to work in rf only
  old_dir_rf[0] = pkt_ptr->dir[0];
  old_dir_rf[1] = pkt_ptr->dir[1];
  old_dir_rf[2] = pkt_ptr->dir[2];

  /// Need to assign a new direction. Assume isotropic emission in the cmf
  double zrand = gsl_rng_uniform(rng);
  double zrand2 = gsl_rng_uniform(rng);

  double mu = -1 + (2.*zrand);
  double phi = zrand2 * 2 * PI;
  double sintheta = sqrt(1. - (mu * mu));

  new_dir_cmf[0] = sintheta * cos(phi);
  new_dir_cmf[1] = sintheta * sin(phi);
  new_dir_cmf[2] = mu;

  //printout("[debug] pkt_ptr->dir in CMF: %g %g %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);

  /// This direction is in the cmf - we want to convert it to the rest
  /// frame - use aberation of angles. We want to convert from cmf to
  /// rest so need -ve velocity.
  get_velocity(pkt_ptr->pos, vel_vec, (-1*(t_current)));
  ///negative time since we want the backwards transformation here

  angle_ab(new_dir_cmf, vel_vec, dummy_dir);
  pkt_ptr->dir[0] = dummy_dir[0];
  pkt_ptr->dir[1] = dummy_dir[1];
  pkt_ptr->dir[2] = dummy_dir[2];
  //printout("[debug] pkt_ptr->dir in RF: %g %g %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);

  /// Check unit vector.
  #ifdef DEBUG_ON
    if (fabs(vec_len(pkt_ptr->dir) - 1) > 1.e-8)
    {
      printout("[fatal] do_ma: Not a unit vector. Abort.\n");
      abort();
    }
  #endif

  /// Finally we want to put in the rest frame energy and frequency. And record
  /// that it's now a r-pkt.
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


  ///for polarization want to get a direction which is perpendiculat to scattering plane
  //printout("rf pkt direction before scattering %g, %g, %g\n",old_dir_rf[0],old_dir_rf[1],old_dir_rf[2]);
  //printout("rf pkt direction after scattering %g, %g, %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);

  cross_prod(old_dir_rf,pkt_ptr->dir,dir_perp);
  //printout("old pol reference dir %g, %g, %g\n",pkt_ptr->pol_dir[0],pkt_ptr->pol_dir[1],pkt_ptr->pol_dir[2]);
  //printout("new pol reference dir %g, %g, %g\n",dir_perp[0],dir_perp[1],dir_perp[2]);
  vec_norm(dir_perp,dir_perp);
  //printout("new pol reference dir after normalisation %g, %g, %g\n",dir_perp[0],dir_perp[1],dir_perp[2]);

  ///now we want to know the angle that rotated the coordinate system for polarization from the stored reference direction pol_dir to the new one defined by dir_perp. This comes from their dor and cross products. The angle is "gamma", following Hillier 91
  double cos_gamma = dot(pkt_ptr->pol_dir,dir_perp);
  cross_prod(dir_perp,pkt_ptr->pol_dir,dummy_dir);
  double sin_gamma = dot(dummy_dir,old_dir_rf);
  //printout("cos_gamma %g, sin_gamma %g\n",cos_gamma,sin_gamma);

  if ((fabs(((cos_gamma*cos_gamma) + (sin_gamma*sin_gamma))-1.0)) > 1.e-6)
  {
    printout("Polarization rotation angle misbehaving. %g %g %g\n", cos_gamma, sin_gamma, ((cos_gamma*cos_gamma) + (sin_gamma*sin_gamma)));
  }

  ///transform Q and U by rotation
  double sin_2gamma = 2* sin_gamma*cos_gamma;
  double cos_2gamma = 2*(cos_gamma*cos_gamma) - 1.0;

  //printout("stokes qu old, %g %g\n",pkt_ptr->stokes_qu[0],pkt_ptr->stokes_qu[1]);
  double Qold = (pkt_ptr->stokes_qu[0] * cos_2gamma) + (pkt_ptr->stokes_qu[1] * sin_2gamma);
  double Uold = (pkt_ptr->stokes_qu[1] * cos_2gamma) - (pkt_ptr->stokes_qu[0] * sin_2gamma);
  //printout("stokes qu old rotated, %g %g\n",pkt_ptr->stokes_qu[0],pkt_ptr->stokes_qu[1]);

  ///now apply the scattering matrix
  ///conceptually, I is 1.0; we want to get Inew, Qnew and Unew then divide Qnew and Unew by Inew

  mu = dot(old_dir_rf,pkt_ptr->dir);
  double mu2 = mu * mu;
  double Inew = 0.75 * ((1. + mu2) + ((mu2 - 1.0)*Qold));
  double Qnew = 0.75 * ((mu2 - 1.0) + ((mu2 + 1.0)*Qold));
  double Unew = 1.5 * mu * Uold;

  pkt_ptr->stokes_qu[0] = Qnew / Inew;
  pkt_ptr->stokes_qu[1] = Unew / Inew;
  //printout("mu, mu2 %g %g\n",mu,mu2);
  //printout("stokes qu new, %g %g\n",pkt_ptr->stokes_qu[0],pkt_ptr->stokes_qu[1]);

  pkt_ptr->pol_dir[0] = dir_perp[0];
  pkt_ptr->pol_dir[1] = dir_perp[1];
  pkt_ptr->pol_dir[2] = dir_perp[2];
 }
