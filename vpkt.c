#include "sn3d.h"
#include "exspec.h"
#include <string.h>
#include "atomic.h"
#include "boundary.h"
#include "ltepop.h"
#include "move.h"
#include "rpkt.h"
#include "vectors.h"
#include "vpkt.h"


/*******************************************************/
int
rlc_emiss_vpkt(PKT *pkt_ptr, double t_current, int bin, double *obs, int realtype)
{
  double vel_vec[3],old_dir_cmf[3],obs_cmf[3],vel_rev[3];
  double sdist;
  double kap_cont;
  int snext;
  bool end_packet;
  double t_future,t_arrive;
  double ldist;
  double nutrans;
  int element, ion, upper, lower;
  double A_ul,B_ul,B_lu;
  double n_u,n_l,tau,t_line;
  int mgi;
  double Qi,Ui,Qold,Uold,Inew,Qnew,Unew,I,Q,U,pn,prob;
  double mu,i1,i2,cos2i1,sin2i1,cos2i2,sin2i2;
  double ref1[3],ref2[3];

  int bin_range;

  PKT dummy;
  dummy = *pkt_ptr;
  PKT *dummy_ptr;
  dummy_ptr = &dummy;


  end_packet = false;
  tau = 0;
  sdist = 0;
  ldist = 0;
  nutrans = 0;
  t_future = t_current;

  dummy_ptr->dir[0] = obs[0] ;
  dummy_ptr->dir[1] = obs[1] ;
  dummy_ptr->dir[2] = obs[2] ;

  nvpkt++;          /* increment the number of virtual packet in the given timestep */

  get_velocity(pkt_ptr->pos, vel_vec, t_current);

  /* rf frequency and energy */
  dummy_ptr->nu_rf = dummy_ptr->nu_cmf / doppler(dummy_ptr->dir, vel_vec);
  dummy_ptr->e_rf = dummy_ptr->e_cmf * dummy_ptr->nu_rf /dummy_ptr->nu_cmf;


  Qi = dummy_ptr->stokes[1];
  Ui = dummy_ptr->stokes[2];


  // ------------ SCATTERING EVENT: dipole function --------------------

  if (realtype==1)
  {

      // Transform Stokes Parameters from the RF to the CMF

      frame_transform(pkt_ptr->dir,&Qi,&Ui,vel_vec,old_dir_cmf);

      // Need to rotate Stokes Parameters in the scattering plane

      angle_ab(dummy_ptr->dir, vel_vec, obs_cmf);

      meridian(old_dir_cmf,ref1,ref2);

      /* This is the i1 angle of Bulla+2015, obtained by computing the angle between the
         reference axes ref1 and ref2 in the meridian frame and the corresponding axes
         ref1_sc and ref2_sc in the scattering plane. */
      i1 = rot_angle(old_dir_cmf,obs_cmf,ref1,ref2);
      cos2i1 = cos(2 * i1) ;
      sin2i1 = sin(2 * i1) ;

      Qold = Qi * cos2i1 - Ui * sin2i1;
      Uold = Qi * sin2i1 + Ui * cos2i1;


      // Scattering

      mu = dot(old_dir_cmf,obs_cmf);

      pn=3./(16.*PI)*( 1+pow(mu,2.) + ( pow(mu,2.) - 1 ) * Qold );

      Inew = 0.75 * ( (mu * mu + 1.0) + Qold * (mu * mu - 1.0) ) ;
      Qnew = 0.75 * ( (mu * mu - 1.0) + Qold * (mu * mu + 1.0) ) ;
      Unew = 1.5 * mu * Uold ;

      Qnew = Qnew/Inew ;
      Unew = Unew/Inew ;
      I = Inew/Inew ;


      // Need to rotate Stokes Parameters out of the scattering plane to the meridian frame

      meridian(obs_cmf,ref1,ref2);

      /* This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
         reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
         meridian frame. NB: we need to add PI to transform THETA to i2 */
      i2 = PI + rot_angle(obs_cmf,old_dir_cmf,ref1,ref2);
      cos2i2 = cos(2 * i2) ;
      sin2i2 = sin(2 * i2) ;

      Q = Qnew * cos2i2 + Unew * sin2i2;
      U = - Qnew * sin2i2 + Unew * cos2i2;


      // Transform Stokes Parameters from the CMF to the RF

      vel_rev[0] = - vel_vec[0] ;
      vel_rev[1] = - vel_vec[1] ;
      vel_rev[2] = - vel_vec[2] ;

      frame_transform(obs_cmf,&Q,&U,vel_rev,obs);

  }


  // ------------ MACROATOM and KPKT: isotropic emission --------------------


  if  ( realtype == 2 || realtype == 3 ) {

      I=1; Q=0; U=0;
      pn = 1 / (4*PI) ;

  }


  // --------- compute the optical depth to boundary ----------------

  //printout("next_trans = %d \t nu_cmf = %g \n",dummy_ptr->next_trans,dummy_ptr->nu_cmf);

  mgi = cell[dummy_ptr->where].modelgridindex;

  while (!end_packet) {

      ldist = 0 ;

      sdist = boundary_cross(dummy_ptr, t_future, &snext); /* distance to the next cell */

      //printout("before \n");
      calculate_kappa_vpkt_cont(dummy_ptr, t_future);
      //printout("after \n");
      kap_cont = kappa_rpkt_cont[tid].total;

      tau += kap_cont * sdist * t_current * t_current * t_current / (t_future * t_future * t_future);

      //printout("--- New cell\ntau = %g ldist = %g \t sdist = %g \n",tau,ldist,sdist);

      /* kill vpkt with high optical depth */
      if (tau > tau_max_vpkt) return(0) ;

      while ( ldist < sdist ) {

          //printout("next_trans = %d \t nutrans = %g \t",dummy_ptr->next_trans,nutrans);

          nutrans = closest_transition(dummy_ptr);

          element = mastate[tid].element;
          ion = mastate[tid].ion;
          upper = mastate[tid].level;
          lower = mastate[tid].activatedfromlevel;
          A_ul = mastate[tid].einstein;

          if (nutrans > 0) {

              if (dummy_ptr->nu_cmf < nutrans) ldist = 0;
              else ldist = CLIGHT * t_current * (dummy_ptr->nu_cmf/nutrans - 1);

              if (ldist < 0.) printout("[warning] get_event: ldist < 0 %g\n",ldist);

              if (ldist > sdist) {      /* exit the while loop if you reach the boundary; go back to the previous transition to start next cell with the excluded line */

                  dummy_ptr->next_trans -= 1 ;
                  //printout("ldist > sdist : line in the next cell\n");
                  break ;
              }

              t_line = t_current + ldist / CLIGHT ;

              B_ul = CLIGHTSQUAREDOVERTWOH / pow(nutrans,3) * A_ul;
              B_lu = stat_weight(element,ion,upper)/stat_weight(element,ion,lower) * B_ul;

              n_u = calculate_exclevelpop(mgi,element,ion,upper);
              n_l = calculate_exclevelpop(mgi,element,ion,lower);

              tau += (B_lu*n_l - B_ul*n_u) * HCLIGHTOVERFOURPI * t_line;

              //printout("tau = %g ldist = %g \t sdist = %g \n",tau,ldist,sdist);


              /* kill vpkt with high optical depth */
              if (tau > tau_max_vpkt) return(0) ;


          }

      }


      // virtual packet is still at the starting position
      // move it to cell boundary and go to next cell
      //printf("I'm changing cell. I'm going from nu_cmf = %.e ",dummy_ptr->nu_cmf);

      t_future += (sdist / CLIGHT_PROP);
      move_pkt(dummy_ptr, sdist, t_future);

      //printout("About to change vpkt cell\n");
      change_cell_vpkt(dummy_ptr, snext, &end_packet, t_future);
      //printout("Completed change vpkt cell\n");

      //printout("dummy->nu_cmf = %g \n",dummy_ptr->nu_cmf);
      mgi = cell[dummy_ptr->where].modelgridindex;
      // break if you reach an empty cell
      if (mgi == MMODELGRID) break;

      /* kill vpkt with pass through a thick cell */
      if (modelgrid[mgi].thick == 1) return(0);


  }


  /* increment the number of escaped virtual packet in the given timestep */
  #ifdef _OPENMP
  #pragma omp critical
  #endif
  {
    if (realtype==1) nvpkt_esc1++;
    else if (realtype==2) nvpkt_esc2++;
    else if (realtype==3) nvpkt_esc3++;
  }

  // -------------- final stokes vector ---------------

  prob = pn * exp( - tau );

  I = I * prob ;
  Q = Q * prob ;
  U = U * prob ;

  dummy_ptr->stokes[0] = I;
  dummy_ptr->stokes[1] = Q;
  dummy_ptr->stokes[2] = U;

  if (I!=I || Q!=Q || U!=U) printout("Nan Number!! %g %g %g %g %g %g %g %g \n",I,Q,U,pn,tau,mu,i1,i2);

  /* bin on fly and produce file with spectrum */

  t_arrive = t_current - (dot(pkt_ptr->pos, dummy_ptr->dir)/CLIGHT_PROP);

  add_to_vspecpol(dummy_ptr,bin,t_arrive);


  // vpkt grid

  if (vgrid_flag==1) {

      for (bin_range=0;bin_range<Nrange_grid;bin_range++) {

          if ( dummy_ptr->nu_rf > nu_grid_min[bin_range] && dummy_ptr->nu_rf < nu_grid_max[bin_range] ) { // Frequency selection

              if ( t_arrive > tmin_grid && t_arrive < tmax_grid ) { // Time selection

                  add_to_vpkt_grid(dummy_ptr, vel_vec, bin_range, bin, obs);

              }
          }
      }
  }




  return(0);

}








///****************************************************************************
double rot_angle(double *n1, double *n2, double *ref1, double *ref2)
{
/* ------------- Rotation angle from the scattering plane --------------------------------------------- */
/* -------- We need to rotate Stokes Parameters to (or from) the scattering plane from (or to) -------- */
/* -------- the meridian frame such that Q=1 is in the scattering plane and along ref1 ---------------- */

    double ref1_sc[3],cos_stokes_rot_1, cos_stokes_rot_2, i;

    // ref1_sc is the ref1 axis in the scattering plane ref1 = n1 x ( n1 x n2 )
    ref1_sc[0] = n1[0] * dot(n1,n2) - n2[0];
    ref1_sc[1] = n1[1] * dot(n1,n2) - n2[1];
    ref1_sc[2] = n1[2] * dot(n1,n2) - n2[2];
    vec_norm(ref1_sc,ref1_sc);

    cos_stokes_rot_1 = dot(ref1_sc,ref1);
    cos_stokes_rot_2 = dot(ref1_sc,ref2);

    if (cos_stokes_rot_1<-1) cos_stokes_rot_1=-1;
    if (cos_stokes_rot_1>1) cos_stokes_rot_1=1;

    if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 > 0)) i = acos(cos_stokes_rot_1);
    if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 < 0)) i = 2 * acos(-1.) - acos(cos_stokes_rot_1);
    if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 < 0)) i = acos(-1.) + acos(fabs(cos_stokes_rot_1));
    if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 > 0)) i = acos(-1.) - acos(fabs(cos_stokes_rot_1));
    if (cos_stokes_rot_1 == 0) i = acos(-1.)/2.;
    if (cos_stokes_rot_2 == 0) i = 0.0 ;

    if (!isfinite(i)) printf("Warning NaN: %3.6f \t %3.6f \t %3.6f \n",cos_stokes_rot_1,cos_stokes_rot_2,acos(cos_stokes_rot_1));

    return i;
}





/* ----------------------- Routine to compute the meridian frame axes ref1 and ref2 ----------------------------------------*/

void meridian(double *n, double *ref1, double *ref2)
{
    // for ref_1 use (from triple product rule)

    ref1[0] = -1. * n[0] * n[2]/ sqrt(n[0]*n[0] + n[1]*n[1]);
    ref1[1] = -1. * n[1] * n[2]/ sqrt(n[0]*n[0] + n[1]*n[1]);
    ref1[2] = (1 - (n[2] * n[2]))/ sqrt(n[0]*n[0] + n[1]*n[1]);

    // for ref_2 use vector product of n_cmf with ref1

    ref2[0] = n[2] * ref1[1] - n[1] * ref1[2];
    ref2[1] = n[0] * ref1[2] - n[2] * ref1[0];
    ref2[2] = n[1] * ref1[0] - n[0] * ref1[1];


}


/* ----------------------- Routine to transform the Stokes Parameters from RF to CMF ----------------------------------------*/

void frame_transform(double *n_rf, double *Q, double *U, double *v, double *n_cmf)
{
    double rot_angle,cos2rot_angle,sin2rot_angle,p,Q0,U0;
    double ref1[3],ref2[3],e_rf[3],e_cmf[3];
    double e_cmf_ref1, e_cmf_ref2, theta_rot;

    // Meridian frame in the RF
    meridian(n_rf,ref1,ref2);

    Q0 = *Q;
    U0 = *U;

    // Compute polarisation (which is invariant)
    p = sqrt(Q0*Q0+U0*U0);

    // We want to compute the angle between ref1 and the electric field
    rot_angle=0;

    if (p > 0) {

        cos2rot_angle = Q0/p;
        sin2rot_angle = U0/p;

        if ((cos2rot_angle > 0) && (sin2rot_angle > 0)) rot_angle = acos(Q0 / p) / 2. ;
        if ((cos2rot_angle < 0) && (sin2rot_angle > 0)) rot_angle = (acos(-1.) - acos(fabs(Q0 / p))) / 2. ;
        if ((cos2rot_angle < 0) && (sin2rot_angle < 0)) rot_angle = (acos(-1.) + acos(fabs(Q0 / p))) / 2. ;
        if ((cos2rot_angle > 0) && (sin2rot_angle < 0)) rot_angle = (2. * acos(-1.) - acos(fabs(Q0 / p))) / 2. ;
        if (cos2rot_angle == 0) {
            rot_angle = 0.25 * acos(-1);
            if (U0 < 0) rot_angle = 0.75 * acos(-1);
        }
        if (sin2rot_angle == 0) {
            rot_angle = 0.0;
            if (Q0 < 0) rot_angle = 0.5 * acos(-1);

        }

    }

    // Define electric field by linear combination of ref1 and ref2 (using the angle just computed)
    e_rf[0] =  cos(rot_angle) * ref1[0] - sin(rot_angle) * ref2[0];
    e_rf[1] =  cos(rot_angle) * ref1[1] - sin(rot_angle) * ref2[1];
    e_rf[2] =  cos(rot_angle) * ref1[2] - sin(rot_angle) * ref2[2];

    // Aberration
    angle_ab(n_rf,v,n_cmf);

    // Lorentz transformation of E
    lorentz(e_rf,n_rf,v,e_cmf);

    // Meridian frame in the CMF
    meridian(n_cmf,ref1,ref2);

    // Projection of E onto ref1 and ref2
    e_cmf_ref1 = e_cmf[0] * ref1[0] + e_cmf[1] * ref1[1] + e_cmf[2] * ref1[2];
    e_cmf_ref2 = e_cmf[0] * ref2[0] + e_cmf[1] * ref2[1] + e_cmf[2] * ref2[2];

    // Compute the angle between ref1 and the electric field
    if ((e_cmf_ref1 > 0) && (e_cmf_ref2 < 0)) theta_rot = acos(e_cmf_ref1);
    if ((e_cmf_ref1 < 0) && (e_cmf_ref2 < 0)) theta_rot = acos(-1.) - acos(fabs(e_cmf_ref1));
    if ((e_cmf_ref1 < 0) && (e_cmf_ref2 > 0)) theta_rot = acos(-1.) + acos(fabs(e_cmf_ref1));
    if ((e_cmf_ref1 > 0) && (e_cmf_ref2 > 0)) theta_rot = 2*acos(-1.) - acos(e_cmf_ref1);
    if (e_cmf_ref1 == 0) theta_rot = acos(-1.)/2. ;
    if (e_cmf_ref2 == 0) theta_rot = 0.0 ;
    if (e_cmf_ref1 > 1) theta_rot = 0.0 ;
    if (e_cmf_ref1 < -1) theta_rot = acos(-1.) ;

    // Compute Stokes Parameters in the CMF
    *Q = cos(2 * theta_rot ) * p ;
    *U = sin(2 * theta_rot ) * p ;

}


/* ----------------------- Lorentz transformations from RF to CMF --------------------------------------------- */

void lorentz(double *e_rf, double *n_rf, double *v, double *e_cmf)
{
    double beta[3],e_par[3], e_perp[3], b_rf[3], b_par[3], b_perp[3], vsqr, gamma_rel, v_cr_b[3], v_cr_e[3], b_cmf[3];

    beta[0] = v[0] / CLIGHT ;
    beta[1] = v[1] / CLIGHT ;
    beta[2] = v[2] / CLIGHT ;
    vsqr = dot(beta,beta);

    gamma_rel = 1./(sqrt(1 - vsqr));


    e_par[0] = (e_rf[0]*beta[0] + e_rf[1]*beta[1] + e_rf[2]*beta[2]) * beta[0] / (vsqr);
    e_par[1] = (e_rf[0]*beta[0] + e_rf[1]*beta[1] + e_rf[2]*beta[2]) * beta[1] / (vsqr);
    e_par[2] = (e_rf[0]*beta[0] + e_rf[1]*beta[1] + e_rf[2]*beta[2]) * beta[2] / (vsqr);

    e_perp[0] = e_rf[0] - e_par[0];
    e_perp[1] = e_rf[1] - e_par[1];
    e_perp[2] = e_rf[2] - e_par[2];

    b_rf[0]=n_rf[1]*e_rf[2] - n_rf[2]*e_rf[1];
    b_rf[1]=n_rf[2]*e_rf[0] - n_rf[0]*e_rf[2];
    b_rf[2]=n_rf[0]*e_rf[1] - n_rf[1]*e_rf[0];

    b_par[0] = (b_rf[0]*beta[0] + b_rf[1]*beta[1] + b_rf[2]*beta[2]) * beta[0] / (vsqr);
    b_par[1] = (b_rf[0]*beta[0] + b_rf[1]*beta[1] + b_rf[2]*beta[2]) * beta[1] / (vsqr);
    b_par[2] = (b_rf[0]*beta[0] + b_rf[1]*beta[1] + b_rf[2]*beta[2]) * beta[2] / (vsqr);

    b_perp[0] = b_rf[0] - b_par[0];
    b_perp[1] = b_rf[1] - b_par[1];
    b_perp[2] = b_rf[2] - b_par[2];



    v_cr_b[0]=beta[1]*b_rf[2] - beta[2]*b_rf[1];
    v_cr_b[1]=beta[2]*b_rf[0] - beta[0]*b_rf[2];
    v_cr_b[2]=beta[0]*b_rf[1] - beta[1]*b_rf[0];

    v_cr_e[0]=beta[1]*e_rf[2] - beta[2]*e_rf[1];
    v_cr_e[1]=beta[2]*e_rf[0] - beta[0]*e_rf[2];
    v_cr_e[2]=beta[0]*e_rf[1] - beta[1]*e_rf[0];


    e_cmf[0] = e_par[0] + gamma_rel * (e_perp[0] + v_cr_b[0]);
    e_cmf[1] = e_par[1] + gamma_rel * (e_perp[1] + v_cr_b[1]);
    e_cmf[2] = e_par[2] + gamma_rel * (e_perp[2] + v_cr_b[2]);

    b_cmf[0] = b_par[0] + gamma_rel * (b_perp[0] - v_cr_e[0]);
    b_cmf[1] = b_par[1] + gamma_rel * (b_perp[1] - v_cr_e[1]);
    b_cmf[2] = b_par[2] + gamma_rel * (b_perp[2] - v_cr_e[2]);

    vec_norm(e_cmf,e_cmf);
    vec_norm(b_cmf,b_cmf);


}







/**********************************************************************/
/*Routine to add a packet to the outcoming spectrum.*/
int add_to_vspecpol(PKT *pkt_ptr, int bin, double t_arrive)
{
    /** Need to (1) decide which time bin to put it in and (2) which frequency bin. */

    double deltai,deltaq,deltau;
    int nt, nnu;

    /// Put this into the time grid.
    if (t_arrive > tmin && t_arrive < tmax)
    {

        nt = (log(t_arrive) - log(tmin)) / dlogt;
        if (pkt_ptr->nu_rf > nu_min_r && pkt_ptr->nu_rf < nu_max_r)
        {
            nnu = (log(pkt_ptr->nu_rf) - log(nu_min_r)) /  dlognu;

            deltai = pkt_ptr->stokes[0]*pkt_ptr->e_rf / vstokes_i[nt][bin].delta_t / delta_freq_vspec[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs * 4 * PI ;
            deltaq = pkt_ptr->stokes[1]*pkt_ptr->e_rf / vstokes_i[nt][bin].delta_t / delta_freq_vspec[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs * 4 * PI ;
            deltau = pkt_ptr->stokes[2]*pkt_ptr->e_rf / vstokes_i[nt][bin].delta_t / delta_freq_vspec[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs * 4 * PI ;

            vstokes_i[nt][bin].flux[nnu] += deltai;
            vstokes_q[nt][bin].flux[nnu] += deltaq;
            vstokes_u[nt][bin].flux[nnu] += deltau;

        }


    }

    return(0);

}



/**********************************************************************/
void init_vspecpol(void)
{
  for (int bin = 0;bin < MOBS; bin++)
  {
    /** start by setting up the time and frequency bins. */
    /** it is all done in terms of a logarithmic spacing in both t and nu - get the
     step sizes first. */

    dlogt = (log(tmax) - log(tmin))/VMTBINS;
    dlognu = (log(nu_max_r) - log(nu_min_r))/VMNUBINS;

    for (int n = 0; n < VMTBINS; n++)
    {
      vstokes_i[n][bin].lower_time = exp( log(tmin) + (n * (dlogt)));
      vstokes_i[n][bin].delta_t = exp( log(tmin) + ((n+1) * (dlogt))) - vstokes_i[n][bin].lower_time;

      for (int m = 0; m < VMNUBINS; m++)
      {
        lower_freq_vspec[m] = exp( log(nu_min_r) + (m * (dlognu)));
        delta_freq_vspec[m] = exp( log(nu_min_r) + ((m+1) * (dlognu))) - lower_freq_vspec[m];

        vstokes_i[n][bin].flux[m] = 0.0;
        vstokes_q[n][bin].flux[m] = 0.0;
        vstokes_u[n][bin].flux[m] = 0.0;
      }
    }
  }
}



/*******************************************************/
int write_vspecpol(FILE *specpol_file)
{
    for (int bin = 0; bin < Nobs; bin++)
    {

        fprintf(specpol_file, "%g ", 0.0);

        for (int l = 0; l < 3;l++) {
            for (int p = 0; p < VMTBINS; p++)
              fprintf(specpol_file, "%g ", (vstokes_i[p][bin].lower_time + (vstokes_i[p][bin].delta_t/2))/DAY);
        }

        fprintf(specpol_file, "\n");

        for (int m = 0; m < VMNUBINS; m++)
        {

            fprintf(specpol_file, "%g ", (lower_freq_vspec[m]+(delta_freq_vspec[m]/2)));

            // Stokes I
            for (int p = 0; p < VMTBINS; p++)
              fprintf(specpol_file, "%g ", vstokes_i[p][bin].flux[m]);

            // Stokes Q
            for (int p = 0; p < VMTBINS; p++)
              fprintf(specpol_file, "%g ", vstokes_q[p][bin].flux[m]);

            // Stokes U
            for (int p = 0; p < VMTBINS; p++)
              fprintf(specpol_file, "%g ", vstokes_u[p][bin].flux[m]);

            fprintf(specpol_file, "\n");

        }

        /*
         fclose(specpol_file);
         fclose(emissionpol_file);
         */


    }

    return 0;
}



/*******************************************************/
int read_vspecpol(FILE *specpol_file)
{
    int n,m,j,p,l,bin;
    float a,b,c;

    for (bin = 0; bin < Nobs; bin++) {

        // Initialise times and frequencies
        dlogt = (log(tmax) - log(tmin))/VMTBINS;
        dlognu = (log(nu_max_r) - log(nu_min_r))/VMNUBINS;

        for (n = 0; n < VMTBINS; n++)
        {
            vstokes_i[n][bin].lower_time = exp( log(tmin) + (n * (dlogt)));
            vstokes_i[n][bin].delta_t = exp( log(tmin) + ((n+1) * (dlogt))) - vstokes_i[n][bin].lower_time;

            for (m=0; m < VMNUBINS; m++)
            {
                lower_freq_vspec[m] = exp( log(nu_min_r) + (m * (dlognu)));
                delta_freq_vspec[m] = exp( log(nu_min_r) + ((m+1) * (dlognu))) - lower_freq_vspec[m];

            }
        }

        // Initialise I,Q,U fluxes (from temporary files)
        fscanf(specpol_file, "%g ", &a);

        for (l=0;l<3;l++) {
            for (p = 0; p < VMTBINS; p++) fscanf(specpol_file, "%g ", &b);
        }

        fscanf(specpol_file, "\n");

        for (j=0; j < VMNUBINS; j++)
        {

            fscanf(specpol_file, "%g ", &c);

            // Stokes I
            for (p = 0; p < VMTBINS; p++) fscanf(specpol_file, "%lg ", &vstokes_i[p][bin].flux[j]);

            // Stokes Q
            for (p = 0; p < VMTBINS; p++) fscanf(specpol_file, "%lg ", &vstokes_q[p][bin].flux[j]);

            // Stokes U
            for (p = 0; p < VMTBINS; p++) fscanf(specpol_file, "%lg ", &vstokes_u[p][bin].flux[j]);

            fscanf(specpol_file, "\n");

        }

        /*
         fclose(specpol_file);
         fclose(emissionpol_file);
         */


    }


    return(0);
}



//********************************** VPKT GRID *********************************************


/**********************************************************************/
void init_vpkt_grid(void)
{

    int n,m,bin_range,bin;

    double ybin, zbin;

    ybin = 2 * vmax / NY_VGRID ;
    zbin = 2 * vmax / NZ_VGRID ;


    for (bin=0;bin<MOBS;bin++) {

        for (bin_range=0;bin_range<MRANGE_GRID;bin_range++) {

            for (n = 0; n < NY_VGRID; n++) {

                for (m = 0; m < NZ_VGRID; m++) {

                    vgrid_i[n][m].flux[bin_range][bin] = 0.0;
                    vgrid_q[n][m].flux[bin_range][bin] = 0.0;
                    vgrid_u[n][m].flux[bin_range][bin] = 0.0;

                    vgrid_i[n][m].yvel[bin_range][bin] = vmax - (n+0.5) * ybin ;
                    vgrid_i[n][m].zvel[bin_range][bin] = vmax - (m+0.5) * zbin ;

                }
            }
        }
    }


}


/**********************************************************************/
/*Routine to add a packet to the outcoming spectrum.*/
int add_to_vpkt_grid(PKT *dummy_ptr, double *vel, int bin_range, int bin, double *obs)
{

    double vref1,vref2;
    double ybin, zbin;
    double nx,ny,nz;
    int nt,mt;

    // Observer orientation

    nx = obs[0] ;
    ny = obs[1] ;
    nz = obs[2] ;

    // Packet velocity

    /* if nobs = x , vref1 = vy and vref2 = vz */
    if (nx==1) {

        vref1 = vel[1] ;
        vref2 = vel[2] ;

    }

    /* if nobs = x , vref1 = vy and vref2 = vz */
    else if (nx==-1) {

        vref1 = - vel[1] ;
        vref2 = - vel[2] ;

    }

    /* Rotate velocity into projected area seen by the observer (see notes) */
    else {

        // Rotate velocity from (x,y,z) to (n_obs,ref1,ref2) so that x correspond to n_obs (see notes)
        vref1 = - ny * vel[0]  +  ( nx + nz * nz / (1 + nx) ) * vel[1]  -  ny * nz * (1 - nx) / sqrt(1 - nx * nx) * vel[2] ;
        vref2 = - nz * vel[0]  -  ny * nz * (1 - nx) / sqrt(1 - nx * nx) * vel[1]  +  ( nx + ny * ny / (1 + nx) ) * vel[2] ;

    }


    // Outside the grid
    if (fabs(vref1) >= vmax || fabs(vref2) >= vmax) return(0);


    // Bin size
    ybin = 2 * vmax / NY_VGRID ;
    zbin = 2 * vmax / NZ_VGRID ;

    // Grid cell
    nt = ( vmax - vref1 ) / ybin ;
    mt = ( vmax - vref2 ) / zbin ;

    // Add contribution
    if ( dummy_ptr->nu_rf > nu_grid_min[bin_range] && dummy_ptr->nu_rf < nu_grid_max[bin_range] ) {

        vgrid_i[nt][mt].flux[bin_range][bin] += dummy_ptr->stokes[0] * dummy_ptr->e_rf ;
        vgrid_q[nt][mt].flux[bin_range][bin] += dummy_ptr->stokes[1] * dummy_ptr->e_rf ;
        vgrid_u[nt][mt].flux[bin_range][bin] += dummy_ptr->stokes[2] * dummy_ptr->e_rf ;

    }


    return(0);

}


/*******************************************************/
int write_vpkt_grid(FILE *vpkt_grid_file)
{
    int n,m,bin_range,bin;


    for (bin=0;bin<Nobs;bin++) {

        for (bin_range=0;bin_range<Nrange_grid;bin_range++) {

             for (n = 0; n < NY_VGRID; n++) {

                 for (m = 0; m < NZ_VGRID; m++) {

                     fprintf(vpkt_grid_file, "%g ", vgrid_i[n][m].yvel[bin_range][bin]);
                     fprintf(vpkt_grid_file, "%g ", vgrid_i[n][m].zvel[bin_range][bin]);

                     fprintf(vpkt_grid_file, "%g ", vgrid_i[n][m].flux[bin_range][bin]);
                     fprintf(vpkt_grid_file, "%g ", vgrid_q[n][m].flux[bin_range][bin]);
                     fprintf(vpkt_grid_file, "%g ", vgrid_u[n][m].flux[bin_range][bin]);

                     fprintf(vpkt_grid_file, "\n");
                 }
             }
         }
     }


    return(0);
}


/*******************************************************/
int read_vpkt_grid(FILE *vpkt_grid_file)
{
    int n,m,bin_range,bin;

    for (bin=0;bin<Nobs;bin++) {

        for (bin_range=0;bin_range<Nrange_grid;bin_range++) {

            for (n = 0; n < NY_VGRID; n++) {

                for (m = 0; m < NZ_VGRID; m++) {

                  fscanf(vpkt_grid_file, "%lg ", &vgrid_i[n][m].yvel[bin_range][bin]);
                  fscanf(vpkt_grid_file, "%lg ", &vgrid_i[n][m].zvel[bin_range][bin]);

                  fscanf(vpkt_grid_file, "%lg ", &vgrid_i[n][m].flux[bin_range][bin]);
                  fscanf(vpkt_grid_file, "%lg ", &vgrid_q[n][m].flux[bin_range][bin]);
                  fscanf(vpkt_grid_file, "%lg ", &vgrid_u[n][m].flux[bin_range][bin]);

                  fscanf(vpkt_grid_file, "\n");
                }
            }
        }
    }



    return(0);
}
