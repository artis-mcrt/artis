#include "sn3d.h"
#include "vpkt.h"
#include "rpkt.h"
#include "boundary.h"
#include "ltepop.h"
#include "vectors.h"
#include "update_grid.h"
#include "atomic.h"
#include <cstring>


struct vspecpol
{
  double flux[VMNUBINS];
  float lower_time;
  float delta_t;
};


struct vspecpol **vstokes_i;
struct vspecpol **vstokes_q;
struct vspecpol **vstokes_u;

float lower_freq_vspec[VMNUBINS];
float delta_freq_vspec[VMNUBINS];

// --------- INPUT PARAMETERS -----------

int Nobs;  // Number of observer directions
int Nspectra; // Number of virtual packet spectra per observer direction (total + elements switched off)
double *nz_obs_vpkt;
double *phiobs;
double tmin_vspec_input;
double tmax_vspec_input;
int Nrange;

double numin_vspec_input[MRANGE];
double numax_vspec_input[MRANGE];
double cell_is_optically_thick_vpkt;
double tau_max_vpkt;
double *exclude;
double *tau_vpkt;

// --------- VPKT GRID -----------

struct vgrid
{
  double *flux[MRANGE_GRID];
  double *yvel[MRANGE_GRID];
  double *zvel[MRANGE_GRID];
};

struct vgrid vgrid_i[NY_VGRID][NZ_VGRID];
struct vgrid vgrid_q[NY_VGRID][NZ_VGRID];
struct vgrid vgrid_u[NY_VGRID][NZ_VGRID];

int Nrange_grid;
double tmin_grid;
double tmax_grid;
double nu_grid_min[MRANGE_GRID];
double nu_grid_max[MRANGE_GRID];
int vgrid_flag;
double dlogt_vspec;
double dlognu_vspec;

int realtype;

// number of virtual packets in a given timestep
int nvpkt;

// number of escaped virtual packet in a given timestep (with tau < tau_max)
int nvpkt_esc1; // electron scattering event
int nvpkt_esc2; // kpkt deactivation
int nvpkt_esc3; // macroatom deactivation


void rlc_emiss_vpkt(PKT *pkt_ptr, double t_current, int bin, double *obs, int realtype)
{
  double vel_vec[3],old_dir_cmf[3],obs_cmf[3],vel_rev[3];
  double s_cont;
  double kap_cont,kap_cont_nobf,kap_cont_noff,kap_cont_noes;
  int snext;
  double t_arrive;
  int element, ion, upper, lower;
  double A_ul,B_ul,B_lu;
  double n_u,n_l,t_line;
  int mgi;
  double Qold,Uold,Inew,Qnew,Unew,Itmp,Qtmp,Utmp,I,Q,U,pn,prob;
  double mu,i1,i2,cos2i1,sin2i1,cos2i2,sin2i2;
  double ref1[3],ref2[3];
  int anumber,tau_flag=0;

  int bin_range;

  PKT dummy;
  dummy = *pkt_ptr;
  PKT *dummy_ptr;
  dummy_ptr = &dummy;

  bool end_packet = false;
  double sdist = 0;
  double ldist = 0;
  double t_future = t_current;

  for (int ind = 0; ind < Nspectra; ind++)
  {
    tau_vpkt[ind] = 0;
  }

  dummy_ptr->dir[0] = obs[0];
  dummy_ptr->dir[1] = obs[1];
  dummy_ptr->dir[2] = obs[2];

  nvpkt++;          // increment the number of virtual packet in the given timestep

  get_velocity(pkt_ptr->pos, vel_vec, t_current);

  // rf frequency and energy
  dummy_ptr->nu_rf = dummy_ptr->nu_cmf / doppler(dummy_ptr->dir, vel_vec);
  dummy_ptr->e_rf = dummy_ptr->e_cmf * dummy_ptr->nu_rf /dummy_ptr->nu_cmf;

  double Qi = dummy_ptr->stokes[1];
  double Ui = dummy_ptr->stokes[2];

  // ------------ SCATTERING EVENT: dipole function --------------------

  if (realtype == 1)
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

    pn = 3./(16.*PI)*( 1+pow(mu,2.) + ( pow(mu,2.) - 1 ) * Qold );

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

  if  (realtype == 2 || realtype == 3)
  {
    I = 1;
    Q = 0;
    U = 0;
    pn = 1 / (4*PI);
  }

  // --------- compute the optical depth to boundary ----------------

  mgi = grid::get_cell_modelgridindex(dummy_ptr->where);

  while (end_packet == false)
  {
    ldist = 0 ;

    /* distance to the next cell */
    sdist = boundary_cross(dummy_ptr, t_future, &snext);
    s_cont = sdist * t_current * t_current * t_current / (t_future * t_future * t_future);

    calculate_kappa_rpkt_cont(dummy_ptr, &globals::kappa_rpkt_cont[tid]);

    kap_cont = globals::kappa_rpkt_cont[tid].total;
    kap_cont_nobf = kap_cont - globals::kappa_rpkt_cont[tid].bf;
    kap_cont_noff = kap_cont - globals::kappa_rpkt_cont[tid].ff;
    kap_cont_noes = kap_cont - globals::kappa_rpkt_cont[tid].es;

    for (int ind = 0; ind < Nspectra; ind++)
    {
      if (exclude[ind] == -2)
        tau_vpkt[ind] += kap_cont_nobf * s_cont;
      else if (exclude[ind] == -3)
        tau_vpkt[ind] += kap_cont_noff * s_cont;
      else if (exclude[ind] == -4)
        tau_vpkt[ind] += kap_cont_noes * s_cont;
      else
        tau_vpkt[ind] += kap_cont * s_cont;
    }

    // kill vpkt with high optical depth
    tau_flag = check_tau(tau_vpkt, &tau_max_vpkt);
    if (tau_flag == 0)
    {
      return;
    }

    while (ldist < sdist)
    {
      //printout("next_trans = %d \t nutrans = %g \t",dummy_ptr->next_trans,nutrans);

      const int lineindex = closest_transition(dummy_ptr->nu_cmf, dummy_ptr->next_trans);

      const double nutrans = globals::linelist[lineindex].nu;

      if (lineindex >= 0)
      {
        element = globals::linelist[lineindex].elementindex;
        ion = globals::linelist[lineindex].ionindex;
        upper = globals::linelist[lineindex].upperlevelindex;
        lower = globals::linelist[lineindex].lowerlevelindex;
        A_ul = globals::linelist[lineindex].einstein_A;

        anumber = get_element(element);

        dummy_ptr->next_trans = lineindex + 1;

        if (dummy_ptr->nu_cmf < nutrans)
        {
          ldist = 0;
        }
        else
        {
          ldist = CLIGHT * t_current * (dummy_ptr->nu_cmf / nutrans - 1);
        }

        if (ldist < 0.) printout("[warning] get_event: ldist < 0 %g\n",ldist);

        if (ldist > sdist)
        {      /* exit the while loop if you reach the boundary; go back to the previous transition to start next cell with the excluded line */

          dummy_ptr->next_trans -= 1 ;
          //printout("ldist > sdist : line in the next cell\n");
          break ;
        }

        t_line = t_current + ldist / CLIGHT ;

        B_ul = CLIGHTSQUAREDOVERTWOH / pow(nutrans, 3) * A_ul;
        B_lu = stat_weight(element,ion,upper) / stat_weight(element,ion,lower) * B_ul;

        n_u = calculate_exclevelpop(mgi,element,ion,upper);
        n_l = calculate_exclevelpop(mgi,element,ion,lower);

        // Check on the element to exclude
        // NB: ldist before need to be computed anyway (I want to move the packets to the
        // line interaction point even if I don't interact)

        for (int ind = 0; ind < Nspectra; ind++)
        {
          // If exclude[ind]==-1, I do not include line opacity
          if (exclude[ind] != -1 && (anumber != exclude[ind]))
          {
            tau_vpkt[ind] += (B_lu*n_l - B_ul*n_u) * HCLIGHTOVERFOURPI * t_line;
          }
        }

        /* kill vpkt with high optical depth */
        tau_flag = check_tau(tau_vpkt, &tau_max_vpkt) ;
        if (tau_flag == 0)
        {
          return;
        }
      }
      else
      {
        dummy_ptr->next_trans = globals::nlines + 1;  ///helper variable to overcome numerical problems after line scattering
      }

    }

    // virtual packet is still at the starting position
    // move it to cell boundary and go to next cell
    //printf("I'm changing cell. I'm going from nu_cmf = %.e ",dummy_ptr->nu_cmf);

    t_future += (sdist / globals::CLIGHT_PROP);
    dummy_ptr->prop_time = t_future;
    move_pkt(dummy_ptr, sdist, t_future);

    //printout("About to change vpkt cell\n");
    change_cell(dummy_ptr, snext, t_future);
    end_packet = (dummy_ptr->type == TYPE_ESCAPE);
    //printout("Completed change vpkt cell\n");

    //printout("dummy->nu_cmf = %g \n",dummy_ptr->nu_cmf);
    mgi = grid::get_cell_modelgridindex(dummy_ptr->where);
    // break if you reach an empty cell
    if (mgi == grid::get_npts_model()) break;

    /* kill vpkt with pass through a thick cell */
    if (grid::modelgrid[mgi].thick == 1)
    {
      return;
    }
  }


  /* increment the number of escaped virtual packet in the given timestep */
  if (realtype==1) safeincrement(nvpkt_esc1);
  else if (realtype==2) safeincrement(nvpkt_esc2);
  else if (realtype==3) safeincrement(nvpkt_esc3);

  // -------------- final stokes vector ---------------

  for (int ind = 0; ind < Nspectra; ind++)
  {
    // printout("bin %d spectrum %d tau_vpkt %g\n", bin, ind, tau_vpkt[ind]);
    prob = pn * exp( - tau_vpkt[ind] ) ;

    Itmp = I * prob;
    Qtmp = Q * prob;
    Utmp = U * prob;

    dummy_ptr->stokes[0] = Itmp;
    dummy_ptr->stokes[1] = Qtmp;
    dummy_ptr->stokes[2] = Utmp;

    if (Itmp!=Itmp || Qtmp!=Qtmp || Utmp!=Utmp) printout("Nan Number!! %g %g %g %g %g %g %g %g \n",Itmp,Qtmp,Utmp,pn,tau_vpkt[ind],mu,i1,i2);

    /* bin on fly and produce file with spectrum */

    t_arrive = t_current - (dot(pkt_ptr->pos, dummy_ptr->dir) / globals::CLIGHT_PROP) ;

    add_to_vspecpol(dummy_ptr,bin,ind,t_arrive);
  }

  // vpkt grid

  if (vgrid_flag == 1)
  {
    prob = pn * exp( - tau_vpkt[0] ) ;

    Itmp = I * prob;
    Qtmp = Q * prob;
    Utmp = U * prob;

    dummy_ptr->stokes[0] = Itmp;
    dummy_ptr->stokes[1] = Qtmp;
    dummy_ptr->stokes[2] = Utmp;

    for (bin_range = 0; bin_range < Nrange_grid; bin_range++)
    {
      if (dummy_ptr->nu_rf > nu_grid_min[bin_range] && dummy_ptr->nu_rf < nu_grid_max[bin_range] )
      { // Frequency selection
        if (t_arrive > tmin_grid && t_arrive < tmax_grid)
        { // Time selection
          add_to_vpkt_grid(dummy_ptr, vel_vec, bin_range, bin, obs);
        }
      }
    }
  }
}


/* Virtual packet is killed when tau reaches tau_max_vpkt for ALL the different setups
E.g. imagine that a packet in the first setup (all elements included) reaches tau = tau_max_vpkt
because of the element Zi. If we remove Zi, tau now could be lower than tau_max_vpkt and could
thus contribute to the spectrum. */
int check_tau(double *tau, double *tau_max)
{
  int count = 0;

  for (int i = 0; i < Nspectra; i++)
  {
    if (tau[i] > *tau_max)
      count += 1;
  }

  if (count == Nspectra)
    return 0;
  else
    return 1;
}


// Routine to add a packet to the outcoming spectrum.
void add_to_vspecpol(PKT *pkt_ptr, int bin, int ind, double t_arrive)
{
  // Need to decide in which (1) time and (2) frequency bin the vpkt is escaping

  const int ind_comb = Nspectra * bin + ind;

  /// Put this into the time grid.
  if (t_arrive > tmin_vspec && t_arrive < tmax_vspec)
  {
    int nt = (log(t_arrive) - log(tmin_vspec)) / dlogt_vspec;
    if (pkt_ptr->nu_rf > numin_vspec && pkt_ptr->nu_rf < numax_vspec)
    {
      const int nnu = (log(pkt_ptr->nu_rf) - log(numin_vspec)) / dlognu_vspec;
      const double pktcontrib = pkt_ptr->e_rf / vstokes_i[nt][ind_comb].delta_t / delta_freq_vspec[nnu] / 4.e12 / PI / PARSEC / PARSEC / globals::nprocs * 4 * PI;

      safeadd(vstokes_i[nt][ind_comb].flux[nnu], pkt_ptr->stokes[0] * pktcontrib);
      safeadd(vstokes_q[nt][ind_comb].flux[nnu], pkt_ptr->stokes[1] * pktcontrib);
      safeadd(vstokes_u[nt][ind_comb].flux[nnu], pkt_ptr->stokes[2] * pktcontrib);
    }
  }
}


void init_vspecpol(void)
{
  vstokes_i = (struct vspecpol **) malloc(VMTBINS * sizeof(struct vspecpol *));
  vstokes_q = (struct vspecpol **) malloc(VMTBINS * sizeof(struct vspecpol *));
  vstokes_u = (struct vspecpol **) malloc(VMTBINS * sizeof(struct vspecpol *));

  const int indexmax = Nspectra * Nobs;
  for (int p = 0; p < VMTBINS; p++)
  {
    vstokes_i[p] = (struct vspecpol *) malloc(indexmax * sizeof(struct vspecpol));
    vstokes_q[p] = (struct vspecpol *) malloc(indexmax * sizeof(struct vspecpol));
    vstokes_u[p] = (struct vspecpol *) malloc(indexmax * sizeof(struct vspecpol));
  }

  for (int ind_comb = 0; ind_comb < indexmax; ind_comb++)
  {
    // start by setting up the time and frequency bins.
    // it is all done interms of a logarithmic spacing in both t and nu - get the
    // step sizes first.

    dlogt_vspec = (log(tmax_vspec) - log(tmin_vspec)) / VMTBINS;
    dlognu_vspec = (log(numax_vspec) - log(numin_vspec)) / VMNUBINS;

    for (int n = 0; n < VMTBINS; n++)
    {
      vstokes_i[n][ind_comb].lower_time = exp(log(tmin_vspec) + (n * (dlogt_vspec)));
      vstokes_i[n][ind_comb].delta_t = exp(log(tmin_vspec) + ((n + 1) * (dlogt_vspec))) - vstokes_i[n][ind_comb].lower_time;

      for (int m = 0; m < VMNUBINS; m++)
      {
        lower_freq_vspec[m] = exp(log(numin_vspec) + (m * (dlognu_vspec)));
        delta_freq_vspec[m] = exp(log(numin_vspec) + ((m + 1) * (dlognu_vspec))) - lower_freq_vspec[m];

        vstokes_i[n][ind_comb].flux[m] = 0.0;
        vstokes_q[n][ind_comb].flux[m] = 0.0;
        vstokes_u[n][ind_comb].flux[m] = 0.0;
      }
    }
  }
}


void write_vspecpol(FILE *specpol_file)
{
  for (int ind_comb = 0; ind_comb < (Nobs * Nspectra); ind_comb++)
  {
    fprintf(specpol_file, "%g ", 0.);

    for (int l = 0; l < 3; l++)
    {
      for (int p = 0; p < VMTBINS; p++)
      {
        fprintf(specpol_file, "%g ", (vstokes_i[p][ind_comb].lower_time + (vstokes_i[p][ind_comb].delta_t / 2.)) / DAY);
      }
    }

    fprintf(specpol_file, "\n");

    for (int m = 0; m < VMNUBINS; m++)
    {
      fprintf(specpol_file, "%g ", (lower_freq_vspec[m] + (delta_freq_vspec[m] / 2.)));

      // Stokes I
      for (int p = 0; p < VMTBINS; p++)
      {
        fprintf(specpol_file, "%g ", vstokes_i[p][ind_comb].flux[m]);
      }

      // Stokes Q
      for (int p = 0; p < VMTBINS; p++)
      {
        fprintf(specpol_file, "%g ", vstokes_q[p][ind_comb].flux[m]);
      }

      // Stokes U
      for (int p = 0; p < VMTBINS; p++)
      {
        fprintf(specpol_file, "%g ", vstokes_u[p][ind_comb].flux[m]);
      }

      fprintf(specpol_file, "\n");
    }

    /*
     fclose(specpol_file);
     fclose(emissionpol_file);
     */

  }
}



void read_vspecpol(int my_rank, int nts)
{
  char filename[1024];

  if (nts % 2 == 0)
    sprintf(filename, "vspecpol_%d_%d_odd.tmp", 0, my_rank);
  else
    sprintf(filename, "vspecpol_%d_%d_even.tmp", 0, my_rank);

  FILE *vspecpol_file = fopen_required(filename, "rb");

  float a,b,c;

  for (int ind_comb = 0; ind_comb < (Nobs * Nspectra); ind_comb++)
  {
    // Initialise times and frequencies
    dlogt_vspec = (log(tmax_vspec) - log(tmin_vspec))/VMTBINS;
    dlognu_vspec = (log(numax_vspec) - log(numin_vspec))/VMNUBINS;

    for (int n = 0; n < VMTBINS; n++)
    {
      vstokes_i[n][ind_comb].lower_time = exp(log(tmin_vspec) + (n * (dlogt_vspec)));
      vstokes_i[n][ind_comb].delta_t = exp(log(tmin_vspec) + ((n+1) * (dlogt_vspec))) - vstokes_i[n][ind_comb].lower_time;

      for (int m = 0; m < VMNUBINS; m++)
      {
        lower_freq_vspec[m] = exp(log(numin_vspec) + (m * (dlognu_vspec)));
        delta_freq_vspec[m] = exp(log(numin_vspec) + ((m + 1) * (dlognu_vspec))) - lower_freq_vspec[m];
      }
    }

    // Initialise I,Q,U fluxes (from temporary files)
    fscanf(vspecpol_file, "%g ", &a);

    for (int l = 0; l < 3; l++) {
      for (int p = 0; p < VMTBINS; p++)
      {
        fscanf(vspecpol_file, "%g ", &b);
      }
    }

    fscanf(vspecpol_file, "\n");

    for (int j = 0; j < VMNUBINS; j++)
    {
      fscanf(vspecpol_file, "%g ", &c);

      // Stokes I
      for (int p = 0; p < VMTBINS; p++) fscanf(vspecpol_file, "%lg ", &vstokes_i[p][ind_comb].flux[j]);

      // Stokes Q
      for (int p = 0; p < VMTBINS; p++) fscanf(vspecpol_file, "%lg ", &vstokes_q[p][ind_comb].flux[j]);

      // Stokes U
      for (int p = 0; p < VMTBINS; p++) fscanf(vspecpol_file, "%lg ", &vstokes_u[p][ind_comb].flux[j]);

      fscanf(vspecpol_file, "\n");
    }
  }

  fclose(vspecpol_file);
}


void init_vpkt_grid(void)
{
  const double ybin = 2 * globals::vmax / NY_VGRID;
  const double zbin = 2 * globals::vmax / NZ_VGRID;

  for (int n = 0; n < NY_VGRID; n++)
  {
    for (int m = 0; m < NZ_VGRID; m++)
    {
      for (int bin_range = 0; bin_range < MRANGE_GRID; bin_range++)
      {
        vgrid_i[n][m].flux[bin_range] = (double *) malloc(Nobs * sizeof(double));
        vgrid_i[n][m].yvel[bin_range] = (double *) malloc(Nobs * sizeof(double));
        vgrid_i[n][m].zvel[bin_range] = (double *) malloc(Nobs * sizeof(double));

        vgrid_q[n][m].flux[bin_range] = (double *) malloc(Nobs * sizeof(double));
        vgrid_q[n][m].yvel[bin_range] = (double *) malloc(Nobs * sizeof(double));
        vgrid_q[n][m].zvel[bin_range] = (double *) malloc(Nobs * sizeof(double));

        vgrid_u[n][m].flux[bin_range] = (double *) malloc(Nobs * sizeof(double));
        vgrid_u[n][m].yvel[bin_range] = (double *) malloc(Nobs * sizeof(double));
        vgrid_u[n][m].zvel[bin_range] = (double *) malloc(Nobs * sizeof(double));

        vgrid_i[n][m].flux[bin_range] = (double *) malloc(Nobs * sizeof(double));
        for (int bin = 0; bin < Nobs; bin++)
        {
          vgrid_i[n][m].flux[bin_range][bin] = 0.0;
          vgrid_q[n][m].flux[bin_range][bin] = 0.0;
          vgrid_u[n][m].flux[bin_range][bin] = 0.0;

          vgrid_i[n][m].yvel[bin_range][bin] = globals::vmax - (n+0.5) * ybin;
          vgrid_i[n][m].zvel[bin_range][bin] = globals::vmax - (m+0.5) * zbin;
        }
      }
    }
  }
}


// Routine to add a packet to the outcoming spectrum.
void add_to_vpkt_grid(PKT *dummy_ptr, double *vel, int bin_range, int bin, double *obs)
{
  double vref1,vref2;
  double ybin, zbin;
  double nx,ny,nz;
  int nt,mt;

  // Observer orientation

  nx = obs[0];
  ny = obs[1];
  nz = obs[2];

  // Packet velocity

  /* if nobs = x , vref1 = vy and vref2 = vz */
  if (nx == 1)
  {
    vref1 = vel[1];
    vref2 = vel[2];
  }
  /* if nobs = x , vref1 = vy and vref2 = vz */
  else if (nx == -1)
  {
    vref1 = - vel[1];
    vref2 = - vel[2];
  }

  // Rotate velocity into projected area seen by the observer (see notes)
  else
  {
    // Rotate velocity from (x,y,z) to (n_obs,ref1,ref2) so that x correspond to n_obs (see notes)
    vref1 = - ny * vel[0]  +  ( nx + nz * nz / (1 + nx) ) * vel[1]  -  ny * nz * (1 - nx) / sqrt(1 - nx * nx) * vel[2];
    vref2 = - nz * vel[0]  -  ny * nz * (1 - nx) / sqrt(1 - nx * nx) * vel[1]  +  ( nx + ny * ny / (1 + nx) ) * vel[2];
  }

  // Outside the grid
  if (fabs(vref1) >= globals::vmax || fabs(vref2) >= globals::vmax) return;

  // Bin size
  ybin = 2 * globals::vmax / NY_VGRID;
  zbin = 2 * globals::vmax / NZ_VGRID;

  // Grid cell
  nt = ( globals::vmax - vref1 ) / ybin;
  mt = ( globals::vmax - vref2 ) / zbin;

  // Add contribution
  if (dummy_ptr->nu_rf > nu_grid_min[bin_range] && dummy_ptr->nu_rf < nu_grid_max[bin_range])
  {
    safeadd(vgrid_i[nt][mt].flux[bin_range][bin], dummy_ptr->stokes[0] * dummy_ptr->e_rf);
    safeadd(vgrid_q[nt][mt].flux[bin_range][bin], dummy_ptr->stokes[1] * dummy_ptr->e_rf);
    safeadd(vgrid_u[nt][mt].flux[bin_range][bin], dummy_ptr->stokes[2] * dummy_ptr->e_rf);
  }
}


void write_vpkt_grid(FILE *vpkt_grid_file)
{
  for (int bin = 0; bin < Nobs; bin++)
  {
    for (int bin_range = 0; bin_range < Nrange_grid; bin_range++)
    {
      for (int n = 0; n < NY_VGRID; n++)
      {
        for (int m = 0; m < NZ_VGRID; m++)
        {
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
}


void read_vpkt_grid(FILE *vpkt_grid_file)
{
  for (int bin = 0; bin < Nobs; bin++)
  {
    for (int bin_range = 0; bin_range < Nrange_grid; bin_range++)
    {
      for (int n = 0; n < NY_VGRID; n++)
      {
        for (int m = 0; m < NZ_VGRID; m++)
        {
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
}


void read_parameterfile_vpkt(void)
{
  FILE *input_file = fopen_required("vpkt.txt", "r");

  // Nobs
  fscanf(input_file, "%d", &Nobs);

  printout("vpkt.txt: Nobs %d directions\n", Nobs);

  // nz_obs_vpkt. Cos(theta) to the observer. A list in the case of many observers
  nz_obs_vpkt = (double *) malloc(Nobs * sizeof(double));
  for (int i = 0; i < Nobs; i++)
  {
    fscanf(input_file, "%lg", &nz_obs_vpkt[i]);

    if (fabs(nz_obs_vpkt[i]) > 1)
    {
      printout("Wrong observer direction \n");
      exit(0);
    }
    else if (nz_obs_vpkt[i] == 1)
    {
      nz_obs_vpkt[i] = 0.9999;
    }
    else if (nz_obs_vpkt[i] == -1)
    {
      nz_obs_vpkt[i] = -0.9999;
    }
  }

  // phi to the observer (degrees). A list in the case of many observers
  phiobs = (double *) malloc(Nobs * sizeof(double));
  for (int i = 0; i < Nobs; i++)
  {
    double phi_degrees = 0.;
    fscanf(input_file, "%lg \n", &phi_degrees);
    phiobs[i] = phi_degrees * PI / 180.;

    printout("vpkt.txt:   direction %d costheta %g phi %g (%g degrees)\n", i, nz_obs_vpkt[i], phiobs[i], phi_degrees);
  }

  // Nspectra opacity choices (i.e. Nspectra spectra for each observer)
  int nspectra_customlist_flag;
  fscanf(input_file, "%d ", &nspectra_customlist_flag);

  if (nspectra_customlist_flag != 1)
  {
    Nspectra = 1;
    exclude = (double *) malloc(Nspectra * sizeof(double));

    exclude[0] = 0;
  }
  else
  {
    fscanf(input_file, "%d ", &Nspectra);
    exclude = (double *) malloc(Nspectra * sizeof(double));

    for (int i = 0; i < Nspectra; i++)
    {
      fscanf(input_file, "%lg ", &exclude[i]);

      // The first number should be equal to zero!
      assert_always(exclude[0] == 0); // The first spectrum should allow for all opacities (exclude[i]=0)
    }
  }

  printout("vpkt.txt: Nspectra %d per observer\n", Nspectra);
  tau_vpkt = (double *) malloc(Nspectra * sizeof(double));

  // time window. If dum4=1 it restrict vpkt to time windown (dum5,dum6)
  int override_tminmax = 0;
  double vspec_tmin_in_days = 0.;
  double vspec_tmax_in_days = 0.;
  fscanf(input_file, "%d %lg %lg \n", &override_tminmax, &vspec_tmin_in_days, &vspec_tmax_in_days);

  printout("vpkt: compiled with tmin_vspec %.1fd tmax_vspec %1.fd VMTBINS %d\n", tmin_vspec / DAY, tmax_vspec / DAY, VMTBINS);
  if (override_tminmax == 1)
  {
    tmin_vspec_input = vspec_tmin_in_days * DAY;
    tmax_vspec_input = vspec_tmax_in_days * DAY;
    printout("vpkt.txt: tmin_vspec_input %.1fd, tmax_vspec_input %.1fd\n", tmin_vspec_input / DAY, tmax_vspec_input / DAY);
  }
  else
  {
    tmin_vspec_input = tmin_vspec;
    tmax_vspec_input = tmax_vspec;
    printout("vpkt.txt: tmin_vspec_input %.1fd, tmax_vspec_input %.1fd (inherited from tmin_vspec and tmax_vspec)\n", tmin_vspec_input / DAY, tmax_vspec_input / DAY);
  }

  assert_always(tmax_vspec_input >= tmin_vspec);
  assert_always(tmax_vspec_input <= tmax_vspec);

  // frequency window. dum4 restrict vpkt to a frequency range, dum5 indicates the number of ranges,
  // followed by a list of ranges (dum6,dum7)
  int flag_custom_freq_ranges = 0;
  fscanf(input_file, "%d ", &flag_custom_freq_ranges);

  printout("vpkt: compiled with VMNUBINS %d\n", VMNUBINS);
  assert_always(numax_vspec > numin_vspec);
  printout("vpkt: compiled with numax_vspec %g lambda_min %g Å\n", numax_vspec, 1e8 * CLIGHT / numax_vspec);
  printout("vpkt: compiled with numin_vspec %g lambda_max %g Å\n", numin_vspec, 1e8 * CLIGHT / numin_vspec);
  if (flag_custom_freq_ranges == 1)
  {
    fscanf(input_file, "%d ", &Nrange);
    assert_always(Nrange <= MRANGE);

    printout("vpkt.txt: Nrange %d frequency intervals per spectrum per observer\n", Nrange);

    for (int i = 0; i < Nrange; i++)
    {
      double lmin_vspec_input = 0.;
      double lmax_vspec_input = 0.;
      fscanf(input_file, "%lg %lg", &lmin_vspec_input, &lmax_vspec_input);

      numin_vspec_input[i] = CLIGHT / (lmax_vspec_input * 1e-8);
      numax_vspec_input[i] = CLIGHT / (lmin_vspec_input * 1e-8);
      printout("vpkt.txt:   range %d lambda [%g, %g] Angstroms\n", i, 1e8 * CLIGHT / numax_vspec_input[i], 1e8 * CLIGHT / numin_vspec_input[i]);
    }
  }
  else
  {
    Nrange = 1;

    numin_vspec_input[0] = numin_vspec;
    numax_vspec_input[0] = numax_vspec;

    printout("vpkt.txt: Nrange 1 frequency interval (inherited from numin_vspec and numax_vspec)\n");
    const int i = 0;
    printout("vpkt.txt:   range %d lambda [%g, %g] Angstroms\n", i, 1e8 * CLIGHT / numax_vspec_input[i], 1e8 * CLIGHT / numin_vspec_input[i]);
  }

  // if dum7=1, vpkt are not created when cell optical depth is larger than cell_is_optically_thick_vpkt
  int overrride_thickcell_tau = 0;
  fscanf(input_file, "%d %lg \n", &overrride_thickcell_tau, &cell_is_optically_thick_vpkt);

  if (overrride_thickcell_tau == 1)
  {
    printout("vpkt.txt: cell_is_optically_thick_vpkt %lg\n", cell_is_optically_thick_vpkt);
  }
  else
  {
    cell_is_optically_thick_vpkt = globals::cell_is_optically_thick;
    printout("vpkt.txt: cell_is_optically_thick_vpkt %lg (inherited from cell_is_optically_thick)\n", cell_is_optically_thick_vpkt);
  }

  // Maximum optical depth. If a vpkt reaches dum7 is thrown away
  fscanf(input_file, "%lg \n", &tau_max_vpkt);
  printout("vpkt.txt: tau_max_vpkt %g\n", tau_max_vpkt);

  // Produce velocity grid map if =1
  fscanf(input_file, "%d \n", &vgrid_flag);
  printout("vpkt.txt: velocity grid map %s\n", (vgrid_flag == 1) ? "ENABLED" : "DISABLED");

  if (vgrid_flag == 1)
  {
    double tmin_grid_in_days;
    double tmax_grid_in_days;
    // Specify time range for velocity grid map
    fscanf(input_file, "%lg %lg \n", &tmin_grid_in_days, &tmax_grid_in_days);
    tmin_grid = tmin_grid_in_days * DAY;
    tmax_grid = tmax_grid_in_days * DAY;

    printout("vpkt.txt: velocity grid time range tmin_grid %gd tmax_grid %gd\n", tmin_grid / DAY, tmax_grid / DAY);

    // Specify wavelength range: number of intervals (dum9) and limits (dum10,dum11)
    fscanf(input_file, "%d ", &Nrange_grid);

    printout("vpkt.txt: velocity grid frequency intervals %d\n", Nrange_grid);

    assert_always(Nrange_grid <= MRANGE_GRID);

    for (int i = 0; i < Nrange_grid; i++)
    {
      double range_lambda_min = 0.;
      double range_lambda_max = 0.;
      fscanf(input_file, "%lg %lg", &range_lambda_min, &range_lambda_max);

      nu_grid_max[i] = CLIGHT / (range_lambda_min * 1e-8);
      nu_grid_min[i] = CLIGHT / (range_lambda_max * 1e-8);

      printout("vpkt.txt:   velgrid range %d lambda [%g, %g] Angstroms\n", i, 1e8 * CLIGHT / nu_grid_max[i], 1e8 * CLIGHT / nu_grid_min[i]);
    }
  }

  fclose(input_file);
}


__host__ __device__
int vpkt_call_estimators(PKT *pkt_ptr, double t_current, int realtype)
{
  double obs[3];
  int vflag = 0;

  double vel_vec[3];
  get_velocity(pkt_ptr->pos, vel_vec, t_current);

  // Cut on vpkts
  int mgi = grid::get_cell_modelgridindex(pkt_ptr->where);

  if (grid::modelgrid[mgi].thick != 0)
  {
    return 0;
  }

  /* this is just to find the next_trans value when is set to 0 (avoid doing that in the vpkt routine for each observer) */
  if (pkt_ptr->next_trans == 0)
  {
    const int lineindex = closest_transition(pkt_ptr->nu_cmf, pkt_ptr->next_trans);  ///returns negative
    if (lineindex < 0)
    {
      pkt_ptr->next_trans = lineindex + 1;
    }
  }

  for (int bin = 0; bin < Nobs; bin++)
  {
    /* loop over different observers */

    obs[0] = sqrt(1 - nz_obs_vpkt[bin] * nz_obs_vpkt[bin]) * cos(phiobs[bin]);
    obs[1] = sqrt(1 - nz_obs_vpkt[bin] * nz_obs_vpkt[bin]) * sin(phiobs[bin]);
    obs[2] = nz_obs_vpkt[bin];

    const double t_arrive = t_current - (dot(pkt_ptr->pos, obs) / globals::CLIGHT_PROP);

    if (t_arrive >= tmin_vspec_input && t_arrive <= tmax_vspec_input)
    {
      // time selection

      for (int i = 0; i < Nrange; i++)
      {
        // Loop over frequency ranges

        if (pkt_ptr->nu_cmf / doppler(obs, vel_vec) > numin_vspec_input[i] && pkt_ptr->nu_cmf / doppler(obs, vel_vec) < numax_vspec_input[i])
        {
          // frequency selection

          rlc_emiss_vpkt(pkt_ptr, t_current, bin, obs, realtype);

          vflag = 1;

          // Need to update the starting cell for next observer
          // If previous vpkt reached tau_lim, change_cell (and then update_cell) hasn't been called
          mgi = grid::get_cell_modelgridindex(pkt_ptr->where);
          cellhistory_reset(mgi, false);
        }
      }
    }
  }

  // we just used the opacity variables for v-packets. We need to reset them for the original r packet
  calculate_kappa_rpkt_cont(pkt_ptr, &globals::kappa_rpkt_cont[tid]);

  return vflag;
}



double rot_angle(double *n1, double *n2, double *ref1, double *ref2)
{
/* ------------- Rotation angle from the scattering plane --------------------------------------------- */
/* -------- We need to rotate Stokes Parameters to (or from) the scattering plane from (or to) -------- */
/* -------- the meridian frame such that Q=1 is in the scattering plane and along ref1 ---------------- */

  double ref1_sc[3], i;

  // ref1_sc is the ref1 axis in the scattering plane ref1 = n1 x ( n1 x n2 )
  ref1_sc[0] = n1[0] * dot(n1,n2) - n2[0];
  ref1_sc[1] = n1[1] * dot(n1,n2) - n2[1];
  ref1_sc[2] = n1[2] * dot(n1,n2) - n2[2];
  vec_norm(ref1_sc, ref1_sc);

  double cos_stokes_rot_1 = dot(ref1_sc,ref1);
  double cos_stokes_rot_2 = dot(ref1_sc,ref2);

  if (cos_stokes_rot_1 < -1) cos_stokes_rot_1 = -1;
  if (cos_stokes_rot_1 > 1) cos_stokes_rot_1 = 1;

  if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 > 0)) i = acos(cos_stokes_rot_1);
  if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 < 0)) i = 2 * acos(-1.) - acos(cos_stokes_rot_1);
  if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 < 0)) i = acos(-1.) + acos(fabs(cos_stokes_rot_1));
  if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 > 0)) i = acos(-1.) - acos(fabs(cos_stokes_rot_1));
  if (cos_stokes_rot_1 == 0) i = acos(-1.) / 2.;
  if (cos_stokes_rot_2 == 0) i = 0.0 ;

  if (!std::isfinite(i))
    printout("Warning NaN: %3.6f \t %3.6f \t %3.6f \n", cos_stokes_rot_1, cos_stokes_rot_2, acos(cos_stokes_rot_1));

  return i;
}



// Routine to compute the meridian frame axes ref1 and ref2
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


// Routine to transform the Stokes Parameters from RF to CMF
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
  rot_angle = 0;

  if (p > 0)
  {
    cos2rot_angle = Q0/p;
    sin2rot_angle = U0/p;

    if ((cos2rot_angle > 0) && (sin2rot_angle > 0)) rot_angle = acos(Q0 / p) / 2.;
    if ((cos2rot_angle < 0) && (sin2rot_angle > 0)) rot_angle = (acos(-1.) - acos(fabs(Q0 / p))) / 2.;
    if ((cos2rot_angle < 0) && (sin2rot_angle < 0)) rot_angle = (acos(-1.) + acos(fabs(Q0 / p))) / 2.;
    if ((cos2rot_angle > 0) && (sin2rot_angle < 0)) rot_angle = (2. * acos(-1.) - acos(fabs(Q0 / p))) / 2.;
    if (cos2rot_angle == 0)
    {
      rot_angle = 0.25 * acos(-1);
      if (U0 < 0) rot_angle = 0.75 * acos(-1);
    }
    if (sin2rot_angle == 0)
    {
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
  if (e_cmf_ref1 == 0) theta_rot = acos(-1.)/2.;
  if (e_cmf_ref2 == 0) theta_rot = 0.0;
  if (e_cmf_ref1 > 1) theta_rot = 0.0;
  if (e_cmf_ref1 < -1) theta_rot = acos(-1.);

  // Compute Stokes Parameters in the CMF
  *Q = cos(2 * theta_rot ) * p;
  *U = sin(2 * theta_rot ) * p;
}


/* ----------------------- Lorentz transformations from RF to CMF --------------------------------------------- */
void lorentz(double *e_rf, double *n_rf, double *v, double *e_cmf)
{
  double beta[3], e_par[3], e_perp[3], b_rf[3], b_par[3], b_perp[3], vsqr, gamma_rel, v_cr_b[3], v_cr_e[3], b_cmf[3];

  beta[0] = v[0] / CLIGHT;
  beta[1] = v[1] / CLIGHT;
  beta[2] = v[2] / CLIGHT;
  vsqr = dot(beta, beta);

  gamma_rel = 1. / (sqrt(1 - vsqr));

  e_par[0] = (e_rf[0]*beta[0] + e_rf[1]*beta[1] + e_rf[2]*beta[2]) * beta[0] / (vsqr);
  e_par[1] = (e_rf[0]*beta[0] + e_rf[1]*beta[1] + e_rf[2]*beta[2]) * beta[1] / (vsqr);
  e_par[2] = (e_rf[0]*beta[0] + e_rf[1]*beta[1] + e_rf[2]*beta[2]) * beta[2] / (vsqr);

  e_perp[0] = e_rf[0] - e_par[0];
  e_perp[1] = e_rf[1] - e_par[1];
  e_perp[2] = e_rf[2] - e_par[2];

  b_rf[0] = n_rf[1]*e_rf[2] - n_rf[2]*e_rf[1];
  b_rf[1] = n_rf[2]*e_rf[0] - n_rf[0]*e_rf[2];
  b_rf[2] = n_rf[0]*e_rf[1] - n_rf[1]*e_rf[0];

  b_par[0] = (b_rf[0]*beta[0] + b_rf[1]*beta[1] + b_rf[2]*beta[2]) * beta[0] / (vsqr);
  b_par[1] = (b_rf[0]*beta[0] + b_rf[1]*beta[1] + b_rf[2]*beta[2]) * beta[1] / (vsqr);
  b_par[2] = (b_rf[0]*beta[0] + b_rf[1]*beta[1] + b_rf[2]*beta[2]) * beta[2] / (vsqr);

  b_perp[0] = b_rf[0] - b_par[0];
  b_perp[1] = b_rf[1] - b_par[1];
  b_perp[2] = b_rf[2] - b_par[2];

  v_cr_b[0] = beta[1]*b_rf[2] - beta[2]*b_rf[1];
  v_cr_b[1] = beta[2]*b_rf[0] - beta[0]*b_rf[2];
  v_cr_b[2] = beta[0]*b_rf[1] - beta[1]*b_rf[0];

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