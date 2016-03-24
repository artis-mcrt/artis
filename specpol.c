#include "sn3d.h"
#include "exspec.h"
#include "specpol.h"
#include "vectors.h"

/*******************************************************/
int write_specpol(FILE *specpol_file, FILE *emissionpol_file, FILE *absorptionpol_file)
{
  int i,m,p,l;

  fprintf(specpol_file, "%g ", 0.0);

  for (l=0;l<3;l++) {
    for (p = 0; p < ntbins; p++) fprintf(specpol_file, "%g ", (stokes_i[p].lower_time + (stokes_i[p].delta_t/2))/DAY);
  }

  fprintf(specpol_file, "\n");

  for (m=0; m < nnubins; m++)
  {

    fprintf(specpol_file, "%g ", ((stokes_i[0].lower_freq[m]+(stokes_i[0].delta_freq[m]/2))));

    // Stokes I
    for (p = 0; p < ntbins; p++)
    {
      fprintf(specpol_file, "%g ", stokes_i[p].flux[m]);

      if (do_emission_res == 1)
      {
        for (i = 0; i < 2*nelements*maxion+1; i++) fprintf(emissionpol_file, "%g ", stokes_i[p].stat[m].emission[i]);
        fprintf(emissionpol_file, "\n");
        for (i = 0; i < nelements*maxion; i++) fprintf(absorptionpol_file, "%g ", stokes_i[p].stat[m].absorption[i]);
        fprintf(absorptionpol_file, "\n");
      }
    }


    // Stokes Q
    for (p = 0; p < ntbins; p++)
    {
        fprintf(specpol_file, "%g ", stokes_q[p].flux[m]);

        if (do_emission_res == 1)
        {
            for (i = 0; i < 2*nelements*maxion+1; i++) fprintf(emissionpol_file, "%g ", stokes_q[p].stat[m].emission[i]);
            fprintf(emissionpol_file, "\n");
            for (i = 0; i < nelements*maxion; i++) fprintf(absorptionpol_file, "%g ", stokes_q[p].stat[m].absorption[i]);
            fprintf(absorptionpol_file, "\n");
        }
    }


    // Stokes U
    for (p = 0; p < ntbins; p++)
    {
        fprintf(specpol_file, "%g ", stokes_u[p].flux[m]);

        if (do_emission_res == 1)
        {
            for (i = 0; i < 2*nelements*maxion+1; i++) fprintf(emissionpol_file, "%g ", stokes_u[p].stat[m].emission[i]);
            fprintf(emissionpol_file, "\n");
            for (i = 0; i < nelements*maxion; i++) fprintf(absorptionpol_file, "%g ", stokes_u[p].stat[m].absorption[i]);
            fprintf(absorptionpol_file, "\n");
        }
    }

    fprintf(specpol_file, "\n");
  }

  /*
  fclose(specpol_file);
  fclose(emissionpol_file);
  */
  return(0);
}


/**********************************************************************/
void init_specpol(void)
{
  int n,m,i;

  if (nnubins > MNUBINS)
  {
    printout("Too many frequency bins in spectrum - reducing.\n");
    nnubins = MNUBINS;
  }
  if (ntbins > MTBINS)
  {
    printout("Too many time bins in spectrum - reducing.\n");
    ntbins = MTBINS;
  }

  /** start by setting up the time and frequency bins. */
  /** it is all done interms of a logarithmic spacing in both t and nu - get the
  step sizes first. */
  ///Should be moved to input.c or exspec.c
  dlogt = (log(tmax) - log(tmin))/ntbins;
  dlognu = (log(nu_max_r) - log(nu_min_r))/nnubins;

  for (n = 0; n < ntbins; n++)
  {
    stokes_i[n].lower_time = exp( log(tmin) + (n * (dlogt)));
    stokes_i[n].delta_t = exp( log(tmin) + ((n+1) * (dlogt))) - stokes_i[n].lower_time;
    for (m=0; m < nnubins; m++)
    {
      stokes_i[n].lower_freq[m] = exp( log(nu_min_r) + (m * (dlognu)));
      stokes_i[n].delta_freq[m] = exp( log(nu_min_r) + ((m+1) * (dlognu))) - stokes_i[n].lower_freq[m];

      stokes_i[n].flux[m] = 0.0;
      stokes_q[n].flux[m] = 0.0;
      stokes_u[n].flux[m] = 0.0;

      for (i = 0; i < 2*nelements*maxion+1; i++) {
        stokes_i[n].stat[m].emission[i] = 0;
        stokes_q[n].stat[m].emission[i] = 0;
        stokes_u[n].stat[m].emission[i] = 0;
      }

      for (i = 0; i < nelements*maxion; i++) {
        stokes_i[n].stat[m].absorption[i] = 0;
        stokes_q[n].stat[m].absorption[i] = 0;
        stokes_u[n].stat[m].absorption[i] = 0;
      }

    }
  }
}


int gather_specpol(int depth)
{
  void init_specpol();
  int add_to_spec(EPKT *pkt_ptr);
  //void read_packets(FILE *packets_file);
  int p;
  EPKT *pkt_ptr;
  //int i,n,m,p;
  //PKT *pkt_ptr;
  //int add_to_spec();
  double vcut,vem;

  /// Set up the spectrum grid and initialise the bins to zero.
  init_specpol();

  if (depth < 0)
  {
    /// Do not extract depth-dependent spectra
    /// Now add the energy of all the escaping packets to the
    /// appropriate bins.
    for (p = 0; p < nepkts; p++)
    {
      pkt_ptr = &epkts[p];
      add_to_specpol(pkt_ptr);

    }
  }
  else
  {
    /// Extract depth-dependent spectra
    /// Set velocity cut
    if (depth < 9) vcut=(depth+1.)*vmax/10.;
    else vcut=100*vmax;  /// Make sure that all escaping packets are taken for
                        /// depth=9 . For 2d and 3d models the corners of a
                        /// simulation box contain material with v>vmax.

    /// Now add the energy of all the escaping packets to the
    /// appropriate bins.
    for (p = 0; p < nepkts; p++)
    {
      pkt_ptr = &epkts[p];
      vem = sqrt(pow(pkt_ptr->em_pos[0],2)+pow(pkt_ptr->em_pos[1],2)+pow(pkt_ptr->em_pos[2],2))/pkt_ptr->em_time;
      //printout("vem %g, vcut %g, vmax %g, time %d\n",vem,vcut,vmax,pkt_ptr->em_time);
      if (vem < vcut) add_to_specpol(pkt_ptr);
    }
  }

  return(0);
}


/**********************************************************************/
/*Routine to add a packet to the outcoming spectrum.*/
int add_to_specpol(EPKT *pkt_ptr)
{
  /** Need to (1) decide which time bin to put it in and (2) which frequency bin. */


  double dot();
  double t_arrive,deltai,deltaq,deltau;
  int nt, nnu;
  int at,et,element,ion,nproc;


  /// Put this into the time grid.
  t_arrive = pkt_ptr->arrive_time;
  if (t_arrive > tmin && t_arrive < tmax)
  {
    nt = (log(t_arrive) - log(tmin)) / dlogt;
    if (pkt_ptr->nu_rf > nu_min_r && pkt_ptr->nu_rf < nu_max_r)
    {
      nnu = (log(pkt_ptr->nu_rf) - log(nu_min_r)) /  dlognu;
      deltai = pkt_ptr->stokes[0]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs;
      deltaq = pkt_ptr->stokes[1]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs;
      deltau = pkt_ptr->stokes[2]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs;

      stokes_i[nt].flux[nnu] += deltai;
      stokes_q[nt].flux[nnu] += deltaq;
      stokes_u[nt].flux[nnu] += deltau;

      et = pkt_ptr->emissiontype;
      if (et >= 0)
      {
        /// bb-emission
        element = linelist[et].elementindex;
        ion = linelist[et].ionindex;
        nproc = element*maxion+ion;
      }
      else if (et == -9999999)
      {
        /// ff-emission
        nproc = 2*nelements*maxion;
      }
      else
      {
        /// bf-emission
        et = -1*et - 1;
        element = bflist[et].elementindex;
        ion = bflist[et].ionindex;
        nproc = nelements*maxion + element*maxion+ion;
      }
      stokes_i[nt].stat[nnu].emission[nproc] += deltai;
      stokes_q[nt].stat[nnu].emission[nproc] += deltaq;
      stokes_u[nt].stat[nnu].emission[nproc] += deltau;

      nnu = (log(pkt_ptr->absorptionfreq) - log(nu_min_r)) /  dlognu;
      if (nnu >= 0 && nnu < MNUBINS)
      {
        deltai = pkt_ptr->stokes[0]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs;
        deltaq = pkt_ptr->stokes[1]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs;
        deltau = pkt_ptr->stokes[2]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs;

        at = pkt_ptr->absorptiontype;
        if (at >= 0)
        {
          /// bb-emission
          element = linelist[at].elementindex;
          ion = linelist[at].ionindex;
          nproc = element*maxion+ion;
          stokes_i[nt].stat[nnu].absorption[nproc] += deltai;
          stokes_q[nt].stat[nnu].absorption[nproc] += deltaq;
          stokes_u[nt].stat[nnu].absorption[nproc] += deltau;
        }
      }
    }
  }

  return(0);

}



/***********************************************/
int gather_specpol_res(int current_abin)
{
  //void read_packets(FILE *packets_file);
  //int i,n,m,nn,p;
  int add_to_specpol_res(EPKT *pkt_ptr, int current_abin);
  void init_specpol(void);
  EPKT *pkt_ptr;
  int p;


  /// Set up the spectrum grid and initialise the bins to zero.
  init_specpol();

  /// Now add the energy of all the escaping packets to the
  /// appropriate bins.
  for (p = 0; p < nepkts; p++)
  {
    pkt_ptr = &epkts[p];
    add_to_specpol_res(pkt_ptr,current_abin);
  }

  return(0);
}



/**********************************************************************/
/**Routine to add a packet to the outcoming spectrum.*/
int add_to_specpol_res(EPKT *pkt_ptr, int current_abin)
{
  /* Need to (1) decide which time bin to put it in and (2) which frequency bin. */

  /* Time bin - we know that it escaped at "escape_time". However, we have to allow
     for travel time. Use the formula in Leon's paper.
     The extra distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
  */

  double t_arrive,deltai,deltaq,deltau;
  int nt, nnu, na;
  int thetabin, phibin;
  double vec1[3], vec2[3], vec3[3], xhat[3];
  double costheta, cosphi, testphi;
  int at,et,element,ion,nproc;

  xhat[0]=1.0;
  xhat[1]=0;
  xhat[2]=0;



  /// Angle resolved case: need to work out the correct angle bin
  costheta = dot(pkt_ptr->dir, syn_dir);
  thetabin = ((costheta + 1.0) * sqrt(MABINS) / 2.0);
  cross_prod(pkt_ptr->dir, syn_dir, vec1);
  cross_prod(xhat, syn_dir, vec2);
  cosphi = dot(vec1,vec2)/vec_len(vec1)/vec_len(vec2);

  cross_prod(vec2, syn_dir, vec3);
  testphi = dot(vec1,vec3);

  if (testphi > 0)
  {
    phibin = (acos(cosphi) /2. / PI * sqrt(MABINS));
  }
  else
  {
    phibin = ((acos(cosphi) + PI) /2. / PI * sqrt(MABINS));
  }
  na = (thetabin*sqrt(MABINS)) + phibin;

  /// Add only packets which escape to the current angle bin
  if (na == current_abin)
  {
    /// Put this into the time grid.
    t_arrive = pkt_ptr->arrive_time;
    if (t_arrive > tmin && t_arrive < tmax)
    {
      nt = (log(t_arrive) - log(tmin)) / dlogt;

      /// and the correct frequency bin.
      if (pkt_ptr->nu_rf > nu_min_r && pkt_ptr->nu_rf < nu_max_r)
      {
        nnu = (log(pkt_ptr->nu_rf) - log(nu_min_r)) /  dlognu;
        deltai = pkt_ptr->stokes[0]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC * MABINS / nprocs;
        deltaq = pkt_ptr->stokes[1]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC * MABINS / nprocs;
        deltau = pkt_ptr->stokes[2]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC * MABINS / nprocs;
        stokes_i[nt].flux[nnu] += deltai;
        stokes_q[nt].flux[nnu] += deltaq;
        stokes_u[nt].flux[nnu] += deltau;

        if (do_emission_res == 1)
        {
          et = pkt_ptr->emissiontype;
          if (et >= 0)
          {
            /// bb-emission
            element = linelist[et].elementindex;
            ion = linelist[et].ionindex;
            nproc = element*maxion+ion;
          }
          else if (et == -9999999)
          {
            /// ff-emission
            nproc = 2*nelements*maxion;
          }
          else
          {
            /// bf-emission
            et = -1*et - 1;
            element = bflist[et].elementindex;
            ion = bflist[et].ionindex;
            nproc = nelements*maxion + element*maxion+ion;
          }
          stokes_i[nt].stat[nnu].emission[nproc] += deltai;
          stokes_q[nt].stat[nnu].emission[nproc] += deltaq;
          stokes_u[nt].stat[nnu].emission[nproc] += deltau;

          nnu = (log(pkt_ptr->absorptionfreq) - log(nu_min_r)) /  dlognu;
          if (nnu >= 0 && nnu < MNUBINS)
          {
            deltai = pkt_ptr->stokes[0]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC * MABINS / nprocs;
            deltaq = pkt_ptr->stokes[1]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC * MABINS / nprocs;
            deltau = pkt_ptr->stokes[2]*pkt_ptr->e_rf / stokes_i[nt].delta_t / stokes_i[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC * MABINS / nprocs;

            at = pkt_ptr->absorptiontype;
            if (at >= 0)
            {
              /// bb-emission
              element = linelist[at].elementindex;
              ion = linelist[at].ionindex;
              nproc = element*maxion+ion;
              stokes_i[nt].stat[nnu].absorption[nproc] += deltai;
              stokes_q[nt].stat[nnu].absorption[nproc] += deltaq;
              stokes_u[nt].stat[nnu].absorption[nproc] += deltau;
            }
          }
        }
      }
    }
  }

  return(0);
}
