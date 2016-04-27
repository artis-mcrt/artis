#include "sn3d.h"
#include "exspec.h"
#include "spectrum.h"
#include "vectors.h"

// private functions
void init_spectrum(void);
int add_to_spec(const EPKT *pkt_ptr);
int add_to_spec_res(EPKT *pkt_ptr, int current_abin);


/*int make_spectrum()
/// Routine to make a MC spectrum from the packets.
{
  gather_spectrum(0);
  write_spectrum();

  return 0;
}*/


/*******************************************************/
int write_spectrum(FILE *spec_file, FILE *emission_file, FILE *absorption_file)
{
  //FILE *spec_file,*emission_file;

  /*
  float dum1, dum2;
  /// The spectra are now done - just need to print them out.
  if (file_set == 1)
  {
    if ((spec_file = fopen("spec.out", "r")) == NULL)
    {
      printout("Cannot open spec_file.txt.\n");
      exit(0);
    }
    fscanf(spec_file, "%g ", &dum1);

    for (p = 0; p < ntbins; p++)
    {
      fscanf(spec_file, "%g ", &dum1);
    }

    for (m=0; m < nnubins; m++)
    {
      fscanf(spec_file, "%g ", &dum1);

      for (p = 0; p < ntbins; p++)
      {
        fscanf(spec_file, "%g ",
                &dum2);
        spectra[p].flux[m] += dum2;
      }
    }
    fclose(spec_file);
  }


  if ((emission_file = fopen("emission.out", "w")) == NULL)
  {
    printf("Cannot open emission file\n");
    exit(0);
  }
  /// write to file
  if ((spec_file = fopen("spec.out", "w+")) == NULL){
    printout("Cannot open spec_file.txt.\n");
    exit(0);
  }
  */

  fprintf(spec_file, "%g ", 0.0);
  for (int p = 0; p < ntbins; p++)
  {
    /// ????????????????????????????????????????????????????????????????????????????????????????????????
    /// WHY HERE OTHER CALCULATION OF "SPECTRA.MIDTIME" THAN FOR time_step.mid ?????????????????????????
    fprintf(spec_file, "%g ", (spectra[p].lower_time + (spectra[p].delta_t/2))/DAY);
  }
  fprintf(spec_file, "\n");

  for (int m = 0; m < nnubins; m++)
  {
    fprintf(spec_file, "%g ", ((spectra[0].lower_freq[m]+(spectra[0].delta_freq[m]/2))));

    for (int p = 0; p < ntbins; p++)
    {
      fprintf(spec_file, "%g ", spectra[p].flux[m]);
      if (do_emission_res == 1)
      {
        for (int i = 0; i < 2*nelements*maxion+1; i++)
          fprintf(emission_file, "%g ", spectra[p].stat[m].emission[i]);
        fprintf(emission_file, "\n");
        for (int i = 0; i < nelements*maxion; i++)
          fprintf(absorption_file, "%g ", spectra[p].stat[m].absorption[i]);
        fprintf(absorption_file, "\n");
      }
    }
    fprintf(spec_file, "\n");
  }

  /*
  fclose(spec_file);
  fclose(emission_file);
  */
  return 0;
}


/**********************************************************************/
void init_spectrum(void)
{
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

  for (int n = 0; n < ntbins; n++)
  {
    spectra[n].lower_time = exp( log(tmin) + (n * (dlogt)));
    spectra[n].delta_t = exp( log(tmin) + ((n+1) * (dlogt))) - spectra[n].lower_time;
    for (int m = 0; m < nnubins; m++)
    {
      spectra[n].lower_freq[m] = exp( log(nu_min_r) + (m * (dlognu)));
      spectra[n].delta_freq[m] = exp( log(nu_min_r) + ((m+1) * (dlognu))) - spectra[n].lower_freq[m];
      spectra[n].flux[m] = 0.0;
      for (int i = 0; i < 2*nelements*maxion+1; i++)
        spectra[n].stat[m].emission[i] = 0;  ///added
      for (int i = 0; i < nelements*maxion; i++)
        spectra[n].stat[m].absorption[i] = 0;  ///added
    }
  }
}

int gather_spectrum(int depth)
{
  //void read_packets(FILE *packets_file);
  EPKT *pkt_ptr;
  //int i,n,m,p;
  //PKT *pkt_ptr;
  //int add_to_spec();

  /// Set up the spectrum grid and initialise the bins to zero.
  init_spectrum();

  if (depth < 0)
  {
    /// Do not extract depth-dependent spectra
    /// Now add the energy of all the escaping packets to the
    /// appropriate bins.
    for (int p = 0; p < nepkts; p++)
    {
      pkt_ptr = &epkts[p];
      add_to_spec(pkt_ptr);
    }
  }
  else
  {
    /// Extract depth-dependent spectra
    /// Set velocity cut
    double vcut;
    if (depth < 9)
      vcut = (depth+1.)*vmax/10.;
    else
      vcut = 100*vmax;     /// Make sure that all escaping packets are taken for
                        /// depth=9 . For 2d and 3d models the corners of a
                        /// simulation box contain material with v>vmax.

    /// Now add the energy of all the escaping packets to the
    /// appropriate bins.
    for (int p = 0; p < nepkts; p++)
    {
      pkt_ptr = &epkts[p];
      double vem = sqrt(pow(pkt_ptr->em_pos[0],2)+pow(pkt_ptr->em_pos[1],2)+pow(pkt_ptr->em_pos[2],2))/pkt_ptr->em_time;
      //printout("vem %g, vcut %g, vmax %g, time %d\n",vem,vcut,vmax,pkt_ptr->em_time);
      if (vem < vcut) add_to_spec(pkt_ptr);
    }
  }

  return 0;
}


/**********************************************************************/
/*Routine to add a packet to the outcoming spectrum.*/
int add_to_spec(const EPKT *pkt_ptr)
{
  /** Need to (1) decide which time bin to put it in and (2) which frequency bin. */
  int at,element,ion,nproc;

  /// Put this into the time grid.
  double t_arrive = pkt_ptr->arrive_time;
  if (t_arrive > tmin && t_arrive < tmax)
  {
    int nt = (log(t_arrive) - log(tmin)) / dlogt;
    if (pkt_ptr->nu_rf > nu_min_r && pkt_ptr->nu_rf < nu_max_r)
    {
      int nnu = (log(pkt_ptr->nu_rf) - log(nu_min_r)) /  dlognu;
      double deltaE = pkt_ptr->e_rf / spectra[nt].delta_t / spectra[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs;
      spectra[nt].flux[nnu] += deltaE;

      int et = pkt_ptr->emissiontype;
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
      spectra[nt].stat[nnu].emission[nproc] += deltaE;

      nnu = (log(pkt_ptr->absorptionfreq) - log(nu_min_r)) /  dlognu;
      if (nnu >= 0 && nnu < MNUBINS)
      {
        deltaE = pkt_ptr->e_rf / spectra[nt].delta_t / spectra[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC / nprocs;
        at = pkt_ptr->absorptiontype;
        if (at >= 0)
        {
          /// bb-emission
          element = linelist[at].elementindex;
          ion = linelist[at].ionindex;
          nproc = element*maxion+ion;
          spectra[nt].stat[nnu].absorption[nproc] += deltaE;
        }
      }
    }
  }

  return 0;
}



/***********************************************/
int gather_spectrum_res(int current_abin)
{
  //void read_packets(FILE *packets_file);
  //int i,n,m,nn,p;
  EPKT *pkt_ptr;

  /// Set up the spectrum grid and initialise the bins to zero.
  init_spectrum();

  /// Now add the energy of all the escaping packets to the
  /// appropriate bins.
  for (int p = 0; p < nepkts; p++)
  {
    pkt_ptr = &epkts[p];
    add_to_spec_res(pkt_ptr,current_abin);
  }

  return 0;
}


/**********************************************************************/
/**Routine to add a packet to the outcoming spectrum.*/
int add_to_spec_res(EPKT *pkt_ptr, int current_abin)
{
  /* Need to (1) decide which time bin to put it in and (2) which frequency bin. */

  /* Time bin - we know that it escaped at "escape_time". However, we have to allow
     for travel time. Use the formula in Leon's paper.
     The extra distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
  */

  double t_arrive,deltaE;
  int nt, nnu;
  int phibin;
  double vec1[3], vec2[3], vec3[3], xhat[3];
  int at,et,element,ion,nproc;

  xhat[0] = 1.0;
  xhat[1] = 0;
  xhat[2] = 0;

  /// Angle resolved case: need to work out the correct angle bin
  double costheta = dot(pkt_ptr->dir, syn_dir);
  double thetabin = ((costheta + 1.0) * sqrt(MABINS) / 2.0);
  cross_prod(pkt_ptr->dir, syn_dir, vec1);
  cross_prod(xhat, syn_dir, vec2);
  double cosphi = dot(vec1,vec2)/vec_len(vec1)/vec_len(vec2);

  cross_prod(vec2, syn_dir, vec3);
  double testphi = dot(vec1,vec3);

  if (testphi > 0)
  {
    phibin = (acos(cosphi) /2. / PI * sqrt(MABINS));
  }
  else
  {
    phibin = ((acos(cosphi) + PI) /2. / PI * sqrt(MABINS));
  }
  int na = (thetabin*sqrt(MABINS)) + phibin;

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
        deltaE = pkt_ptr->e_rf / spectra[nt].delta_t / spectra[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC * MABINS / nprocs;
        spectra[nt].flux[nnu] += deltaE;

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
          spectra[nt].stat[nnu].emission[nproc] += deltaE;

          nnu = (log(pkt_ptr->absorptionfreq) - log(nu_min_r)) /  dlognu;
          if (nnu >= 0 && nnu < MNUBINS)
          {
            deltaE = pkt_ptr->e_rf / spectra[nt].delta_t / spectra[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC * MABINS / nprocs;
            at = pkt_ptr->absorptiontype;
            if (at >= 0)
            {
              /// bb-emission
              element = linelist[at].elementindex;
              ion = linelist[at].ionindex;
              nproc = element*maxion+ion;
              spectra[nt].stat[nnu].absorption[nproc] += deltaE;
            }
          }
        }
      }
    }
  }

  return 0;
}
