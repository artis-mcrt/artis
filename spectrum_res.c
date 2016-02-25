#include "sn3d.h"
#include "exspec.h"

/* Routine to make a MC spectrum from the packets. */

int make_spectrum_res()
{
  int gather_spectrum_res();
  int write_spectrum_res();

  gather_spectrum_res(0);
  write_spectrum_res();
  return 0;
}


void init_spectrum_res()
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

  /* start by setting up the time and frequency bins. */
  /* it is all done interms of a logarithmic spacing in both t and nu - get the
     step sizes first. */
  ///Should be moved to input.c or exspec.c
  dlogt = (log(tmax) - log(tmin))/ntbins;
  dlognu = (log(nu_max_r) - log(nu_min_r))/nnubins;

  for (int nn = 0; nn < MABINS; nn++)
  {
    for (int n = 0; n < ntbins; n++)
    {
      spectra_res[n][nn].lower_time = exp( log(tmin) + (n * (dlogt)));
      spectra_res[n][nn].delta_t = exp( log(tmin) + ((n+1) * (dlogt))) - spectra_res[n][nn].lower_time;
      for (int m = 0; m < nnubins; m++)
      {
        spectra_res[n][nn].lower_freq[m] = exp( log(nu_min_r) + (m * (dlognu)));
        spectra_res[n][nn].delta_freq[m] = exp( log(nu_min_r) + ((m+1) * (dlognu))) - spectra_res[n][nn].lower_freq[m];
        spectra_res[n][nn].flux[m] = 0.0;

        if (do_emission_res == 1)
        {
          for (int i = 0; i < 2*nelements*maxion+1; i++)
            spectra_res[n][nn].emission[m].count[i] = 0;
        }
      }
    }
  }
}

/***********************************************/
int gather_spectrum_res(int my_rank)
{
  int add_to_spec_res();

  /// The grids are now set up. Now we loop over all the packets, check if they made it out or not,
  /// and if they did we add their rest frame energy to the appropriate cell.

  /// And figure out the escaping packets of those.
  for (int p = 0; p < npkts; p++)
  {
    PKT *pkt_ptr;
    pkt_ptr = &pkt[p];
    if (pkt_ptr->type == TYPE_ESCAPE && pkt_ptr->escape_type == TYPE_RPKT)
    {
      /// It made it out.
      add_to_spec_res(pkt_ptr);
    }
  }

  return 0;
}


/*************************************************************************/
int write_spectrum_res()
{
  FILE *spec_file;
  FILE *emission_file;
  float dum1, dum2;

  /// The spectra are now done - just need to print them out.
  if (file_set == 1)
  {
    if ((spec_file = fopen("spec_res.out", "r")) == NULL)
    {
      printout("Cannot open spec_file.txt.\n");
      exit(0);
    }
    fscanf(spec_file, "%g ", &dum1);

    for (int p = 0; p < ntbins; p++)
    {
      fscanf(spec_file, "%g ", &dum1);
    }

    for (int nn = 0; nn < MABINS; nn++)
    {
      for (int m = 0; m < nnubins; m++)
      {
        fscanf(spec_file, "%g ", &dum1);

        for (int p = 0; p < ntbins; p++)
        {
          fscanf(spec_file, "%g ", &dum2);
          spectra_res[p][nn].flux[m] += dum2;
        }
      }
    }
    fclose(spec_file);
  }

  if ((spec_file = fopen("spec_res.out", "w+")) == NULL)
  {
    printout("Cannot open spec_file.txt.\n");
    exit(0);
  }

  if (do_emission_res == 1)
  {
    if ((emission_file = fopen("emission_res.out", "w")) == NULL)
    {
      printf("Cannot open emission file\n");
      exit(0);
    }
  }

  fprintf(spec_file, "%g ", 0.0);

  for (int p = 0; p < ntbins; p++)
  {
    fprintf(spec_file, "%g ", (spectra_res[p][0].lower_time + (spectra_res[0][p].delta_t/2))/DAY);
  }

  fprintf(spec_file, "\n");

  for (int nn = 0; nn < MABINS; nn++)
  {
    for (int m = 0; m < nnubins; m++)
    {
      fprintf(spec_file, "%g ", ((spectra_res[0][0].lower_freq[m]+(spectra_res[0][0].delta_freq[m]/2))));

      for (int p = 0; p < ntbins; p++)
      {
        fprintf(spec_file, "%g ", spectra_res[p][nn].flux[m]);

        if (do_emission_res == 1)
        {
          for (int i = 0; i < 2*nelements*maxion+1; i++)
            fprintf(emission_file, "%d ", spectra_res[p][nn].emission[m].count[i]);
          fprintf(emission_file, "\n");
        }
      }
    fprintf(spec_file, "\n");
    }

  }

  fclose(spec_file);
  if (do_emission_res == 1)
    fclose(emission_file);

  return 0;
}

/**********************************************************************/
/*Routine to add a packet to the outcoming spectrum.*/
int add_to_spec_res(PKT *pkt_ptr)
{
  /* Need to (1) decide which time bin to put it in and (2) which frequency bin. */

  /* Time bin - we know that it escaped at "escape_time". However, we have to allow
     for travel time. Use the formula in Leon's paper.
     The extra distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
  */

  int i, nt, nnu, na;
  int thetabin, phibin;
  double vec1[3], vec2[3], vec3[3], xhat[3];
  double costheta, cosphi, testphi;
  int et,element,ion,nproc;

  xhat[0]=1.0;
  xhat[1]=0;
  xhat[2]=0;

  double t_arrive = pkt_ptr->escape_time - (dot(pkt_ptr->pos, pkt_ptr->dir)/CLIGHT_PROP);

  /* Put this into the time grid. */

  if (t_arrive > tmin && t_arrive < tmax)
  {
    nt = (log(t_arrive) - log(tmin)) / dlogt;

    /* for angle resolved case, need to work out the correct angle bin */

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

    if (pkt_ptr->nu_rf > nu_min_r && pkt_ptr->nu_rf < nu_max_r)
    {
      nnu = (log(pkt_ptr->nu_rf) - log(nu_min_r)) /  dlognu;
      spectra_res[nt][na].flux[nnu] += pkt_ptr->e_rf / spectra_res[nt][na].delta_t / spectra_res[nt][na].delta_freq[nnu] / 4.e12 / PI / PARSEC /PARSEC * MABINS / nprocs;

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

        spectra_res[nt][na].emission[nnu].count[nproc] += 1;
      }
    }
  }

  return 0;
}
