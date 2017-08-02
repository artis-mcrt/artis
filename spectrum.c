#include "sn3d.h"
#include "exspec.h"
#include "atomic.h"
#include "spectrum.h"
#include "vectors.h"

/*int make_spectrum()
/// Routine to make a MC spectrum from the packets.
{
  gather_spectrum(0);
  write_spectrum();

  return 0;
}*/

#define TRACE_EMISSION_REGION_ON

#ifdef TRACE_EMISSION_REGION_ON
  #define traceemiss_nulower (CLIGHT / (5700e-8))  // in Angstroms
  #define traceemiss_nuupper (CLIGHT / (4200e-8))  // in Angstroms
  #define traceemiss_timestepmin 34
  #define traceemiss_timestepmax 57

  typedef struct emissioncontrib
  {
    double fluxcontrib;
    double emiss_velocity_sum;
    int lineindex;   // this will be important when the list gets sorted
  } emissioncontrib;

  struct emissioncontrib *traceemisscontributions;
  double traceemiss_totalflux = 0.;
  static int compare_emisscontrib(const void *p1, const void *p2)
  {
    const struct emissioncontrib *elem1 = p1;
    const struct emissioncontrib *elem2 = p2;

   if (elem1->fluxcontrib < elem2->fluxcontrib)
      return 1;
   else if (elem1->fluxcontrib > elem2->fluxcontrib)
      return -1;
   else
      return 0;
  }
#endif



void write_spectrum(FILE *spec_file, FILE *emission_file, FILE *trueemission_file, FILE *absorption_file)
{
  //FILE *spec_file,*emission_file;

  /*
  float dum1, dum2;
  /// The spectra are now done - just need to print them out.
  if (file_set)
  {
    if ((spec_file = fopen("spec.out", "r")) == NULL)
    {
      printout("Cannot open spec_file.txt.\n");
      abort();
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
    abort();
  }
  /// write to file
  if ((spec_file = fopen("spec.out", "w+")) == NULL){
    printout("Cannot open spec_file.txt.\n");
    abort();
  }
  */

  fprintf(spec_file, "%g ", 0.0);
  for (int p = 0; p < ntbins; p++)
  {
    /// ????????????????????????????????????????????????????????????????????????????????????????????????
    /// WHY HERE OTHER CALCULATION OF "SPECTRA.MIDTIME" THAN FOR time_step.mid ?????????????????????????
    fprintf(spec_file, "%g ", (spectra[p].lower_time + (spectra[p].delta_t/2)) / DAY);
  }
  fprintf(spec_file, "\n");

  for (int m = 0; m < nnubins; m++)
  {
    fprintf(spec_file, "%g ", ((spectra[0].lower_freq[m] + (spectra[0].delta_freq[m] / 2))));

    for (int p = 0; p < ntbins; p++)
    {
      fprintf(spec_file, "%g ", spectra[p].flux[m]);
      if (do_emission_res == 1)
      {
        for (int i = 0; i < 2 * nelements * maxion + 1; i++)
          fprintf(emission_file, "%g ", spectra[p].stat[m].emission[i]);
        fprintf(emission_file, "\n");

        for (int i = 0; i < 2 * nelements * maxion + 1; i++)
          fprintf(trueemission_file, "%g ", spectra[p].stat[m].trueemission[i]);
        fprintf(trueemission_file, "\n");

        for (int i = 0; i < nelements * maxion; i++)
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
}


static int columnindex_from_emissiontype(const int et)
{
  if (et >= 0)
  {
    /// bb-emission
    const int element = linelist[et].elementindex;
    const int ion = linelist[et].ionindex;
    return element * maxion + ion;
  }
  else if (et == -9999999)
  {
    /// ff-emission
    return 2 * nelements * maxion;
  }
  else
  {
    /// bf-emission
    const int et_new = -1 * et - 1;
    const int element = bflist[et_new].elementindex;
    const int ion = bflist[et_new].ionindex;
    return nelements * maxion + element * maxion + ion;
  }
}


static void add_to_spec(const EPKT *pkt_ptr)
// Routine to add a packet to the outgoing spectrum.
{
  // Need to (1) decide which time bin to put it in and (2) which frequency bin.

  /// Put this into the time grid.
  const double t_arrive = pkt_ptr->arrive_time;
  if (t_arrive > tmin && t_arrive < tmax)
  {
    const int nt = (log(t_arrive) - log(tmin)) / dlogt;
    if (pkt_ptr->nu_rf > nu_min_r && pkt_ptr->nu_rf < nu_max_r)
    {
      int nnu = (log(pkt_ptr->nu_rf) - log(nu_min_r)) /  dlognu;
      const double deltaE = pkt_ptr->e_rf / spectra[nt].delta_t / spectra[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC / PARSEC / nprocs;
      spectra[nt].flux[nnu] += deltaE;

      const int nproc = columnindex_from_emissiontype(pkt_ptr->emissiontype);
      spectra[nt].stat[nnu].emission[nproc] += deltaE;

      const int truenproc = columnindex_from_emissiontype(pkt_ptr->trueemissiontype);
      spectra[nt].stat[nnu].trueemission[truenproc] += deltaE;

      #ifdef TRACE_EMISSION_REGION_ON
      const int et = pkt_ptr->trueemissiontype;
      if (et >= 0)
      {
        if (nt >= traceemiss_timestepmin && nt <= traceemiss_timestepmax)
        {
          if (pkt_ptr->nu_rf >= traceemiss_nulower && pkt_ptr->nu_rf <= traceemiss_nuupper)
          {
            traceemisscontributions[et].fluxcontrib += deltaE;
            traceemiss_totalflux += deltaE;

            double vel_vec[3];
            get_velocity(pkt_ptr->em_pos, vel_vec, pkt_ptr->em_time);
            traceemisscontributions[et].emiss_velocity_sum += vec_len(vel_vec) * deltaE;

            // printout("packet in range, Z=%d ion_stage %d upperlevel %4d lowerlevel %4d fluxcontrib %g linecontrib %g index %d nlines %d\n",
                    //  get_element(element), get_ionstage(element, ion), linelist[et].upperlevelindex, linelist[et].lowerlevelindex, deltaE, traceemisscontributions[et].fluxcontrib, et, nlines);
          }
        }
      }
      #endif

      nnu = (log(pkt_ptr->absorptionfreq) - log(nu_min_r)) /  dlognu;
      if (nnu >= 0 && nnu < MNUBINS)
      {
        const double deltaE_absorption = pkt_ptr->e_rf / spectra[nt].delta_t / spectra[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC / PARSEC / nprocs;
        const int at = pkt_ptr->absorptiontype;
        if (at >= 0)
        {
          /// bb-emission
          const int element = linelist[at].elementindex;
          const int ion = linelist[at].ionindex;
          spectra[nt].stat[nnu].absorption[element * maxion + ion] += deltaE_absorption;
        }
      }
    }
  }
}


static void init_spectrum(void)
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

  // start by setting up the time and frequency bins.
  // it is all done interms of a logarithmic spacing in both t and nu - get the
  // step sizes first.
  ///Should be moved to input.c or exspec.c
  dlogt = (log(tmax) - log(tmin)) / ntbins;
  dlognu = (log(nu_max_r) - log(nu_min_r)) / nnubins;

  for (int n = 0; n < ntbins; n++)
  {
    spectra[n].lower_time = exp(log(tmin) + (n * (dlogt)));
    spectra[n].delta_t = exp(log(tmin) + ((n + 1) * (dlogt))) - spectra[n].lower_time;
    for (int m = 0; m < nnubins; m++)
    {
      spectra[n].lower_freq[m] = exp(log(nu_min_r) + (m * (dlognu)));
      spectra[n].delta_freq[m] = exp(log(nu_min_r) + ((m + 1) * (dlognu))) - spectra[n].lower_freq[m];
      spectra[n].flux[m] = 0.0;
      for (int i = 0; i < 2 * nelements * maxion + 1; i++)
      {
        spectra[n].stat[m].emission[i] = 0;
        spectra[n].stat[m].trueemission[i] = 0;
      }
      for (int i = 0; i < nelements * maxion; i++)
        spectra[n].stat[m].absorption[i] = 0;
    }
  }
}


void gather_spectrum(int depth)
{
  /// Set up the spectrum grid and initialise the bins to zero.
  init_spectrum();

  #ifdef TRACE_EMISSION_REGION_ON
  traceemiss_totalflux = 0.;
  traceemisscontributions = malloc(nlines*sizeof(emissioncontrib));
  for (int i = 0; i < nlines; i++)
  {
    traceemisscontributions[i].fluxcontrib = 0.;
    traceemisscontributions[i].emiss_velocity_sum = 0.;
    traceemisscontributions[i].lineindex = i; // this will be important when the list gets sorted
  }
  #endif

  if (depth < 0)
  {
    /// Do not extract depth-dependent spectra
    /// Now add the energy of all the escaping packets to the
    /// appropriate bins.
    for (int p = 0; p < nepkts; p++)
    {
      add_to_spec(&epkts[p]);
    }
  }
  else
  {
    /// Extract depth-dependent spectra
    /// Set velocity cut
    double vcut;
    if (depth < 9)
      vcut = (depth + 1.) * vmax / 10.;
    else
      vcut = 100 * vmax;     /// Make sure that all escaping packets are taken for
                             /// depth=9 . For 2d and 3d models the corners of a
                             /// simulation box contain material with v>vmax.

    /// Now add the energy of all the escaping packets to the
    /// appropriate bins.
    for (int p = 0; p < nepkts; p++)
    {
      const double vem = sqrt(pow(epkts[p].em_pos[0],2) + pow(epkts[p].em_pos[1],2) + pow(epkts[p].em_pos[2],2)) / (epkts[p].em_time);
      //printout("vem %g, vcut %g, vmax %g, time %d\n",vem,vcut,vmax,pkt_ptr->em_time);
      if (vem < vcut)
        add_to_spec(&epkts[p]);
    }
  }

  #ifdef TRACE_EMISSION_REGION_ON
  qsort(traceemisscontributions, nlines, sizeof(emissioncontrib), compare_emisscontrib);
  printout("Top line emission contributions in the range lambda [%5.1f, %5.1f] timestep [%d, %d] (flux %g)\n",
           1e8 * CLIGHT / traceemiss_nuupper, 1e8 * CLIGHT / traceemiss_nulower, traceemiss_timestepmin, traceemiss_timestepmax,
           traceemiss_totalflux);

  // display the top entries of the sorted list
  const int nlines_limited = nlines;
  if (nlines > 50)
    nlines = 50;
  printout("%17s %5s %9s %5s %5s %8s %8s %6s %6s %7s\n", "flux", "Z", "ion_stage", "upper", "lower", "coll_str", "A", "forbid", "lambda", "<v_rad>");
  for (int i = 0; i < nlines_limited; i++)
  {
    const double fluxcontrib = traceemisscontributions[i].fluxcontrib;
    if (fluxcontrib / traceemiss_totalflux > 0.01 || i < 5) // lines that contribute at least 5% of the flux, with a minimum of 5 lines and max of 50
    {
      const int lineindex = traceemisscontributions[i].lineindex;
      const int element = linelist[lineindex].elementindex;
      const int ion = linelist[lineindex].ionindex;
      const double linelambda = 1e8 * CLIGHT / linelist[lineindex].nu;
      // flux-weighted average radial velocity of emission in km/s
      const double v_rad = traceemisscontributions[i].emiss_velocity_sum / traceemisscontributions[i].fluxcontrib / 1e5;
      printout("%7.2e (%5.1f%%) %5d %9d %5d %5d %8.1f %8.2e %6d %6.1f %7.1f\n",
               fluxcontrib, 100 * fluxcontrib / traceemiss_totalflux, get_element(element),
               get_ionstage(element, ion), linelist[lineindex].upperlevelindex, linelist[lineindex].lowerlevelindex,
               linelist[lineindex].coll_str, einstein_spontaneous_emission(lineindex), linelist[lineindex].forbidden,
               linelambda, v_rad);
     }
     else
      break;
  }

  free(traceemisscontributions);
  #endif
}


static void add_to_spec_res(EPKT *pkt_ptr, int current_abin)
// Routine to add a packet to the outgoing spectrum.
{
  /* Need to (1) decide which time bin to put it in and (2) which frequency bin. */

  /* Time bin - we know that it escaped at "escape_time". However, we have to allow
     for travel time. Use the formula in Leon's paper.
     The extra distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
  */

  double xhat[3] = {1.0, 0.0, 0.0};

  /// Angle resolved case: need to work out the correct angle bin
  const double costheta = dot(pkt_ptr->dir, syn_dir);
  const double thetabin = ((costheta + 1.0) * sqrt(MABINS) / 2.0);
  double vec1[3];
  cross_prod(pkt_ptr->dir, syn_dir, vec1);
  double vec2[3];
  cross_prod(xhat, syn_dir, vec2);
  const double cosphi = dot(vec1,vec2) / vec_len(vec1) / vec_len(vec2);

  double vec3[3];
  cross_prod(vec2, syn_dir, vec3);
  const double testphi = dot(vec1,vec3);

  int phibin;
  if (testphi > 0)
  {
    phibin = (acos(cosphi) / 2. / PI * sqrt(MABINS));
  }
  else
  {
    phibin = ((acos(cosphi) + PI) / 2. / PI * sqrt(MABINS));
  }
  const int na = (thetabin*sqrt(MABINS)) + phibin;

  /// Add only packets which escape to the current angle bin
  if (na == current_abin)
  {
    /// Put this into the time grid.
    const double t_arrive = pkt_ptr->arrive_time;
    if (t_arrive > tmin && t_arrive < tmax)
    {
      const int nt = (log(t_arrive) - log(tmin)) / dlogt;

      /// and the correct frequency bin.
      if (pkt_ptr->nu_rf > nu_min_r && pkt_ptr->nu_rf < nu_max_r)
      {
        int nnu = (log(pkt_ptr->nu_rf) - log(nu_min_r)) /  dlognu;
        double deltaE = pkt_ptr->e_rf / spectra[nt].delta_t / spectra[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC / PARSEC * MABINS / nprocs;
        spectra[nt].flux[nnu] += deltaE;

        if (do_emission_res == 1)
        {
          const int nproc = columnindex_from_emissiontype(pkt_ptr->emissiontype);
          spectra[nt].stat[nnu].emission[nproc] += deltaE;

          const int truenproc = columnindex_from_emissiontype(pkt_ptr->trueemissiontype);
          spectra[nt].stat[nnu].trueemission[truenproc] += deltaE;

          nnu = (log(pkt_ptr->absorptionfreq) - log(nu_min_r)) /  dlognu;
          if (nnu >= 0 && nnu < MNUBINS)
          {
            deltaE = pkt_ptr->e_rf / spectra[nt].delta_t / spectra[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC / PARSEC * MABINS / nprocs;
            const int at = pkt_ptr->absorptiontype;
            if (at >= 0)
            {
              /// bb-emission
              const int element = linelist[at].elementindex;
              const int ion = linelist[at].ionindex;
              spectra[nt].stat[nnu].absorption[element * maxion + ion] += deltaE;
            }
          }
        }
      }
    }
  }
}


void gather_spectrum_res(const int current_abin)
{
  /// Set up the spectrum grid and initialise the bins to zero.
  init_spectrum();

  /// Now add the energy of all the escaping packets to the
  /// appropriate bins.
  for (int p = 0; p < nepkts; p++)
  {
    add_to_spec_res(&epkts[p], current_abin);
  }

}
