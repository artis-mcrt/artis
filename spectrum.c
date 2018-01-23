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

// would do this in Python artistools except it currently doesn't understand how to match linelindex to a transition
#define TRACE_EMISSION_ABSORPTION_REGION_ON

#ifdef TRACE_EMISSION_ABSORPTION_REGION_ON
  #define traceemissabs_lambdamin 1000.  // in Angstroms
  #define traceemissabs_lambdamax 10000.
  #define traceemissabs_nulower (1.e8 * CLIGHT / traceemissabs_lambdamax)
  #define traceemissabs_nuupper (1.e8 * CLIGHT / traceemissabs_lambdamin)
  #define traceemissabs_timemin 300.
  #define traceemissabs_timemax 340.

  typedef struct emissionabsorptioncontrib
  {
    double energyemitted;
    double emission_weightedvelocity_sum;
    double energyabsorbed;
    double absorption_weightedvelocity_sum;
    int lineindex;   // this will be important when the list gets sorted
  } emissionabsorptioncontrib;

  struct emissionabsorptioncontrib *traceemissionabsorption;
  double traceemission_totalenergy = 0.;
  double traceabsorption_totalenergy = 0.;
  static int compare_emission(const void *p1, const void *p2)
  {
    const struct emissionabsorptioncontrib *elem1 = p1;
    const struct emissionabsorptioncontrib *elem2 = p2;

   if (elem1->energyemitted < elem2->energyemitted)
      return 1;
   else if (elem1->energyemitted > elem2->energyemitted)
      return -1;
   else
      return 0;
  }

  static int compare_absorption(const void *p1, const void *p2)
  {
    const struct emissionabsorptioncontrib *elem1 = p1;
    const struct emissionabsorptioncontrib *elem2 = p2;

   if (elem1->energyabsorbed < elem2->energyabsorbed)
      return 1;
   else if (elem1->energyabsorbed > elem2->energyabsorbed)
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
    spec_file = fopen_required("spec.out", "r");
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


  emission_file = fopen_required("emission.out", "w");
  /// write to file
  spec_file = fopen_required("spec.out", "w+");
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

      #ifdef TRACE_EMISSION_ABSORPTION_REGION_ON
      const int et = pkt_ptr->trueemissiontype;
      if (et >= 0)
      {
        if (nt >= 34 && nt <= 57)
        // if (t_arrive >= traceemissabs_timemin && t_arrive <= traceemissabs_timemax)
        {
          if (pkt_ptr->nu_rf >= traceemissabs_nulower && pkt_ptr->nu_rf <= traceemissabs_nuupper)
          {
            traceemissionabsorption[et].energyemitted += deltaE;

            traceemissionabsorption[et].emission_weightedvelocity_sum += pkt_ptr->trueemissionvelocity * deltaE;

            traceemission_totalenergy += deltaE;

            // const int element = linelist[et].elementindex;
            // const int ion = linelist[et].ionindex;
            // printout("packet in range, Z=%d ion_stage %d upperlevel %4d lowerlevel %4d fluxcontrib %g linecontrib %g index %d nlines %d\n",
            //          get_element(element), get_ionstage(element, ion), linelist[et].upperlevelindex, linelist[et].lowerlevelindex, deltaE, traceemissionabsorption[et].energyemitted, et, nlines);
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

          #ifdef TRACE_EMISSION_ABSORPTION_REGION_ON
          if (pkt_ptr->nu_rf >= traceemissabs_nulower && pkt_ptr->nu_rf <= traceemissabs_nuupper)
          {
            traceemissionabsorption[at].energyabsorbed += deltaE_absorption;

            double vel_vec[3];
            get_velocity(pkt_ptr->em_pos, vel_vec, pkt_ptr->em_time);
            traceemissionabsorption[at].absorption_weightedvelocity_sum += vec_len(vel_vec) * deltaE_absorption;

            traceabsorption_totalenergy += deltaE_absorption;
          }
          #endif
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

  #ifdef TRACE_EMISSION_ABSORPTION_REGION_ON
  traceemission_totalenergy = 0.;
  traceemissionabsorption = malloc(nlines * sizeof(emissionabsorptioncontrib));
  traceabsorption_totalenergy = 0.;
  for (int i = 0; i < nlines; i++)
  {
    traceemissionabsorption[i].energyemitted = 0.;
    traceemissionabsorption[i].emission_weightedvelocity_sum = 0.;
    traceemissionabsorption[i].energyabsorbed = 0.;
    traceemissionabsorption[i].absorption_weightedvelocity_sum = 0.;
    traceemissionabsorption[i].lineindex = i; // this will be important when the list gets sorted
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

  #ifdef TRACE_EMISSION_ABSORPTION_REGION_ON
  const int maxlinesprinted = 100;

  for (int mode = 0; mode < 2; mode++) // mode is 0 for emission and 1 for absorption
  {
    if (mode == 0)
    {
      qsort(traceemissionabsorption, nlines, sizeof(emissionabsorptioncontrib), compare_emission);
      printout("lambda [%5.1f, %5.1f] nu %g %g\n",
               traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_nulower, traceemissabs_nuupper);

      printout("Top line emission contributions in the range lambda [%5.1f, %5.1f] time [%5.1fd, %5.1fd] (%g erg)\n",
               traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_timemin, traceemissabs_timemax,
               traceemission_totalenergy);
    }
    else
    {
      qsort(traceemissionabsorption, nlines, sizeof(emissionabsorptioncontrib), compare_absorption);
      printout("Top line absorption contributions in the range lambda [%5.1f, %5.1f] time [%5.1fd, %5.1fd] (%g erg)\n",
               traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_timemin, traceemissabs_timemax,
               traceabsorption_totalenergy);
    }

    // display the top entries of the sorted list
    int nlines_limited = nlines;
    if (nlines > maxlinesprinted)
      nlines_limited = maxlinesprinted;
    printout("%17s %4s %9s %5s %5s %8s %8s %4s %7s %7s %7s %7s\n", "energy", "Z", "ion_stage", "upper", "lower", "coll_str", "A", "forb", "lambda", "<v_rad>", "B_lu", "B_ul");
    for (int i = 0; i < nlines_limited; i++)
    {
      double encontrib;
      double totalenergy;
      if (mode == 0)
      {
        encontrib = traceemissionabsorption[i].energyemitted;
        totalenergy = traceemission_totalenergy;
      }
      else
      {
        encontrib = traceemissionabsorption[i].energyabsorbed;
        totalenergy = traceabsorption_totalenergy;
      }
      if (encontrib > 0.) // lines that emit/absorb some energy
      {
        const int lineindex = traceemissionabsorption[i].lineindex;
        const int element = linelist[lineindex].elementindex;
        const int ion = linelist[lineindex].ionindex;
        const double linelambda = 1e8 * CLIGHT / linelist[lineindex].nu;
        // flux-weighted average radial velocity of emission in km/s
        double v_rad;
        if (mode == 0)
          v_rad = traceemissionabsorption[i].emission_weightedvelocity_sum / traceemissionabsorption[i].energyemitted / 1e5;
        else
          v_rad = traceemissionabsorption[i].absorption_weightedvelocity_sum / traceemissionabsorption[i].energyabsorbed / 1e5;

        const int lower = linelist[lineindex].lowerlevelindex;
        const int upper = linelist[lineindex].upperlevelindex;

        const double statweight_target = statw_upper(lineindex);
        const double statweight_lower = statw_lower(lineindex);

        const double nu_trans = (epsilon(element, ion, upper) - epsilon(element, ion, lower)) / H;
        const double A_ul = einstein_spontaneous_emission(lineindex);
        const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
        const double B_lu = statweight_target / statweight_lower * B_ul;

        // const double n_l = get_levelpop(modelgridindex,element,ion,lower);
        // const double n_u = get_levelpop(modelgridindex,element,ion,upper);
        // const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * em_time;

        printout("%7.2e (%5.1f%%) %4d %9d %5d %5d %8.1f %8.2e %4d %7.1f %7.1f %7.1e %7.1e\n",
                 encontrib, 100 * encontrib / totalenergy, get_element(element),
                 get_ionstage(element, ion), linelist[lineindex].upperlevelindex, linelist[lineindex].lowerlevelindex,
                 linelist[lineindex].coll_str, einstein_spontaneous_emission(lineindex), linelist[lineindex].forbidden,
                 linelambda, v_rad, B_lu, B_ul);
       }
       else
        break;
    }
    printout("\n");
  }

  free(traceemissionabsorption);
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
