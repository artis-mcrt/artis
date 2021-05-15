#include "sn3d.h"
#include "exspec.h"
#include "atomic.h"
#include "spectrum.h"
#include "vectors.h"


const bool TRACE_EMISSION_ABSORPTION_REGION_ON = true;

#define traceemissabs_lambdamin 1000.  // in Angstroms
#define traceemissabs_lambdamax 25000.
#define traceemissabs_nulower (1.e8 * CLIGHT / traceemissabs_lambdamax)
#define traceemissabs_nuupper (1.e8 * CLIGHT / traceemissabs_lambdamin)
#define traceemissabs_timemin 320. * DAY
#define traceemissabs_timemax 340. * DAY

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
  const struct emissionabsorptioncontrib *elem1 = (struct emissionabsorptioncontrib *) p1;
  const struct emissionabsorptioncontrib *elem2 = (struct emissionabsorptioncontrib *) p2;

 if (elem1->energyemitted < elem2->energyemitted)
    return 1;
 else if (elem1->energyemitted > elem2->energyemitted)
    return -1;
 else
    return 0;
}

static int compare_absorption(const void *p1, const void *p2)
{
  const struct emissionabsorptioncontrib *elem1 = (struct emissionabsorptioncontrib *) p1;
  const struct emissionabsorptioncontrib *elem2 = (struct emissionabsorptioncontrib *) p2;

 if (elem1->energyabsorbed < elem2->energyabsorbed)
    return 1;
 else if (elem1->energyabsorbed > elem2->energyabsorbed)
    return -1;
 else
    return 0;
}


static void printout_tracemission_stats(void)
{
  const int maxlinesprinted = 500;

  // mode is 0 for emission and 1 for absorption
  for (int mode = 0; mode < 2; mode++)
  {
    if (mode == 0)
    {
      qsort(traceemissionabsorption, globals::nlines, sizeof(emissionabsorptioncontrib), compare_emission);
      printout("lambda [%5.1f, %5.1f] nu %g %g\n",
               traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_nulower, traceemissabs_nuupper);

      printout("Top line emission contributions in the range lambda [%5.1f, %5.1f] time [%5.1fd, %5.1fd] (%g erg)\n",
               traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_timemin / DAY, traceemissabs_timemax / DAY,
               traceemission_totalenergy);
    }
    else
    {
      qsort(traceemissionabsorption, globals::nlines, sizeof(emissionabsorptioncontrib), compare_absorption);
      printout("Top line absorption contributions in the range lambda [%5.1f, %5.1f] time [%5.1fd, %5.1fd] (%g erg)\n",
               traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_timemin / DAY, traceemissabs_timemax / DAY,
               traceabsorption_totalenergy);
    }

    // display the top entries of the sorted list
    int nlines_limited = globals::nlines;
    if (globals::nlines > maxlinesprinted)
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
        const int element = globals::linelist[lineindex].elementindex;
        const int ion = globals::linelist[lineindex].ionindex;
        const double linelambda = 1e8 * CLIGHT / globals::linelist[lineindex].nu;
        // flux-weighted average radial velocity of emission in km/s
        double v_rad;
        if (mode == 0)
          v_rad = traceemissionabsorption[i].emission_weightedvelocity_sum / traceemissionabsorption[i].energyemitted / 1e5;
        else
          v_rad = traceemissionabsorption[i].absorption_weightedvelocity_sum / traceemissionabsorption[i].energyabsorbed / 1e5;

        const int lower = globals::linelist[lineindex].lowerlevelindex;
        const int upper = globals::linelist[lineindex].upperlevelindex;

        const double statweight_target = statw_upper(lineindex);
        const double statweight_lower = statw_lower(lineindex);

        const double nu_trans = (epsilon(element, ion, upper) - epsilon(element, ion, lower)) / H;
        const double A_ul = einstein_spontaneous_emission(lineindex);
        const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
        const double B_lu = statweight_target / statweight_lower * B_ul;

        // const double n_l = calculate_exclevelpop(modelgridindex,element,ion,lower);
        // const double n_u = calculate_exclevelpop(modelgridindex,element,ion,upper);
        // const double tau_sobolev = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * em_time;

        printout("%7.2e (%5.1f%%) %4d %9d %5d %5d %8.1f %8.2e %4d %7.1f %7.1f %7.1e %7.1e\n",
                 encontrib, 100 * encontrib / totalenergy, get_element(element),
                 get_ionstage(element, ion), globals::linelist[lineindex].upperlevelindex, globals::linelist[lineindex].lowerlevelindex,
                 globals::linelist[lineindex].coll_str, einstein_spontaneous_emission(lineindex), globals::linelist[lineindex].forbidden,
                 linelambda, v_rad, B_lu, B_ul);
       }
       else
        break;
    }
    printout("\n");
  }

  free(traceemissionabsorption);
}


void write_spectrum(
  char spec_filename[], char emission_filename[], char trueemission_filename[],
  char absorption_filename[], struct spec *spectra)
{
  FILE *spec_file = fopen_required(spec_filename, "w");

  FILE *emission_file = NULL;
  FILE *trueemission_file = NULL;
  FILE *absorption_file = NULL;

  bool do_emission_res = spectra[0].do_emission_res;

  if (do_emission_res)
  {
    emission_file = fopen_required(emission_filename, "w");
    assert_always(emission_file != NULL);
    trueemission_file = fopen_required(trueemission_filename, "w");
    assert_always(trueemission_file != NULL);
    absorption_file = fopen_required(absorption_filename, "w");
    assert_always(absorption_file != NULL);
    printout("write_spectrum %s, %s, %s, and %s\n", spec_filename, emission_filename, trueemission_filename, absorption_filename);
  }
  else
  {
    printout("write_spectrum %s\n", spec_filename);
  }

  if (TRACE_EMISSION_ABSORPTION_REGION_ON && do_emission_res)
  {
    printout_tracemission_stats();
  }

  fprintf(spec_file, "%g ", 0.0);
  for (int p = 0; p < globals::ntstep; p++)
  {
    /// ????????????????????????????????????????????????????????????????????????????????????????????????
    /// WHY HERE OTHER CALCULATION OF "SPECTRA.MIDTIME" THAN FOR time_step.mid ?????????????????????????
    fprintf(spec_file, "%g ", globals::time_step[p].mid / DAY);
  }
  fprintf(spec_file, "\n");

  for (int m = 0; m < globals::nnubins; m++)
  {
    fprintf(spec_file, "%g ", ((spectra[0].lower_freq[m] + (spectra[0].delta_freq[m] / 2))));

    for (int p = 0; p < globals::ntstep; p++)
    {
      fprintf(spec_file, "%g ", spectra[p].flux[m]);
      if (do_emission_res)
      {
        for (int i = 0; i < 2 * get_nelements() * get_max_nions() + 1; i++)
          fprintf(emission_file, "%g ", spectra[p].stat[m].emission[i]);
        fprintf(emission_file, "\n");

        for (int i = 0; i < 2 * get_nelements() * get_max_nions() + 1; i++)
          fprintf(trueemission_file, "%g ", spectra[p].stat[m].trueemission[i]);
        fprintf(trueemission_file, "\n");

        for (int i = 0; i < get_nelements() * get_max_nions(); i++)
          fprintf(absorption_file, "%g ", spectra[p].stat[m].absorption[i]);
        fprintf(absorption_file, "\n");
      }
    }
    fprintf(spec_file, "\n");
  }

  fclose(spec_file);
  if (do_emission_res)
  {
    fclose(emission_file);
    fclose(trueemission_file);
    fclose(absorption_file);
  }
}


void write_specpol(
  char spec_filename[], char emission_filename[], char absorption_filename[],
  struct spec *stokes_i, struct spec *stokes_q, struct spec *stokes_u)
{
  FILE *specpol_file = fopen_required(spec_filename, "w");
  FILE *emissionpol_file = NULL;
  FILE *absorptionpol_file = NULL;

  bool do_emission_res = stokes_i[0].do_emission_res;

  if (do_emission_res)
  {
    emissionpol_file = fopen_required(emission_filename, "w");
    absorptionpol_file = fopen_required(absorption_filename, "w");
  }

  fprintf(specpol_file, "%g ", 0.0);

  for (int l = 0; l < 3; l++)
  {
    for (int p = 0; p < globals::ntstep; p++)
    {
      fprintf(specpol_file, "%g ", globals::time_step[p].mid / DAY);
    }
  }

  fprintf(specpol_file, "\n");

  for (int m = 0; m < globals::nnubins; m++)
  {
    fprintf(specpol_file, "%g ", ((stokes_i[0].lower_freq[m] + (stokes_i[0].delta_freq[m] / 2))));

    // Stokes I
    for (int p = 0; p < globals::ntstep; p++)
    {
      fprintf(specpol_file, "%g ", stokes_i[p].flux[m]);

      if (do_emission_res)
      {
        for (int i = 0; i < 2 * get_nelements() * get_max_nions() + 1; i++)
        {
          fprintf(emissionpol_file, "%g ", stokes_i[p].stat[m].emission[i]);
        }
        fprintf(emissionpol_file, "\n");

        for (int i = 0; i < get_nelements() * get_max_nions(); i++)
        {
          fprintf(absorptionpol_file, "%g ", stokes_i[p].stat[m].absorption[i]);
        }
        fprintf(absorptionpol_file, "\n");
      }
    }

    // Stokes Q
    for (int p = 0; p < globals::ntstep; p++)
    {
      fprintf(specpol_file, "%g ", stokes_q[p].flux[m]);

      if (do_emission_res)
      {
          for (int i = 0; i < 2 * get_nelements() * get_max_nions() + 1; i++)
          {
            fprintf(emissionpol_file, "%g ", stokes_q[p].stat[m].emission[i]);
          }
          fprintf(emissionpol_file, "\n");

          for (int i = 0; i < get_nelements() * get_max_nions(); i++)
          {
            fprintf(absorptionpol_file, "%g ", stokes_q[p].stat[m].absorption[i]);
          }
          fprintf(absorptionpol_file, "\n");
      }
    }

    // Stokes U
    for (int p = 0; p < globals::ntstep; p++)
    {
      fprintf(specpol_file, "%g ", stokes_u[p].flux[m]);

      if (do_emission_res)
      {
        for (int i = 0; i < 2 * get_nelements() * get_max_nions() + 1; i++)
        {
          fprintf(emissionpol_file, "%g ", stokes_u[p].stat[m].emission[i]);
        }
        fprintf(emissionpol_file, "\n");

        for (int i = 0; i < get_nelements() * get_max_nions(); i++)
        {
          fprintf(absorptionpol_file, "%g ", stokes_u[p].stat[m].absorption[i]);
        }
        fprintf(absorptionpol_file, "\n");
      }
    }

    fprintf(specpol_file, "\n");
  }

  fclose(specpol_file);
  if (do_emission_res)
  {
    fclose(emissionpol_file);
    fclose(absorptionpol_file);
  }
}

static int columnindex_from_emissiontype(const int et)
{
  if (et >= 0)
  {
    /// bb-emission
    const int element = globals::linelist[et].elementindex;
    const int ion = globals::linelist[et].ionindex;
    return element * get_max_nions() + ion;
  }
  else if (et == -9999999)
  {
    /// ff-emission
    return 2 * get_nelements() * get_max_nions();
  }
  else
  {
    /// bf-emission
    const int et_new = -1 - et;
    const int element = globals::bflist[et_new].elementindex;
    const int ion = globals::bflist[et_new].ionindex;
    const int level = globals::bflist[et_new].levelindex;
    const int phixstargetindex = globals::bflist[et_new].phixstargetindex;
    const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);

    assert_always(get_continuumindex(element, ion, level, upperionlevel) == et);

    return get_nelements() * get_max_nions() + element * get_max_nions() + ion;
  }
}


static void add_to_spec(const PKT *const pkt_ptr, const int current_abin, struct spec *spectra, struct spec *stokes_i, struct spec *stokes_q, struct spec *stokes_u)
// Routine to add a packet to the outgoing spectrum.
{
  // Need to (1) decide which time bin to put it in and (2) which frequency bin.

  // specific angle bins contain fewer packets than the full sphere, so must be normalised to match
  const double anglefactor = (current_abin >= 0) ? MABINS : 1.;

  const double t_arrive = get_arrive_time(pkt_ptr);
  if (t_arrive > globals::tmin && t_arrive < globals::tmax && pkt_ptr->nu_rf > globals::nu_min_r && pkt_ptr->nu_rf < globals::nu_max_r)
  {
    const int nt = get_timestep(t_arrive);
    const double nu_min = spectra[nt].nu_min;
    const double nu_max = spectra[nt].nu_max;
    const double dlognu = (log(nu_max) - log(nu_min)) / globals::nnubins;

    const int nnu = (log(pkt_ptr->nu_rf) - log(nu_min)) /  dlognu;

    const double deltaE = pkt_ptr->e_rf / globals::time_step[nt].width / spectra[nt].delta_freq[nnu] / 4.e12 / PI / PARSEC / PARSEC / globals::nprocs * anglefactor;

    spectra[nt].flux[nnu] += deltaE;

    if (stokes_i != NULL)
    {
      stokes_i[nt].flux[nnu] += pkt_ptr->stokes[0] * deltaE;
      // assert_always(spectra[nt].nu_min == stokes_i[nt].nu_min);
      // assert_always(spectra[nt].nu_max == stokes_i[nt].nu_max);
    }
    if (stokes_q != NULL)
    {
      stokes_q[nt].flux[nnu] += pkt_ptr->stokes[1] * deltaE;
      // assert_always(spectra[nt].nu_min == stokes_q[nt].nu_min);
      // assert_always(spectra[nt].nu_max == stokes_q[nt].nu_max);
    }
    if (stokes_u != NULL)
    {
      stokes_u[nt].flux[nnu] += pkt_ptr->stokes[2] * deltaE;
      // assert_always(spectra[nt].nu_min == stokes_u[nt].nu_min);
      // assert_always(spectra[nt].nu_max == stokes_u[nt].nu_max);
    }

    if (spectra[nt].do_emission_res)
    {
      const int nproc = columnindex_from_emissiontype(pkt_ptr->emissiontype);
      spectra[nt].stat[nnu].emission[nproc] += deltaE;

      const int truenproc = columnindex_from_emissiontype(pkt_ptr->trueemissiontype);
      spectra[nt].stat[nnu].trueemission[truenproc] += deltaE;

      if (stokes_i != NULL && stokes_i[nt].do_emission_res)
        stokes_i[nt].stat[nnu].emission[nproc] += pkt_ptr->stokes[0] * deltaE;
      if (stokes_q != NULL && stokes_q[nt].do_emission_res)
        stokes_q[nt].stat[nnu].emission[nproc] += pkt_ptr->stokes[1] * deltaE;
      if (stokes_u != NULL && stokes_u[nt].do_emission_res)
        stokes_u[nt].stat[nnu].emission[nproc] += pkt_ptr->stokes[2] * deltaE;

      if (TRACE_EMISSION_ABSORPTION_REGION_ON && (current_abin == -1))
      {
        const int et = pkt_ptr->trueemissiontype;
        if (et >= 0)
        {
          if (t_arrive >= traceemissabs_timemin && t_arrive <= traceemissabs_timemax)
          {
            if (pkt_ptr->nu_rf >= traceemissabs_nulower && pkt_ptr->nu_rf <= traceemissabs_nuupper)
            {
              traceemissionabsorption[et].energyemitted += deltaE;

              traceemissionabsorption[et].emission_weightedvelocity_sum += pkt_ptr->trueemissionvelocity * deltaE;

              traceemission_totalenergy += deltaE;
            }
          }
        }
      }

      const int nnu_abs = (log(pkt_ptr->absorptionfreq) - log(nu_min)) /  dlognu;
      if (nnu_abs >= 0 && nnu_abs < MNUBINS)
      {
        const double deltaE_absorption = pkt_ptr->e_rf / globals::time_step[nt].width / spectra[nt].delta_freq[nnu_abs] / 4.e12 / PI / PARSEC / PARSEC / globals::nprocs * anglefactor;
        const int at = pkt_ptr->absorptiontype;
        if (at >= 0)
        {
          /// bb-emission
          const int element = globals::linelist[at].elementindex;
          const int ion = globals::linelist[at].ionindex;
          spectra[nt].stat[nnu_abs].absorption[element * get_max_nions() + ion] += deltaE_absorption;

          if (stokes_i != NULL && stokes_i[nt].do_emission_res)
            stokes_i[nt].stat[nnu].absorption[nproc] += pkt_ptr->stokes[0] * deltaE_absorption;
          if (stokes_q != NULL && stokes_q[nt].do_emission_res)
            stokes_q[nt].stat[nnu].absorption[nproc] += pkt_ptr->stokes[1] * deltaE_absorption;
          if (stokes_u != NULL && stokes_u[nt].do_emission_res)
            stokes_u[nt].stat[nnu].absorption[nproc] += pkt_ptr->stokes[2] * deltaE_absorption;

          if (TRACE_EMISSION_ABSORPTION_REGION_ON && t_arrive >= traceemissabs_timemin && t_arrive <= traceemissabs_timemax)
          {
            if ((current_abin == -1) && (pkt_ptr->nu_rf >= traceemissabs_nulower) && (pkt_ptr->nu_rf <= traceemissabs_nuupper))
            {
              traceemissionabsorption[at].energyabsorbed += deltaE_absorption;

              double vel_vec[3];
              get_velocity(pkt_ptr->em_pos, vel_vec, pkt_ptr->em_time);
              traceemissionabsorption[at].absorption_weightedvelocity_sum += vec_len(vel_vec) * deltaE_absorption;

              traceabsorption_totalenergy += deltaE_absorption;
            }
          }
        }
      }
    }
  }
}


void init_spectrum_trace(void)
{
  if (TRACE_EMISSION_ABSORPTION_REGION_ON)
  {
    traceemission_totalenergy = 0.;
    traceemissionabsorption = (struct emissionabsorptioncontrib *) malloc(globals::nlines * sizeof(emissionabsorptioncontrib));
    traceabsorption_totalenergy = 0.;
    for (int i = 0; i < globals::nlines; i++)
    {
      traceemissionabsorption[i].energyemitted = 0.;
      traceemissionabsorption[i].emission_weightedvelocity_sum = 0.;
      traceemissionabsorption[i].energyabsorbed = 0.;
      traceemissionabsorption[i].absorption_weightedvelocity_sum = 0.;
      traceemissionabsorption[i].lineindex = i; // this will be important when the list gets sorted
    }
  }
}


void free_spectra(struct spec *spectra)
{
  for (int n = 0; n < globals::ntstep; n++)
  {
    if (spectra[n].do_emission_res)
    {
      for (int m = 0; m < globals::nnubins; m++)
      {
        free(spectra[n].stat[m].absorption);
        free(spectra[n].stat[m].emission);
        free(spectra[n].stat[m].trueemission);
      }
    }
  }

  free(spectra);
}


void init_spectra(struct spec *spectra, const double nu_min, const double nu_max, const bool do_emission_res)
{
  // start by setting up the time and frequency bins.
  // it is all done interms of a logarithmic spacing in both t and nu - get the
  // step sizes first.
  ///Should be moved to input.c or exspec.c
  const double dlognu = (log(nu_max) - log(nu_min)) / globals::nnubins;

  for (int n = 0; n < globals::ntstep; n++)
  {
    spectra[n].nu_min = nu_min;
    spectra[n].nu_max = nu_max;
    for (int m = 0; m < globals::nnubins; m++)
    {
      spectra[n].lower_freq[m] = exp(log(nu_min) + (m * (dlognu)));
      spectra[n].delta_freq[m] = exp(log(nu_min) + ((m + 1) * (dlognu))) - spectra[n].lower_freq[m];
      spectra[n].flux[m] = 0.0;

      if (do_emission_res)
      {
        for (int i = 0; i < 2 * get_nelements() * get_max_nions() + 1; i++)
        {
          spectra[n].stat[m].emission[i] = 0;
          spectra[n].stat[m].trueemission[i] = 0;
        }

        for (int i = 0; i < get_nelements() * get_max_nions(); i++)
        {
          spectra[n].stat[m].absorption[i] = 0;
        }
      }
    }
  }
}


struct spec *alloc_spectra(const bool do_emission_res)
{
  struct spec *spectra = (struct spec *) calloc(MTBINS, sizeof(struct spec));
  /// Check if enough  memory for spectra has been assigned
  /// and allocate memory for the emission statistics
  if (globals::nnubins > MNUBINS)
  {
    printout("WARNING: Too many frequency bins in spectrum - reducing.\n");
    globals::nnubins = MNUBINS;
  }
  assert_always(globals::ntstep < MTBINS);

  for (int n = 0; n < globals::ntstep; n++)
  {
    spectra[n].do_emission_res = do_emission_res;
    if (spectra[n].do_emission_res)
    {
      for (int m = 0; m < globals::nnubins; m++)
      {
        const int emissioncount = 2 * get_nelements() * get_max_nions() + 1;
        spectra[n].stat[m].absorption = (double *) calloc(get_nelements() * get_max_nions(), sizeof(double));
        spectra[n].stat[m].emission = (double *) calloc(emissioncount, sizeof(double));
        spectra[n].stat[m].trueemission = (double *) calloc(emissioncount, sizeof(double));

        assert_always(spectra[n].stat[m].absorption != NULL);
        assert_always(spectra[n].stat[m].emission != NULL);
        assert_always(spectra[n].stat[m].trueemission != NULL);
      }
    }
  }

  return spectra;
}


void add_to_spec_res(
  const PKT *const pkt_ptr, int current_abin, struct spec *spectra,
  struct spec *stokes_i, struct spec *stokes_q, struct spec *stokes_u)
// Routine to add a packet to the outgoing spectrum.
{
  /* Need to (1) decide which time bin to put it in and (2) which frequency bin. */

  /* Time bin - we know that it escaped at "escape_time". However, we have to allow
     for travel time. Use the formula in Leon's paper.
     The extra distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
  */

  if (current_abin == -1)
  {
    // angle averaged spectrum
    add_to_spec(pkt_ptr, current_abin, spectra, stokes_i, stokes_q, stokes_u);
  }
  else
  {
    double xhat[3] = {1.0, 0.0, 0.0};

    /// Angle resolved case: need to work out the correct angle bin
    const double costheta = dot(pkt_ptr->dir, globals::syn_dir);
    const double thetabin = ((costheta + 1.0) * sqrt(MABINS) / 2.0);
    double vec1[3];
    cross_prod(pkt_ptr->dir, globals::syn_dir, vec1);
    double vec2[3];
    cross_prod(xhat, globals::syn_dir, vec2);
    const double cosphi = dot(vec1,vec2) / vec_len(vec1) / vec_len(vec2);

    double vec3[3];
    cross_prod(vec2, globals::syn_dir, vec3);
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
    const int na = (thetabin * sqrt(MABINS)) + phibin;

    /// Add only packets which escape to the current angle bin
    if (na == current_abin)
    {
      add_to_spec(pkt_ptr, current_abin, spectra, stokes_i, stokes_q, stokes_u);
    }
  }
}
