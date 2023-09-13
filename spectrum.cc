#include "spectrum.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <memory>

#include "atomic.h"
#include "exspec.h"
#include "light_curve.h"
#include "sn3d.h"
#include "vectors.h"

bool TRACE_EMISSION_ABSORPTION_REGION_ON = false;

constexpr double traceemissabs_lambdamin = 1000.;  // in Angstroms
constexpr double traceemissabs_lambdamax = 25000.;
constexpr double traceemissabs_nulower = (1.e8 * CLIGHT / traceemissabs_lambdamax);
constexpr double traceemissabs_nuupper = (1.e8 * CLIGHT / traceemissabs_lambdamin);
constexpr double traceemissabs_timemin = (320. * DAY);
constexpr double traceemissabs_timemax = (340. * DAY);

using emissionabsorptioncontrib = struct emissionabsorptioncontrib {
  double energyemitted;
  double emission_weightedvelocity_sum;
  double energyabsorbed;
  double absorption_weightedvelocity_sum;
  int lineindex;  // this will be important when the list gets sorted
};

static struct emissionabsorptioncontrib *traceemissionabsorption = nullptr;
double traceemission_totalenergy = 0.;
double traceabsorption_totalenergy = 0.;

std::unique_ptr<struct spec> rpkt_spectra = nullptr;

static auto compare_emission(const void *p1, const void *p2) -> int {
  const auto *elem1 = static_cast<const struct emissionabsorptioncontrib *>(p1);
  const auto *elem2 = static_cast<const struct emissionabsorptioncontrib *>(p2);

  if (elem1->energyemitted < elem2->energyemitted) {
    return 1;
  }
  if (elem1->energyemitted > elem2->energyemitted) {
    return -1;
  }
  return 0;
}

static auto compare_absorption(const void *p1, const void *p2) -> int {
  const auto *elem1 = static_cast<const struct emissionabsorptioncontrib *>(p1);
  const auto *elem2 = static_cast<const struct emissionabsorptioncontrib *>(p2);

  if (elem1->energyabsorbed < elem2->energyabsorbed) {
    return 1;
  }
  if (elem1->energyabsorbed > elem2->energyabsorbed) {
    return -1;
  }
  return 0;
}

static void printout_tracemission_stats() {
  const int maxlinesprinted = 500;

  // mode is 0 for emission and 1 for absorption
  for (int mode = 0; mode < 2; mode++) {
    if (mode == 0) {
      qsort(traceemissionabsorption, globals::nlines, sizeof(emissionabsorptioncontrib), compare_emission);
      printout("lambda [%5.1f, %5.1f] nu %g %g\n", traceemissabs_lambdamin, traceemissabs_lambdamax,
               traceemissabs_nulower, traceemissabs_nuupper);

      printout("Top line emission contributions in the range lambda [%5.1f, %5.1f] time [%5.1fd, %5.1fd] (%g erg)\n",
               traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_timemin / DAY,
               traceemissabs_timemax / DAY, traceemission_totalenergy);
    } else {
      qsort(traceemissionabsorption, globals::nlines, sizeof(emissionabsorptioncontrib), compare_absorption);
      printout("Top line absorption contributions in the range lambda [%5.1f, %5.1f] time [%5.1fd, %5.1fd] (%g erg)\n",
               traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_timemin / DAY,
               traceemissabs_timemax / DAY, traceabsorption_totalenergy);
    }

    // display the top entries of the sorted list
    int nlines_limited = globals::nlines;
    if (globals::nlines > maxlinesprinted) {
      nlines_limited = maxlinesprinted;
    }
    printout("%17s %4s %9s %5s %5s %8s %8s %4s %7s %7s %7s %7s\n", "energy", "Z", "ion_stage", "upper", "lower",
             "coll_str", "A", "forb", "lambda", "<v_rad>", "B_lu", "B_ul");
    for (int i = 0; i < nlines_limited; i++) {
      double encontrib = NAN;
      double totalenergy = NAN;
      if (mode == 0) {
        encontrib = traceemissionabsorption[i].energyemitted;
        totalenergy = traceemission_totalenergy;
      } else {
        encontrib = traceemissionabsorption[i].energyabsorbed;
        totalenergy = traceabsorption_totalenergy;
      }
      if (encontrib > 0.)  // lines that emit/absorb some energy
      {
        const int lineindex = traceemissionabsorption[i].lineindex;
        const int element = globals::linelist[lineindex].elementindex;
        const int ion = globals::linelist[lineindex].ionindex;
        const double linelambda = 1e8 * CLIGHT / globals::linelist[lineindex].nu;
        // flux-weighted average radial velocity of emission in km/s
        double v_rad = NAN;
        if (mode == 0) {
          v_rad =
              traceemissionabsorption[i].emission_weightedvelocity_sum / traceemissionabsorption[i].energyemitted / 1e5;
        } else {
          v_rad = traceemissionabsorption[i].absorption_weightedvelocity_sum /
                  traceemissionabsorption[i].energyabsorbed / 1e5;
        }

        const int lower = globals::linelist[lineindex].lowerlevelindex;
        const int upper = globals::linelist[lineindex].upperlevelindex;

        const double statweight_target = stat_weight(element, ion, upper);
        const double statweight_lower = stat_weight(element, ion, lower);

        const double nu_trans = (epsilon(element, ion, upper) - epsilon(element, ion, lower)) / H;
        const double A_ul = einstein_spontaneous_emission(lineindex);
        const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
        const double B_lu = statweight_target / statweight_lower * B_ul;

        const int nupperdowntrans = get_ndowntrans(element, ion, upper);
        auto *downtranslist = globals::elements[element].ions[ion].levels[upper].downtrans;
        auto *downtrans = std::find_if(downtranslist, downtranslist + nupperdowntrans,
                                       [=](const auto &downtrans) { return downtrans.targetlevelindex == lower; });
        assert_always(downtrans != (downtranslist + nupperdowntrans));

        printout("%7.2e (%5.1f%%) %4d %9d %5d %5d %8.1f %8.2e %4d %7.1f %7.1f %7.1e %7.1e\n", encontrib,
                 100 * encontrib / totalenergy, get_atomicnumber(element), get_ionstage(element, ion),
                 globals::linelist[lineindex].upperlevelindex, globals::linelist[lineindex].lowerlevelindex,
                 downtrans->coll_str, einstein_spontaneous_emission(lineindex), downtrans->forbidden, linelambda, v_rad,
                 B_lu, B_ul);
      } else {
        break;
      }
    }
    printout("\n");
  }

  free(traceemissionabsorption);
  traceemissionabsorption = nullptr;
}

static auto get_proccount() -> int
// number of different emission processes (bf and bb for each ion, and free-free)
{
  return 2 * get_nelements() * get_max_nions() + 1;
}

void write_spectrum(const std::string &spec_filename, const char *emission_filename, const char *trueemission_filename,
                    const char *absorption_filename, const struct spec &spectra, int numtimesteps) {
  FILE *spec_file = fopen_required(spec_filename.c_str(), "w");

  FILE *emission_file = nullptr;
  FILE *trueemission_file = nullptr;
  FILE *absorption_file = nullptr;

  bool const do_emission_res = spectra.do_emission_res;

  if (do_emission_res) {
    emission_file = fopen_required(emission_filename, "w");
    assert_always(emission_file != nullptr);
    trueemission_file = fopen_required(trueemission_filename, "w");
    assert_always(trueemission_file != nullptr);
    absorption_file = fopen_required(absorption_filename, "w");
    assert_always(absorption_file != nullptr);
    printout("Writing %s, %s, %s, and %s\n", spec_filename.c_str(), emission_filename, trueemission_filename,
             absorption_filename);
  } else {
    printout("Writing %s\n", spec_filename.c_str());
  }

  if (TRACE_EMISSION_ABSORPTION_REGION_ON && do_emission_res && traceemissionabsorption != nullptr) {
    printout_tracemission_stats();
  }

  assert_always(numtimesteps <= globals::ntstep);

  fprintf(spec_file, "%g ", 0.0);
  for (int p = 0; p < numtimesteps; p++) {
    fprintf(spec_file, "%g ", globals::time_step[p].mid / DAY);
  }
  fprintf(spec_file, "\n");

  const int proccount = get_proccount();
  const int ioncount = get_nelements() * get_max_nions();  // may be higher than the true included ion count
  for (int nnu = 0; nnu < MNUBINS; nnu++) {
    fprintf(spec_file, "%g ", ((spectra.lower_freq[nnu] + (spectra.delta_freq[nnu] / 2))));

    for (int nts = 0; nts < numtimesteps; nts++) {
      fprintf(spec_file, "%g ", spectra.timesteps[nts].flux[nnu]);
      if (do_emission_res) {
        for (int i = 0; i < proccount; i++) {
          fprintf(emission_file, "%g ", spectra.timesteps[nts].emission[nnu * proccount + i]);
        }
        fprintf(emission_file, "\n");

        for (int i = 0; i < proccount; i++) {
          fprintf(trueemission_file, "%g ", spectra.timesteps[nts].trueemission[nnu * proccount + i]);
        }
        fprintf(trueemission_file, "\n");

        for (int i = 0; i < ioncount; i++) {
          fprintf(absorption_file, "%g ", spectra.timesteps[nts].absorption[nnu * ioncount + i]);
        }
        fprintf(absorption_file, "\n");
      }
    }
    fprintf(spec_file, "\n");
  }

  fclose(spec_file);
  if (do_emission_res) {
    fclose(emission_file);
    fclose(trueemission_file);
    fclose(absorption_file);
  }
}

void write_specpol(const std::string &specpol_filename, const std::string &emission_filename,
                   const std::string &absorption_filename, const struct spec *stokes_i, const struct spec *stokes_q,
                   const struct spec *stokes_u) {
  FILE *specpol_file = fopen_required(specpol_filename.c_str(), "w");
  FILE *emissionpol_file = nullptr;
  FILE *absorptionpol_file = nullptr;

  const bool do_emission_res = stokes_i->do_emission_res;

  if (do_emission_res) {
    emissionpol_file = fopen_required(emission_filename.c_str(), "w");
    absorptionpol_file = fopen_required(absorption_filename.c_str(), "w");
    printout("Writing %s, %s, and %s\n", specpol_filename.c_str(), emission_filename.c_str(),
             absorption_filename.c_str());
  } else {
    printout("Writing %s\n", specpol_filename.c_str());
  }

  fprintf(specpol_file, "%g ", 0.0);

  for (int l = 0; l < 3; l++) {
    for (int p = 0; p < globals::ntstep; p++) {
      fprintf(specpol_file, "%g ", globals::time_step[p].mid / DAY);
    }
  }

  fprintf(specpol_file, "\n");

  const int proccount = get_proccount();
  const int ioncount = get_nelements() * get_max_nions();
  for (int m = 0; m < MNUBINS; m++) {
    fprintf(specpol_file, "%g ", ((stokes_i->lower_freq[m] + (stokes_i->delta_freq[m] / 2))));

    // Stokes I
    for (int p = 0; p < globals::ntstep; p++) {
      fprintf(specpol_file, "%g ", stokes_i->timesteps[p].flux[m]);

      if (do_emission_res) {
        for (int i = 0; i < proccount; i++) {
          fprintf(emissionpol_file, "%g ", stokes_i->timesteps[p].emission[m * proccount + i]);
        }
        fprintf(emissionpol_file, "\n");

        for (int i = 0; i < ioncount; i++) {
          fprintf(absorptionpol_file, "%g ", stokes_i->timesteps[p].absorption[m * ioncount + i]);
        }
        fprintf(absorptionpol_file, "\n");
      }
    }

    // Stokes Q
    for (int p = 0; p < globals::ntstep; p++) {
      fprintf(specpol_file, "%g ", stokes_q->timesteps[p].flux[m]);

      if (do_emission_res) {
        for (int i = 0; i < proccount; i++) {
          fprintf(emissionpol_file, "%g ", stokes_q->timesteps[p].emission[m * proccount + i]);
        }
        fprintf(emissionpol_file, "\n");

        for (int i = 0; i < ioncount; i++) {
          fprintf(absorptionpol_file, "%g ", stokes_q->timesteps[p].absorption[m * ioncount + i]);
        }
        fprintf(absorptionpol_file, "\n");
      }
    }

    // Stokes U
    for (int p = 0; p < globals::ntstep; p++) {
      fprintf(specpol_file, "%g ", stokes_u->timesteps[p].flux[m]);

      if (do_emission_res) {
        for (int i = 0; i < proccount; i++) {
          fprintf(emissionpol_file, "%g ", stokes_u->timesteps[p].emission[m * proccount + i]);
        }
        fprintf(emissionpol_file, "\n");

        for (int i = 0; i < ioncount; i++) {
          fprintf(absorptionpol_file, "%g ", stokes_u->timesteps[p].absorption[m * ioncount + i]);
        }
        fprintf(absorptionpol_file, "\n");
      }
    }

    fprintf(specpol_file, "\n");
  }

  fclose(specpol_file);
  if (do_emission_res) {
    fclose(emissionpol_file);
    fclose(absorptionpol_file);
  }
}

static auto columnindex_from_emissiontype(const int et) -> int {
  if (et >= 0) {
    /// bb-emission
    const int element = globals::linelist[et].elementindex;
    const int ion = globals::linelist[et].ionindex;
    return element * get_max_nions() + ion;
  }
  if (et == EMTYPE_FREEFREE) {
    /// ff-emission

    const int contindex = -1 - et;
    assert_always(contindex >= globals::nbfcontinua);  // make sure the special value didn't collide with a real process

    return 2 * get_nelements() * get_max_nions();
  }
  if (et == EMTYPE_NOTSET) {
    return -1;
  }  /// bf-emission
  const int contindex = -1 - et;
  if (globals::nbfcontinua == 0) {
    // assert_always(false);  // if there are no bf processes, we should not get here
    return 2 * get_nelements() * get_max_nions();
  }
  assert_always(contindex < globals::nbfcontinua);
  const int element = globals::bflist[contindex].elementindex;
  const int ion = globals::bflist[contindex].ionindex;
  const int level = globals::bflist[contindex].levelindex;
  const int phixstargetindex = globals::bflist[contindex].phixstargetindex;
  const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);

  assert_always(get_continuumindex(element, ion, level, upperionlevel) == et);

  return get_nelements() * get_max_nions() + element * get_max_nions() + ion;
}

static void add_to_spec(const struct packet *const pkt_ptr, const int current_abin, struct spec &spectra,
                        const struct spec *stokes_i, const struct spec *stokes_q, const struct spec *stokes_u)
// Routine to add a packet to the outgoing spectrum.
{
  // Need to (1) decide which time bin to put it in and (2) which frequency bin.

  // specific angle bins contain fewer packets than the full sphere, so must be normalised to match
  const double anglefactor = (current_abin >= 0) ? MABINS : 1.;

  const double nu_min = spectra.nu_min;
  const double nu_max = spectra.nu_max;
  const double t_arrive = get_arrive_time(pkt_ptr);
  if (t_arrive > globals::tmin && t_arrive < globals::tmax && pkt_ptr->nu_rf > nu_min && pkt_ptr->nu_rf < nu_max) {
    const int nt = get_timestep(t_arrive);
    const double dlognu = (log(nu_max) - log(nu_min)) / MNUBINS;

    const int nnu = static_cast<int>((log(pkt_ptr->nu_rf) - log(nu_min)) / dlognu);
    assert_always(nnu < MNUBINS);

    const double deltaE = pkt_ptr->e_rf / globals::time_step[nt].width / spectra.delta_freq[nnu] / 4.e12 / PI / PARSEC /
                          PARSEC / globals::nprocs * anglefactor;

    spectra.timesteps[nt].flux[nnu] += deltaE;

    if (stokes_i != nullptr) {
      stokes_i->timesteps[nt].flux[nnu] += pkt_ptr->stokes[0] * deltaE;
    }
    if (stokes_q != nullptr) {
      stokes_q->timesteps[nt].flux[nnu] += pkt_ptr->stokes[1] * deltaE;
    }
    if (stokes_u != nullptr) {
      stokes_u->timesteps[nt].flux[nnu] += pkt_ptr->stokes[2] * deltaE;
    }

    if (spectra.do_emission_res) {
      const int proccount = get_proccount();

      const int truenproc = columnindex_from_emissiontype(pkt_ptr->trueemissiontype);
      assert_always(truenproc < proccount);
      if (truenproc >= 0) {
        spectra.timesteps[nt].trueemission[nnu * proccount + truenproc] += deltaE;
      }

      const int nproc = columnindex_from_emissiontype(pkt_ptr->emissiontype);
      assert_always(nproc < proccount);
      if (nproc >= 0) {  // -1 means not set
        spectra.timesteps[nt].emission[nnu * proccount + nproc] += deltaE;

        if (stokes_i != nullptr && stokes_i->do_emission_res) {
          stokes_i->timesteps[nt].emission[nnu * proccount + nproc] += pkt_ptr->stokes[0] * deltaE;
        }
        if (stokes_q != nullptr && stokes_q->do_emission_res) {
          stokes_q->timesteps[nt].emission[nnu * proccount + nproc] += pkt_ptr->stokes[1] * deltaE;
        }
        if (stokes_u != nullptr && stokes_u->do_emission_res) {
          stokes_u->timesteps[nt].emission[nnu * proccount + nproc] += pkt_ptr->stokes[2] * deltaE;
        }
      }

      if (TRACE_EMISSION_ABSORPTION_REGION_ON && (current_abin == -1)) {
        const int et = pkt_ptr->trueemissiontype;
        if (et >= 0) {
          if (t_arrive >= traceemissabs_timemin && t_arrive <= traceemissabs_timemax) {
            if (pkt_ptr->nu_rf >= traceemissabs_nulower && pkt_ptr->nu_rf <= traceemissabs_nuupper) {
              traceemissionabsorption[et].energyemitted += deltaE;

              traceemissionabsorption[et].emission_weightedvelocity_sum += pkt_ptr->trueemissionvelocity * deltaE;

              traceemission_totalenergy += deltaE;
            }
          }
        }
      }

      const int nnu_abs = static_cast<int>((log(pkt_ptr->absorptionfreq) - log(nu_min)) / dlognu);
      if (nnu_abs >= 0 && nnu_abs < MNUBINS) {
        const int ioncount = get_nelements() * get_max_nions();
        const double deltaE_absorption = pkt_ptr->e_rf / globals::time_step[nt].width / spectra.delta_freq[nnu_abs] /
                                         4.e12 / PI / PARSEC / PARSEC / globals::nprocs * anglefactor;
        const int at = pkt_ptr->absorptiontype;
        if (at >= 0) {
          /// bb-emission
          const int element = globals::linelist[at].elementindex;
          const int ion = globals::linelist[at].ionindex;
          spectra.timesteps[nt].absorption[nnu_abs * ioncount + element * get_max_nions() + ion] += deltaE_absorption;

          if (stokes_i != nullptr && stokes_i->do_emission_res) {
            stokes_i->timesteps[nt].absorption[nnu_abs * ioncount + element * get_max_nions() + ion] +=
                pkt_ptr->stokes[0] * deltaE_absorption;
          }
          if (stokes_q != nullptr && stokes_q->do_emission_res) {
            stokes_q->timesteps[nt].absorption[nnu_abs * ioncount + element * get_max_nions() + ion] +=
                pkt_ptr->stokes[1] * deltaE_absorption;
          }
          if (stokes_u != nullptr && stokes_u->do_emission_res) {
            stokes_u->timesteps[nt].absorption[nnu_abs * ioncount + element * get_max_nions() + ion] +=
                pkt_ptr->stokes[2] * deltaE_absorption;
          }

          if (TRACE_EMISSION_ABSORPTION_REGION_ON && t_arrive >= traceemissabs_timemin &&
              t_arrive <= traceemissabs_timemax) {
            if ((current_abin == -1) && (pkt_ptr->nu_rf >= traceemissabs_nulower) &&
                (pkt_ptr->nu_rf <= traceemissabs_nuupper)) {
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

void init_spectrum_trace() {
  if (TRACE_EMISSION_ABSORPTION_REGION_ON) {
    traceemission_totalenergy = 0.;
    traceemissionabsorption =
        static_cast<struct emissionabsorptioncontrib *>(malloc(globals::nlines * sizeof(emissionabsorptioncontrib)));
    traceabsorption_totalenergy = 0.;
    for (int i = 0; i < globals::nlines; i++) {
      traceemissionabsorption[i].energyemitted = 0.;
      traceemissionabsorption[i].emission_weightedvelocity_sum = 0.;
      traceemissionabsorption[i].energyabsorbed = 0.;
      traceemissionabsorption[i].absorption_weightedvelocity_sum = 0.;
      traceemissionabsorption[i].lineindex = i;  // this will be important when the list gets sorted
    }
  }
}

void init_spectra(struct spec &spectra, const double nu_min, const double nu_max, const bool do_emission_res) {
  // start by setting up the time and frequency bins.
  // it is all done interms of a logarithmic spacing in both t and nu - get the
  // step sizes first.

  assert_always(MNUBINS > 0);

  const double dlognu = (log(nu_max) - log(nu_min)) / MNUBINS;

  spectra.nu_min = nu_min;
  spectra.nu_max = nu_max;
  spectra.do_emission_res = do_emission_res;
  assert_always(spectra.lower_freq.get() != nullptr);
  assert_always(spectra.delta_freq.get() != nullptr);
  for (int nnu = 0; nnu < MNUBINS; nnu++) {
    spectra.lower_freq[nnu] = exp(log(nu_min) + (nnu * (dlognu)));
    spectra.delta_freq[nnu] = exp(log(nu_min) + ((nnu + 1) * (dlognu))) - spectra.lower_freq[nnu];
  }

  assert_always(spectra.timesteps != nullptr);

  const int proccount = get_proccount();
  const int ioncount = get_nelements() * get_max_nions();
  for (int nts = 0; nts < globals::ntstep; nts++) {
    assert_always(spectra.timesteps[nts].flux != nullptr);
    for (int nnu = 0; nnu < MNUBINS; nnu++) {
      spectra.timesteps[nts].flux[nnu] = 0.0;
    }

    if (do_emission_res) {
      assert_always(spectra.timesteps[nts].emission != nullptr);
      assert_always(spectra.timesteps[nts].trueemission != nullptr);
      for (int i = 0; i < MNUBINS * proccount; i++) {
        spectra.timesteps[nts].emission[i] = 0;
        spectra.timesteps[nts].trueemission[i] = 0;
      }

      assert_always(spectra.timesteps[nts].absorption != nullptr);
      for (int i = 0; i < MNUBINS * ioncount; i++) {
        spectra.timesteps[nts].absorption[i] = 0;
      }
    }
  }
}

static void alloc_emissionabsorption_spectra(auto &spectra) {
  size_t mem_usage = 0;
  const int proccount = get_proccount();
  spectra->do_emission_res = true;

  mem_usage += globals::ntstep * MNUBINS * get_nelements() * get_max_nions() * sizeof(double);
  spectra->absorptionalltimesteps =
      std::make_unique<double[]>(globals::ntstep * MNUBINS * get_nelements() * get_max_nions());
  assert_always(spectra->absorptionalltimesteps != nullptr);

  mem_usage += 2 * globals::ntstep * MNUBINS * proccount * sizeof(double);
  spectra->emissionalltimesteps = std::make_unique<double[]>(globals::ntstep * MNUBINS * proccount);
  assert_always(spectra->emissionalltimesteps != nullptr);

  spectra->trueemissionalltimesteps = std::make_unique<double[]>(globals::ntstep * MNUBINS * proccount);
  assert_always(spectra->trueemissionalltimesteps != nullptr);

  for (int nts = 0; nts < globals::ntstep; nts++) {
    assert_always(spectra->timesteps[nts].absorption == nullptr);
    assert_always(spectra->timesteps[nts].emission == nullptr);
    assert_always(spectra->timesteps[nts].trueemission == nullptr);

    spectra->timesteps[nts].absorption =
        &spectra->absorptionalltimesteps[nts * MNUBINS * get_nelements() * get_max_nions()];

    spectra->timesteps[nts].emission = &spectra->emissionalltimesteps[nts * MNUBINS * proccount];

    spectra->timesteps[nts].trueemission = &spectra->trueemissionalltimesteps[nts * MNUBINS * proccount];

    assert_always(spectra->timesteps[nts].absorption != nullptr);
    assert_always(spectra->timesteps[nts].emission != nullptr);
    assert_always(spectra->timesteps[nts].trueemission != nullptr);
  }
  printout("[info] mem_usage: allocated set of emission/absorption spectra occupying total of %.3f MB (nnubins %d)\n",
           mem_usage / 1024. / 1024., MNUBINS);
}

auto alloc_spectra(const bool do_emission_res) -> std::unique_ptr<struct spec> {
  size_t mem_usage = 0;
  assert_always(globals::ntstep > 0);
  // std::unique_ptr<struct spec> spectra(new struct spec);
  auto spectra = std::make_unique<struct spec>();
  mem_usage += globals::ntstep * sizeof(struct spec);

  spectra->do_emission_res = false;  // might be set true later by alloc_emissionabsorption_spectra
  spectra->lower_freq = std::make_unique<float[]>(MNUBINS);
  spectra->delta_freq = std::make_unique<float[]>(MNUBINS);

  spectra->timesteps = std::make_unique<struct timestepspec[]>(globals::ntstep);
  mem_usage += globals::ntstep * sizeof(struct timestepspec);

  spectra->fluxalltimesteps = std::make_unique<double[]>(globals::ntstep * MNUBINS);
  mem_usage += globals::ntstep * MNUBINS * sizeof(double);

  assert_always(MNUBINS > 0);
  for (int nts = 0; nts < globals::ntstep; nts++) {
    spectra->timesteps[nts].flux = &spectra->fluxalltimesteps[nts * MNUBINS];

    spectra->timesteps[nts].absorption = nullptr;
    spectra->timesteps[nts].emission = nullptr;
    spectra->timesteps[nts].trueemission = nullptr;
  }

  printout("[info] mem_usage: allocated set of spectra occupying total of %.3f MB (nnubins %d)\n",
           mem_usage / 1024. / 1024., MNUBINS);

  if (do_emission_res) {
    alloc_emissionabsorption_spectra(spectra);
  }

  return spectra;
}

void add_to_spec_res(const struct packet *const pkt_ptr, int current_abin, struct spec &spectra,
                     const struct spec *stokes_i, const struct spec *stokes_q, const struct spec *stokes_u)
// Routine to add a packet to the outgoing spectrum.
{
  // Need to (1) decide which time bin to put it in and (2) which frequency bin.

  // Time bin - we know that it escaped at "escape_time". However, we have to allow
  // for travel time. Use the formula in Leon's paper.
  // The extra distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).

  if (current_abin == -1 || get_escapedirectionbin(pkt_ptr->dir, globals::syn_dir) == current_abin) {
    // either angle average spectrum or packet matches the selected angle bin
    add_to_spec(pkt_ptr, current_abin, spectra, stokes_i, stokes_q, stokes_u);
  }
}

#ifdef MPI_ON
static void mpi_reduce_spectra(int my_rank, struct spec *spectra, int numtimesteps) {
  for (int n = 0; n < numtimesteps; n++) {
    MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : spectra->timesteps[n].flux, spectra->timesteps[n].flux, MNUBINS,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (spectra->do_emission_res) {
      const int proccount = get_proccount();
      MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : spectra->timesteps[n].absorption, spectra->timesteps[n].absorption,
                 MNUBINS * get_nelements() * get_max_nions(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : spectra->timesteps[n].emission, spectra->timesteps[n].emission,
                 MNUBINS * proccount, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : spectra->timesteps[n].trueemission, spectra->timesteps[n].trueemission,
                 MNUBINS * proccount, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  }
}
#endif

void write_partial_lightcurve_spectra(int my_rank, int nts, struct packet *pkts) {
  const time_t time_func_start = time(nullptr);

  std::vector<double> rpkt_light_curve_lum(globals::ntstep, 0.);
  std::vector<double> rpkt_light_curve_lumcmf(globals::ntstep, 0.);
  std::vector<double> gamma_light_curve_lum(globals::ntstep, 0.);
  std::vector<double> gamma_light_curve_lumcmf(globals::ntstep, 0.);

  TRACE_EMISSION_ABSORPTION_REGION_ON = false;

  bool do_emission_res = WRITE_PARTIAL_EMISSIONABSORPTIONSPEC ? globals::do_emission_res : false;

  if (rpkt_spectra == nullptr) {
    rpkt_spectra = alloc_spectra(do_emission_res);
    assert_always(rpkt_spectra != nullptr);
  }

  struct spec *stokes_i = nullptr;
  struct spec *stokes_q = nullptr;
  struct spec *stokes_u = nullptr;

  // the emission resolved spectra are slow to generate, so only allow making them for the final timestep or every n
  if (WRITE_PARTIAL_EMISSIONABSORPTIONSPEC && globals::do_emission_res) {
    if ((nts >= globals::ftstep - 1) || (nts % 5 == 0)) {
      do_emission_res = true;
    }
  }

  init_spectra(*rpkt_spectra, NU_MIN_R, NU_MAX_R, do_emission_res);

  for (int ii = 0; ii < globals::npkts; ii++) {
    if (pkts[ii].type == TYPE_ESCAPE) {
      const int abin = -1;  // all angles
      if (pkts[ii].escape_type == TYPE_RPKT) {
        add_to_lc_res(&pkts[ii], abin, rpkt_light_curve_lum.data(), rpkt_light_curve_lumcmf.data());
        add_to_spec_res(&pkts[ii], abin, *rpkt_spectra, stokes_i, stokes_q, stokes_u);
      } else if (abin == -1 && pkts[ii].escape_type == TYPE_GAMMA) {
        add_to_lc_res(&pkts[ii], abin, gamma_light_curve_lum.data(), gamma_light_curve_lumcmf.data());
      }
    }
  }

  const int numtimesteps = nts + 1;  // only produce spectra and light curves up to one past nts
  assert_always(numtimesteps <= globals::ntstep);

  const time_t time_mpireduction_start = time(nullptr);
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
  mpi_reduce_spectra(my_rank, rpkt_spectra.get(), numtimesteps);
  MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : rpkt_light_curve_lum.data(), rpkt_light_curve_lum.data(), numtimesteps,
             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : rpkt_light_curve_lumcmf.data(), rpkt_light_curve_lumcmf.data(), numtimesteps,
             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : gamma_light_curve_lum.data(), gamma_light_curve_lum.data(), numtimesteps,
             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : gamma_light_curve_lumcmf.data(), gamma_light_curve_lumcmf.data(),
             numtimesteps, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  const time_t time_mpireduction_end = time(nullptr);

  if (my_rank == 0) {
    write_light_curve("light_curve.out", -1, rpkt_light_curve_lum.data(), rpkt_light_curve_lumcmf.data(), numtimesteps);
    write_light_curve("gamma_light_curve.out", -1, gamma_light_curve_lum.data(), gamma_light_curve_lumcmf.data(),
                      numtimesteps);
    write_spectrum("spec.out", "emission.out", "emissiontrue.out", "absorption.out", *rpkt_spectra, numtimesteps);
  }

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  printout("timestep %d: Saving partial light curves and %sspectra took %lds (%lds for MPI reduction)\n", nts,
           do_emission_res ? "emission/absorption " : "", time(nullptr) - time_func_start,
           time_mpireduction_end - time_mpireduction_start);
}