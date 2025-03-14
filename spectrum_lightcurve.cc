#include "spectrum_lightcurve.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <functional>
#include <ios>
#include <string>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "exspec.h"
#include "globals.h"
#include "packet.h"
#include "sn3d.h"
#include "vectors.h"

namespace {

bool TRACE_EMISSION_ABSORPTION_REGION_ON = false;

constexpr double traceemissabs_lambdamin = 1000.;  // in Angstroms
constexpr double traceemissabs_lambdamax = 25000.;
constexpr double traceemissabs_nulower = (1.e8 * CLIGHT / traceemissabs_lambdamax);
constexpr double traceemissabs_nuupper = (1.e8 * CLIGHT / traceemissabs_lambdamin);
constexpr double traceemissabs_timemin = (320. * DAY);
constexpr double traceemissabs_timemax = (340. * DAY);

struct emissionabsorptioncontrib {
  double energyemitted;
  double emission_weightedvelocity_sum;
  double energyabsorbed;
  double absorption_weightedvelocity_sum;
  int lineindex;  // this will be important when the list gets sorted
};

std::vector<emissionabsorptioncontrib> traceemissionabsorption;
double traceemission_totalenergy = 0.;
double traceabsorption_totalenergy = 0.;

Spectra rpkt_spectra;

void printout_tracemission_stats() {
  const int maxlinesprinted = 500;

  // mode is 0 for emission and 1 for absorption
  for (int mode = 0; mode < 2; mode++) {
    if (mode == 0) {
      std::ranges::SORT_OR_STABLE_SORT(traceemissionabsorption,
                                       [](const auto &a, const auto &b) { return a.energyemitted > b.energyemitted; });
      printout("lambda [%5.1f, %5.1f] nu %g %g\n", traceemissabs_lambdamin, traceemissabs_lambdamax,
               traceemissabs_nulower, traceemissabs_nuupper);

      printout("Top line emission contributions in the range lambda [%5.1f, %5.1f] time [%5.1fd, %5.1fd] (%g erg)\n",
               traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_timemin / DAY,
               traceemissabs_timemax / DAY, traceemission_totalenergy);
    } else {
      std::ranges::SORT_OR_STABLE_SORT(traceemissionabsorption, std::ranges::greater{},
                                       &emissionabsorptioncontrib::energyabsorbed);
      printout("Top line absorption contributions in the range lambda [%5.1f, %5.1f] time [%5.1fd, %5.1fd] (%g erg)\n",
               traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_timemin / DAY,
               traceemissabs_timemax / DAY, traceabsorption_totalenergy);
    }

    // display the top entries of the sorted list
    int nlines_limited = globals::nlines;
    if (globals::nlines > maxlinesprinted) {
      nlines_limited = maxlinesprinted;
    }
    printout("%17s %4s %9s %5s %5s %8s %8s %4s %7s %7s %7s %7s\n", "energy", "Z", "ionstage", "upper", "lower",
             "coll_str", "A", "forb", "lambda", "<v_rad>", "B_lu", "B_ul");
    for (int i = 0; i < nlines_limited; i++) {
      double encontrib{NAN};
      double totalenergy{NAN};
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
        double v_rad{NAN};
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
        const double A_ul = globals::linelist[lineindex].einstein_A;
        const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
        const double B_lu = statweight_target / statweight_lower * B_ul;

        const int nupperdowntrans = get_ndowntrans(element, ion, upper);
        const auto *downtranslist = get_downtranslist(element, ion, upper);
        const auto *downtransition = std::find_if(downtranslist, downtranslist + nupperdowntrans,
                                                  [=](const auto &downtr) { return downtr.targetlevelindex == lower; });
        assert_always(downtransition != (downtranslist + nupperdowntrans));

        printout("%7.2e (%5.1f%%) %4d %9d %5d %5d %8.1f %8.2e %4d %7.1f %7.1f %7.1e %7.1e\n", encontrib,
                 100 * encontrib / totalenergy, get_atomicnumber(element), get_ionstage(element, ion),
                 globals::linelist[lineindex].upperlevelindex, globals::linelist[lineindex].lowerlevelindex,
                 downtransition->coll_str, globals::linelist[lineindex].einstein_A,
                 static_cast<int>(downtransition->forbidden), linelambda, v_rad, B_lu, B_ul);
      } else {
        break;
      }
    }
    printout("\n");
  }

  traceemissionabsorption.clear();
}

auto get_proccount() -> int
// number of different emission processes (bf and bb for each ion, and free-free)
{
  return (2 * get_nelements() * get_max_nions()) + 1;
}

auto columnindex_from_emissiontype(const int et) -> int {
  if (et >= 0) {
    // bb-emission
    const int element = globals::linelist[et].elementindex;
    const int ion = globals::linelist[et].ionindex;
    return (element * get_max_nions()) + ion;
  }
  if (et == EMTYPE_FREEFREE) {
    // ff-emission

    const int contindex = -1 - et;
    assert_always(contindex >= globals::nbfcontinua);  // make sure the special value didn't collide with a real process

    return 2 * get_nelements() * get_max_nions();
  }
  if (et == EMTYPE_NOTSET) {
    return -1;
  }  // bf-emission
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

  assert_always(get_emtype_continuum(element, ion, level, upperionlevel) == et);

  return (get_nelements() * get_max_nions()) + (element * get_max_nions()) + ion;
}

[[nodiscard]] auto get_absindex(const int nts, const int nnu_abs, const int element, const int ion) -> ptrdiff_t {
  const auto nelements = get_nelements();
  const auto max_nions = get_max_nions();
  return (static_cast<ptrdiff_t>(nts) * MNUBINS * nelements * max_nions) + (nnu_abs * nelements * max_nions) +
         (element * max_nions) + ion;
}

void add_to_spec(const Packet &pkt, const int current_abin, Spectra &spectra, Spectra *stokes_i, Spectra *stokes_q,
                 Spectra *stokes_u)
// Routine to add a packet to the outgoing spectrum.
{
  // Need to (1) decide which time bin to put it in and (2) which frequency bin.

  // specific angle bins contain fewer packets than the full sphere, so must be normalised to match
  const double anglefactor = (current_abin >= 0) ? MABINS : 1.;

  const double nu_min = spectra.nu_min;
  const double nu_max = spectra.nu_max;
  const double dlognu = (log(nu_max) - log(nu_min)) / MNUBINS;
  const double t_arrive = get_arrive_time(pkt);
  if (t_arrive > globals::tmin && t_arrive < globals::tmax && pkt.nu_rf > nu_min && pkt.nu_rf < nu_max) {
    const int nts = get_timestep(t_arrive);

    const int nnu = static_cast<int>((log(pkt.nu_rf) - log(nu_min)) / dlognu);
    assert_always(nnu < MNUBINS);

    const double deltaE = pkt.e_rf / globals::timesteps[nts].width / spectra.delta_freq[nnu] / 4.e12 / PI / PARSEC /
                          PARSEC / globals::nprocs_exspec * anglefactor;

    const ptrdiff_t fluxindex = (nts * MNUBINS) + nnu;
    spectra.fluxalltimesteps[fluxindex] += deltaE;

    if (stokes_i != nullptr) {
      stokes_i->fluxalltimesteps[fluxindex] += pkt.stokes[0] * deltaE;
    }
    if (stokes_q != nullptr) {
      stokes_q->fluxalltimesteps[fluxindex] += pkt.stokes[1] * deltaE;
    }
    if (stokes_u != nullptr) {
      stokes_u->fluxalltimesteps[fluxindex] += pkt.stokes[2] * deltaE;
    }

    if (spectra.do_emission_res) {
      const int proccount = get_proccount();

      const int truenproc = columnindex_from_emissiontype(pkt.trueemissiontype);
      assert_always(truenproc < proccount);
      if (truenproc >= 0) {
        const auto emindex = (static_cast<ptrdiff_t>(nts) * MNUBINS * proccount) + (nnu * proccount) + truenproc;
        spectra.trueemissionalltimesteps[emindex] += deltaE;
      }

      const int nproc = columnindex_from_emissiontype(pkt.emissiontype);
      assert_always(nproc < proccount);
      if (nproc >= 0) {  // -1 means not set
        const auto emindex = (static_cast<ptrdiff_t>(nts) * MNUBINS * proccount) + (nnu * proccount) + nproc;
        spectra.emissionalltimesteps[emindex] += deltaE;

        if (stokes_i != nullptr && stokes_i->do_emission_res) {
          stokes_i->emissionalltimesteps[emindex] += pkt.stokes[0] * deltaE;
        }
        if (stokes_q != nullptr && stokes_q->do_emission_res) {
          stokes_q->emissionalltimesteps[emindex] += pkt.stokes[1] * deltaE;
        }
        if (stokes_u != nullptr && stokes_u->do_emission_res) {
          stokes_u->emissionalltimesteps[emindex] += pkt.stokes[2] * deltaE;
        }
      }

      if (TRACE_EMISSION_ABSORPTION_REGION_ON && (current_abin == -1)) {
        const int et = pkt.trueemissiontype;
        if (et >= 0) {
          if (t_arrive >= traceemissabs_timemin && t_arrive <= traceemissabs_timemax) {
            if (pkt.nu_rf >= traceemissabs_nulower && pkt.nu_rf <= traceemissabs_nuupper) {
              traceemissionabsorption[et].energyemitted += deltaE;

              traceemissionabsorption[et].emission_weightedvelocity_sum += pkt.trueemissionvelocity * deltaE;

              traceemission_totalenergy += deltaE;
            }
          }
        }
      }

      const int nnu_abs = (pkt.absorptionfreq > 0 && std::isfinite(pkt.absorptionfreq))
                              ? static_cast<int>((log(pkt.absorptionfreq) - log(nu_min)) / dlognu)
                              : -1;
      if (nnu_abs >= 0 && nnu_abs < MNUBINS) {
        const double deltaE_absorption = pkt.e_rf / globals::timesteps[nts].width / spectra.delta_freq[nnu_abs] /
                                         4.e12 / PI / PARSEC / PARSEC / globals::nprocs_exspec * anglefactor;
        const int at = pkt.absorptiontype;
        if (at >= 0) {
          // bb-emission
          const int element = globals::linelist[at].elementindex;
          const int ion = globals::linelist[at].ionindex;
          const auto absindex = get_absindex(nts, nnu_abs, element, ion);
          spectra.absorptionalltimesteps[absindex] += deltaE_absorption;

          if (stokes_i != nullptr && stokes_i->do_emission_res) {
            stokes_i->absorptionalltimesteps[absindex] += pkt.stokes[0] * deltaE_absorption;
          }
          if (stokes_q != nullptr && stokes_q->do_emission_res) {
            stokes_q->absorptionalltimesteps[absindex] += pkt.stokes[1] * deltaE_absorption;
          }
          if (stokes_u != nullptr && stokes_u->do_emission_res) {
            stokes_u->absorptionalltimesteps[absindex] += pkt.stokes[2] * deltaE_absorption;
          }

          if (TRACE_EMISSION_ABSORPTION_REGION_ON && t_arrive >= traceemissabs_timemin &&
              t_arrive <= traceemissabs_timemax) {
            if ((current_abin == -1) && (pkt.nu_rf >= traceemissabs_nulower) && (pkt.nu_rf <= traceemissabs_nuupper)) {
              traceemissionabsorption[at].energyabsorbed += deltaE_absorption;

              const auto vel_vec = get_velocity(pkt.em_pos, pkt.em_time);
              traceemissionabsorption[at].absorption_weightedvelocity_sum += vec_len(vel_vec) * deltaE_absorption;

              traceabsorption_totalenergy += deltaE_absorption;
            }
          }
        }
      }
    }
  }
}

void mpi_reduce_spectra(int my_rank, Spectra &spectra) {
  MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : spectra.fluxalltimesteps.data(), spectra.fluxalltimesteps.data(),
             spectra.fluxalltimesteps.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (spectra.do_emission_res) {
    MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : spectra.absorptionalltimesteps.data(),
               spectra.absorptionalltimesteps.data(), spectra.absorptionalltimesteps.size(), MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);

    MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : spectra.emissionalltimesteps.data(), spectra.emissionalltimesteps.data(),
               spectra.emissionalltimesteps.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : spectra.trueemissionalltimesteps.data(),
               spectra.trueemissionalltimesteps.data(), spectra.trueemissionalltimesteps.size(), MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
  }
}

void write_specpol_param(FILE *specpol_file, FILE *emissionpol_file, FILE *absorptionpol_file, const Spectra &spec,
                         const int nnu, const bool do_emission_res) {
  const int proccount = get_proccount();
  const int ioncount = get_nelements() * get_max_nions();  // may be higher than the true included ion count
  // Stokes I, Q, or U
  for (int nts = 0; nts < globals::ntimesteps; nts++) {
    fprintf(specpol_file, "%g ", spec.fluxalltimesteps[(nts * MNUBINS) + nnu]);

    if (do_emission_res) {
      for (int nproc = 0; nproc < proccount; nproc++) {
        const auto emindex = (static_cast<ptrdiff_t>(nts) * MNUBINS * proccount) + (nnu * proccount) + nproc;
        fprintf(emissionpol_file, "%g ", spec.emissionalltimesteps[emindex]);
      }
      fprintf(emissionpol_file, "\n");

      for (int i = 0; i < ioncount; i++) {
        const auto absindex = get_absindex(nts, nnu, 0, i);
        fprintf(absorptionpol_file, "%g ", spec.absorptionalltimesteps[absindex]);
      }
      fprintf(absorptionpol_file, "\n");
    }
  }
}

}  // anonymous namespace

void write_spectrum(const std::string &spec_filename, const std::string &emission_filename,
                    const std::string &trueemission_filename, const std::string &absorption_filename,
                    const Spectra &spectra, const int numtimesteps) {
  FILE *spec_file = fopen_required(spec_filename, "w");

  FILE *emission_file{};
  FILE *trueemission_file{};
  FILE *absorption_file{};

  const bool do_emission_res = spectra.do_emission_res;

  if (do_emission_res) {
    emission_file = fopen_required(emission_filename, "w");
    assert_always(emission_file != nullptr);
    trueemission_file = fopen_required(trueemission_filename, "w");
    assert_always(trueemission_file != nullptr);
    absorption_file = fopen_required(absorption_filename, "w");
    assert_always(absorption_file != nullptr);
    printout("Writing %s, %s, %s, and %s\n", spec_filename.c_str(), emission_filename.c_str(),
             trueemission_filename.c_str(), absorption_filename.c_str());
  } else {
    printout("Writing %s\n", spec_filename.c_str());
  }

  if (TRACE_EMISSION_ABSORPTION_REGION_ON && do_emission_res && !traceemissionabsorption.empty()) {
    printout_tracemission_stats();
  }

  assert_always(numtimesteps <= globals::ntimesteps);

  fprintf(spec_file, "%g ", 0.0);
  for (int p = 0; p < numtimesteps; p++) {
    fprintf(spec_file, "%g ", globals::timesteps[p].mid / DAY);
  }
  fprintf(spec_file, "\n");

  const int proccount = get_proccount();
  const int ioncount = get_nelements() * get_max_nions();  // may be higher than the true included ion count
  for (int nnu = 0; nnu < MNUBINS; nnu++) {
    fprintf(spec_file, "%g ", ((spectra.lower_freq[nnu] + (spectra.delta_freq[nnu] / 2))));

    for (int nts = 0; nts < numtimesteps; nts++) {
      fprintf(spec_file, "%g ", spectra.fluxalltimesteps[(nts * MNUBINS) + nnu]);
      if (do_emission_res) {
        for (int nproc = 0; nproc < proccount; nproc++) {
          const auto emindex = (static_cast<ptrdiff_t>(nts) * MNUBINS * proccount) + (nnu * proccount) + nproc;
          fprintf(emission_file, "%g ", spectra.emissionalltimesteps[emindex]);
        }
        fprintf(emission_file, "\n");

        for (int truenproc = 0; truenproc < proccount; truenproc++) {
          const auto trueemindex = (static_cast<ptrdiff_t>(nts) * MNUBINS * proccount) + (nnu * proccount) + truenproc;
          fprintf(trueemission_file, "%g ", spectra.trueemissionalltimesteps[trueemindex]);
        }
        fprintf(trueemission_file, "\n");

        for (int i = 0; i < ioncount; i++) {
          const auto absindex = get_absindex(nts, nnu, 0, i);
          fprintf(absorption_file, "%g ", spectra.absorptionalltimesteps[absindex]);
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
                   const std::string &absorption_filename, const Spectra *stokes_i, const Spectra *stokes_q,
                   const Spectra *stokes_u) {
  FILE *specpol_file = fopen_required(specpol_filename, "w");
  FILE *emissionpol_file{};
  FILE *absorptionpol_file{};

  const bool do_emission_res = stokes_i->do_emission_res;

  if (do_emission_res) {
    emissionpol_file = fopen_required(emission_filename, "w");
    absorptionpol_file = fopen_required(absorption_filename, "w");
    printout("Writing %s, %s, and %s\n", specpol_filename.c_str(), emission_filename.c_str(),
             absorption_filename.c_str());
  } else {
    printout("Writing %s\n", specpol_filename.c_str());
  }

  fprintf(specpol_file, "%g ", 0.0);

  for (int l = 0; l < 3; l++) {
    for (int p = 0; p < globals::ntimesteps; p++) {
      fprintf(specpol_file, "%g ", globals::timesteps[p].mid / DAY);
    }
  }

  fprintf(specpol_file, "\n");

  assert_always(stokes_i->lower_freq.size() == stokes_i->delta_freq.size());
  for (ptrdiff_t nnu = 0; nnu < std::ssize(stokes_i->lower_freq); nnu++) {
    fprintf(specpol_file, "%g ", ((stokes_i->lower_freq[nnu] + (stokes_i->delta_freq[nnu] / 2))));

    write_specpol_param(specpol_file, emissionpol_file, absorptionpol_file, *stokes_i, nnu, do_emission_res);
    write_specpol_param(specpol_file, emissionpol_file, absorptionpol_file, *stokes_q, nnu, do_emission_res);
    write_specpol_param(specpol_file, emissionpol_file, absorptionpol_file, *stokes_u, nnu, do_emission_res);

    fprintf(specpol_file, "\n");
  }

  fclose(specpol_file);
  if (do_emission_res) {
    fclose(emissionpol_file);
    fclose(absorptionpol_file);
  }
}

void init_spectrum_trace() {
  if (TRACE_EMISSION_ABSORPTION_REGION_ON) {
    traceemission_totalenergy = 0.;
    resize_exactly(traceemissionabsorption, globals::nlines);
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

// resize and initialize the spectra object
void init_spectra(Spectra &spectra, const double nu_min, const double nu_max, const bool do_emission_res) {
  // setup the time and frequency bins using a logarithmic spacing in both t and nu

  assert_always(MNUBINS > 0);
  size_t mem_usage = 0;
  const double dlognu = (log(nu_max) - log(nu_min)) / MNUBINS;

  spectra.nu_min = nu_min;
  spectra.nu_max = nu_max;
  spectra.do_emission_res = do_emission_res;
  const bool print_memusage =
      (spectra.fluxalltimesteps.empty() || (do_emission_res && spectra.absorptionalltimesteps.empty()));

  for (ptrdiff_t nnu = 0; nnu < MNUBINS; nnu++) {
    spectra.lower_freq[nnu] = exp(log(nu_min) + (nnu * (dlognu)));
    spectra.delta_freq[nnu] = exp(log(nu_min) + ((nnu + 1) * (dlognu))) - spectra.lower_freq[nnu];
  }

  spectra.do_emission_res = do_emission_res;  // might be set true later by alloc_emissionabsorption_spectra

  resize_exactly(spectra.fluxalltimesteps, globals::ntimesteps * MNUBINS);
  std::ranges::fill(spectra.fluxalltimesteps, 0.0);

  mem_usage += globals::ntimesteps * sizeof(Spectra);
  mem_usage += globals::ntimesteps * MNUBINS * sizeof(double);

  if (do_emission_res) {
    const int proccount = get_proccount();

    mem_usage += globals::ntimesteps * MNUBINS * get_nelements() * get_max_nions() * sizeof(double);
    mem_usage += 2 * globals::ntimesteps * MNUBINS * proccount * sizeof(double);

    resize_exactly(spectra.absorptionalltimesteps, globals::ntimesteps * MNUBINS * get_nelements() * get_max_nions());
    resize_exactly(spectra.emissionalltimesteps, globals::ntimesteps * MNUBINS * proccount);
    resize_exactly(spectra.trueemissionalltimesteps, globals::ntimesteps * MNUBINS * proccount);

    std::ranges::fill(spectra.absorptionalltimesteps, 0.0);
    std::ranges::fill(spectra.emissionalltimesteps, 0.0);
    std::ranges::fill(spectra.trueemissionalltimesteps, 0.0);

    if (print_memusage) {
      printout("[info] mem_usage: set of emission/absorption spectra occupy %.3f MB (nnubins %d)\n",
               mem_usage / 1024. / 1024., MNUBINS);
    }

  } else {
    spectra.absorptionalltimesteps.clear();
    spectra.emissionalltimesteps.clear();
    spectra.trueemissionalltimesteps.clear();

    if (print_memusage) {
      printout("[info] mem_usage: set of spectra occupy %.3f MB (nnubins %d)\n", mem_usage / 1024. / 1024., MNUBINS);
    }
  }
}

// Add a packet to the outgoing spectrum.
void add_to_spec_res(const Packet &pkt, const int current_abin, Spectra &spectra, Spectra *stokes_i, Spectra *stokes_q,
                     Spectra *stokes_u) {
  if (current_abin == -1 || get_escapedirectionbin(pkt.dir) == current_abin) {
    // either angle average spectrum or packet matches the selected angle bin
    add_to_spec(pkt, current_abin, spectra, stokes_i, stokes_q, stokes_u);
  }
}

void write_partial_lightcurve_spectra(const int my_rank, const int nts, const Packet *pkts) {
  const auto time_func_start = std::time(nullptr);

  std::vector<double> rpkt_light_curve_lum(globals::ntimesteps, 0.);
  std::vector<double> rpkt_light_curve_lumcmf(globals::ntimesteps, 0.);
  std::vector<double> gamma_light_curve_lum(globals::ntimesteps, 0.);
  std::vector<double> gamma_light_curve_lumcmf(globals::ntimesteps, 0.);

  TRACE_EMISSION_ABSORPTION_REGION_ON = false;

  bool do_emission_res = false;

  // the emission resolved spectra are slow to generate, so only allow making them for the final timestep or every n
  if (WRITE_PARTIAL_EMISSIONABSORPTIONSPEC && globals::do_emission_res) {
    do_emission_res = ((nts >= globals::timestep_finish - 1) || (nts % 5 == 0));
  }

  init_spectra(rpkt_spectra, NU_MIN_R, NU_MAX_R, do_emission_res);

  for (int ii = 0; ii < globals::npkts; ii++) {
    if (pkts[ii].type == TYPE_ESCAPE) {
      const int abin = -1;  // all angles
      if (pkts[ii].escape_type == TYPE_RPKT) {
        add_to_lc_res(pkts[ii], abin, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
        add_to_spec_res(pkts[ii], abin, rpkt_spectra, nullptr, nullptr, nullptr);
      } else if (abin == -1 && pkts[ii].escape_type == TYPE_GAMMA) {
        add_to_lc_res(pkts[ii], abin, gamma_light_curve_lum, gamma_light_curve_lumcmf);
      }
    }
  }

  const int numtimesteps = nts + 1;  // only produce spectra and light curves up to one past nts
  assert_always(numtimesteps <= globals::ntimesteps);

  const auto time_mpireduction_start = std::time(nullptr);
  MPI_Barrier(MPI_COMM_WORLD);
  mpi_reduce_spectra(my_rank, rpkt_spectra);
  MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : rpkt_light_curve_lum.data(), rpkt_light_curve_lum.data(), numtimesteps,
             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : rpkt_light_curve_lumcmf.data(), rpkt_light_curve_lumcmf.data(), numtimesteps,
             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : gamma_light_curve_lum.data(), gamma_light_curve_lum.data(), numtimesteps,
             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : gamma_light_curve_lumcmf.data(), gamma_light_curve_lumcmf.data(),
             numtimesteps, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  const auto time_mpireduction_end = std::time(nullptr);

  if (my_rank == 0) {
    write_light_curve("light_curve.out", -1, rpkt_light_curve_lum, rpkt_light_curve_lumcmf, numtimesteps);
    write_light_curve("gamma_light_curve.out", -1, gamma_light_curve_lum, gamma_light_curve_lumcmf, numtimesteps);
    write_spectrum("spec.out", "emission.out", "emissiontrue.out", "absorption.out", rpkt_spectra, numtimesteps);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  printout("timestep %d: Saving partial light curves and %sspectra took %lds (%lds for MPI reduction)\n", nts,
           do_emission_res ? "emission/absorption " : "", std::time(nullptr) - time_func_start,
           time_mpireduction_end - time_mpireduction_start);
}

void write_light_curve(const std::string &lc_filename, const int current_abin,
                       const std::vector<double> &light_curve_lum, const std::vector<double> &light_curve_lumcmf,
                       const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntimesteps);

  auto lc_file = fstream_required(lc_filename, std::ios::out | std::ios::trunc);

  printout("Writing %s\n", lc_filename.c_str());

  char linebuffer[1024];

  // Print out the UVOIR bolometric light curve.
  for (int nts = 0; nts < numtimesteps; nts++) {
    assert_always(snprintf(linebuffer, sizeof(linebuffer), "%g %g %g", globals::timesteps[nts].mid / DAY,
                           (light_curve_lum[nts] / LSUN),
                           (light_curve_lumcmf[nts] / LSUN)) < static_cast<int>(sizeof(linebuffer)));
    lc_file << linebuffer << '\n';
  }

  if (current_abin == -1) {
    // Now print out the gamma ray deposition rate in the same file.
    for (int m = 0; m < numtimesteps; m++) {
      assert_always(snprintf(linebuffer, sizeof(linebuffer), "%g %g %g", globals::timesteps[m].mid / DAY,
                             (globals::timesteps[m].gamma_dep / LSUN / globals::timesteps[m].width),
                             (globals::timesteps[m].cmf_lum / globals::timesteps[m].width / LSUN)) <
                    static_cast<int>(sizeof(linebuffer)));
      lc_file << linebuffer << '\n';
    }
  }
}

// add a packet to the outgoing light-curve.
void add_to_lc_res(const Packet &pkt, const int current_abin, std::vector<double> &light_curve_lum,
                   std::vector<double> &light_curve_lumcmf) {
  if (current_abin == -1) {
    // Put this into the time grid
    const double arrive_time = get_arrive_time(pkt);
    if (arrive_time > globals::tmin && arrive_time < globals::tmax) {
      const int nts = get_timestep(arrive_time);
      atomicadd(light_curve_lum[nts], pkt.e_rf / globals::timesteps[nts].width / globals::nprocs_exspec);
    }

    const double inverse_gamma = std::sqrt(1. - (globals::vmax * globals::vmax / CLIGHTSQUARED));

    // Now do the cmf light curve.
    // t_arrive = pkt.escape_time * sqrt(1. - (vmax*vmax/CLIGHTSQUARED));
    const double arrive_time_cmf = pkt.escape_time * inverse_gamma;

    if (arrive_time_cmf > globals::tmin && arrive_time_cmf < globals::tmax) {
      const int nts = get_timestep(arrive_time_cmf);
      atomicadd(light_curve_lumcmf[nts],
                pkt.e_cmf / globals::timesteps[nts].width / globals::nprocs_exspec / inverse_gamma);
    }

    return;
  }
  if (get_escapedirectionbin(pkt.dir) == current_abin) {
    // Add only packets which escape to the current angle bin
    const double t_arrive = get_arrive_time(pkt);
    if (t_arrive > globals::tmin && t_arrive < globals::tmax) {
      const int nts = get_timestep(t_arrive);
      atomicadd(light_curve_lum[nts], pkt.e_rf / globals::timesteps[nts].width * MABINS / globals::nprocs_exspec);
    }
  }
}
