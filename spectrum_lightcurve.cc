#include "spectrum_lightcurve.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <ios>
#include <iterator>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

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

// the other "atomicadd" function is atomic only for multithreaded modes (STDPAR or OpenMP), but here we need it to be
// atomic for node-shared memory between processes even in single-threaded mode
template <typename T, typename U>
constexpr void atomicadd_always(T &var, U &&val) {
  std::atomic_ref<T>(var).fetch_add(std::forward<U>(val), std::memory_order_relaxed);
}

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

        const auto upper_uniquelevelindex = get_uniquelevelindex(element, ion, upper);

        const auto alltrans_startdown = get_alltrans_startdown(upper_uniquelevelindex);
        const auto ndowntrans = get_ndowntrans(upper_uniquelevelindex);
        int downtransid = -1;
        for (int alltransindex = alltrans_startdown; alltransindex < alltrans_startdown + ndowntrans; alltransindex++) {
          if (globals::alltrans.targetlevelindex[alltransindex] == lower) {
            downtransid = alltransindex;
            break;
          }
        }
        assert_always(downtransid != -1);

        printout("%7.2e (%5.1f%%) %4d %9d %5d %5d %8.1f %8.2e %4d %7.1f %7.1f %7.1e %7.1e\n", encontrib,
                 100 * encontrib / totalenergy, get_atomicnumber(element), get_ionstage(element, ion),
                 globals::linelist[lineindex].upperlevelindex, globals::linelist[lineindex].lowerlevelindex,
                 globals::alltrans.coll_str[downtransid], globals::linelist[lineindex].einstein_A,
                 static_cast<int>(globals::alltrans.forbidden[downtransid]), linelambda, v_rad, B_lu, B_ul);
      } else {
        break;
      }
    }
    printout("\n");
  }

  traceemissionabsorption.clear();
}

// number of different emission processes (bf and bb for each ion, and free-free)
auto get_proccount() -> int { return (2 * get_nelements() * get_max_nions()) + 1; }

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
  }
  // bf-emission
  const int bfindex = -1 - et;
  if (globals::nbfcontinua == 0) {
    // assert_always(false);  // if there are no bf processes, we should not get here
    return 2 * get_nelements() * get_max_nions();
  }
  assert_always(bfindex < globals::nbfcontinua);
  const int element = globals::bflist[bfindex].elementindex;
  const int ion = globals::bflist[bfindex].ionindex;
  const int level = globals::bflist[bfindex].levelindex;
  const int phixstargetindex = globals::bflist[bfindex].phixstargetindex;
  const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);

  assert_always(get_emtype_continuum(element, ion, level, upperionlevel) == et);

  return (get_nelements() * get_max_nions()) + (element * get_max_nions()) + ion;
}

[[nodiscard]] auto get_absindex(const ptrdiff_t nts, const ptrdiff_t nnu_abs) -> ptrdiff_t {
  const ptrdiff_t nelements = get_nelements();
  const ptrdiff_t max_nions = get_max_nions();
  return (nnu_abs * globals::ntimesteps * nelements * max_nions) + (nts * nelements * max_nions);
}

void write_specpol_param(std::ostream &specpol_file, std::ostream &emissionpol_file, std::ostream &absorptionpol_file,
                         const Spectra &spec, const int nnu, const bool do_emission_absorption) {
  const int proccount = get_proccount();
  const int ioncount = get_nelements() * get_max_nions();  // may be higher than the true included ion count
  // Stokes I, Q, or U
  const auto ntimesteps = static_cast<ptrdiff_t>(globals::ntimesteps);
  for (ptrdiff_t nts = 0; nts < ntimesteps; nts++) {
    specpol_file << spec.fluxalltimesteps[(nnu * ntimesteps) + nts] << ' ';

    if (do_emission_absorption) {
      for (int nproc = 0; nproc < proccount; nproc++) {
        const auto emindex = (nnu * ntimesteps * proccount) + (nts * proccount) + nproc;
        emissionpol_file << spec.emissionalltimesteps[emindex] << ' ';
      }
      emissionpol_file << '\n';

      for (int i = 0; i < ioncount; i++) {
        absorptionpol_file << spec.absorptionalltimesteps[get_absindex(nts, nnu) + i] << ' ';
      }
      absorptionpol_file << '\n';
    }
  }
}

void write_partial_lightcurve_spectra_dirbin(const int nts, std::span<const Packet> pkts,
                                             const bool do_emission_absorption, const int dirbin) {
  thread_local static std::vector<double> rpkt_light_curve_lum;
  thread_local static std::vector<double> rpkt_light_curve_lumcmf;
  thread_local static std::vector<double> gamma_light_curve_lum;
  thread_local static std::vector<double> gamma_light_curve_lumcmf;
  resize_exactly(rpkt_light_curve_lum, globals::ntimesteps);
  resize_exactly(rpkt_light_curve_lumcmf, globals::ntimesteps);
  resize_exactly(gamma_light_curve_lum, globals::ntimesteps);
  resize_exactly(gamma_light_curve_lumcmf, globals::ntimesteps);
  std::ranges::fill(rpkt_light_curve_lum, 0.);
  std::ranges::fill(rpkt_light_curve_lumcmf, 0.);
  std::ranges::fill(gamma_light_curve_lum, 0.);
  std::ranges::fill(gamma_light_curve_lumcmf, 0.);

  TRACE_EMISSION_ABSORPTION_REGION_ON = false;

  init_spectra(rpkt_spectra, NU_MIN_R, NU_MAX_R, do_emission_absorption);

  MPI_Barrier(globals::mpi_comm_node);
#if defined REPRODUCIBLE && REPRODUCIBLE
  for (int node_rank = 0; node_rank < globals::node_nprocs; node_rank++) {
    // do one rank at a time to keep the results reproducible (instead of simultaneous atomic adds to shared memory)
    if (node_rank == globals::rank_in_node) {
#else
  {
    {
#endif
      for (int ii = 0; ii < globals::npkts; ii++) {
        if (pkts[ii].type == TYPE_ESCAPE) {
          if (pkts[ii].escape_type == TYPE_RPKT) {
            add_to_lc_res(pkts[ii], dirbin, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
            add_to_spec_res(pkts[ii], dirbin, rpkt_spectra, nullptr, nullptr, nullptr);
          } else if (dirbin == -1 && pkts[ii].escape_type == TYPE_GAMMA) {
            add_to_lc_res(pkts[ii], dirbin, gamma_light_curve_lum, gamma_light_curve_lumcmf);
          }
        }
      }
    }
    MPI_Barrier(globals::mpi_comm_node);
  }

  const int numtimesteps = nts + 1;  // only produce spectra and light curves up to one past nts
  assert_always(numtimesteps <= globals::ntimesteps);

  MPI_Barrier(MPI_COMM_WORLD);
  if (globals::rank_in_node == 0) {
    MPI_Reduce_safe(rpkt_spectra.fluxalltimesteps, MPI_SUM, 0, globals::mpi_comm_internode);
    if (rpkt_spectra.do_emission_absorption) {
      MPI_Reduce_safe(rpkt_spectra.absorptionalltimesteps, MPI_SUM, 0, globals::mpi_comm_internode);
      MPI_Reduce_safe(rpkt_spectra.emissionalltimesteps, MPI_SUM, 0, globals::mpi_comm_internode);
      MPI_Reduce_safe(rpkt_spectra.trueemissionalltimesteps, MPI_SUM, 0, globals::mpi_comm_internode);
    }
  }
  MPI_Reduce_safe(rpkt_light_curve_lum, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce_safe(rpkt_light_curve_lumcmf, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce_safe(gamma_light_curve_lum, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce_safe(gamma_light_curve_lumcmf, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (globals::my_rank == 0) {
    if (dirbin == -1) {
      write_light_curve("light_curve.out", dirbin, rpkt_light_curve_lum, rpkt_light_curve_lumcmf, numtimesteps);
      write_light_curve("gamma_light_curve.out", dirbin, gamma_light_curve_lum, gamma_light_curve_lumcmf, numtimesteps);
      write_spectrum("spec.out", "emission.out", "emissiontrue.out", "absorption.out", rpkt_spectra, numtimesteps);
    } else {
      if (!std::filesystem::exists(outdir_resfiles)) {
        std::filesystem::create_directory(outdir_resfiles);
      }
      write_light_curve(std::format("{}light_curve_res_{:02d}.out", outdir_resfiles, dirbin), dirbin,
                        rpkt_light_curve_lum, rpkt_light_curve_lumcmf, numtimesteps);
      write_spectrum(std::format("{}spec_res_{:02d}.out", outdir_resfiles, dirbin),
                     std::format("{}emission_res_{:02d}.out", outdir_resfiles, dirbin),
                     std::format("{}emissiontrue_res_{:02d}.out", outdir_resfiles, dirbin),
                     std::format("{}absorption_res_{:02d}.out", outdir_resfiles, dirbin), rpkt_spectra, numtimesteps);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

}  // anonymous namespace

void write_spectrum(const std::string &spec_filename, const std::string &emission_filename,
                    const std::string &trueemission_filename, const std::string &absorption_filename,
                    const Spectra &spectra, const int numtimesteps) {
  auto spec_file = fstream_required(spec_filename, std::ios::out | std::ios::trunc);

  const bool do_emission_absorption = spectra.do_emission_absorption;
  std::fstream emission_file{};
  std::fstream trueemission_file{};
  std::fstream absorption_file{};

  if (do_emission_absorption) {
    emission_file = fstream_required(emission_filename, std::ios::out | std::ios::trunc);
    trueemission_file = fstream_required(trueemission_filename, std::ios::out | std::ios::trunc);
    absorption_file = fstream_required(absorption_filename, std::ios::out | std::ios::trunc);
  }

  if (TRACE_EMISSION_ABSORPTION_REGION_ON && do_emission_absorption && !traceemissionabsorption.empty()) {
    printout_tracemission_stats();
  }

  assert_always(numtimesteps <= globals::ntimesteps);

  spec_file << "0 ";
  for (int p = 0; p < numtimesteps; p++) {
    spec_file << globals::timesteps[p].mid / DAY << ' ';
  }
  spec_file << '\n';

  const auto ntimesteps_all = static_cast<ptrdiff_t>(globals::ntimesteps);

  const int proccount = get_proccount();
  const int ioncount = get_nelements() * get_max_nions();  // may be higher than the true included ion count
  for (ptrdiff_t nubin = 0; nubin < MNUBINS; nubin++) {
    spec_file << (spectra.lower_freq[nubin] + (spectra.delta_freq[nubin] / 2)) << ' ';

    for (ptrdiff_t nts = 0; nts < numtimesteps; nts++) {
      spec_file << spectra.fluxalltimesteps[(nubin * ntimesteps_all) + nts] << ' ';

      if (do_emission_absorption) {
        const auto emindex_nts_nubin = (nubin * ntimesteps_all * proccount) + (nts * proccount);
        for (int nproc = 0; nproc < proccount; nproc++) {
          emission_file << spectra.emissionalltimesteps[emindex_nts_nubin + nproc] << ' ';
        }
        emission_file << '\n';

        for (int truenproc = 0; truenproc < proccount; truenproc++) {
          trueemission_file << spectra.trueemissionalltimesteps[emindex_nts_nubin + truenproc] << ' ';
        }
        trueemission_file << '\n';

        for (int i = 0; i < ioncount; i++) {
          absorption_file << spectra.absorptionalltimesteps[get_absindex(nts, nubin) + i] << ' ';
        }
        absorption_file << '\n';
      }
    }
    spec_file << '\n';
  }
}

void write_specpol(const std::string &specpol_filename, const std::string &emission_filename,
                   const std::string &absorption_filename, const Spectra *stokes_i, const Spectra *stokes_q,
                   const Spectra *stokes_u) {
  auto specpol_file = fstream_required(specpol_filename, std::ios::out | std::ios::trunc);
  std::fstream emissionpol_file{};
  std::fstream absorptionpol_file{};

  const bool do_emission_absorption = stokes_i->do_emission_absorption;

  if (do_emission_absorption) {
    emissionpol_file = fstream_required(emission_filename, std::ios::out | std::ios::trunc);
    absorptionpol_file = fstream_required(absorption_filename, std::ios::out | std::ios::trunc);
    printout("Writing %s, %s, and %s\n", specpol_filename.c_str(), emission_filename.c_str(),
             absorption_filename.c_str());
  } else {
    printout("Writing %s\n", specpol_filename.c_str());
  }

  specpol_file << 0.0 << ' ';

  for (int l = 0; l < 3; l++) {
    for (int p = 0; p < globals::ntimesteps; p++) {
      specpol_file << globals::timesteps[p].mid / DAY << ' ';
    }
  }

  specpol_file << '\n';

  assert_always(std::ssize(stokes_i->delta_freq) == MNUBINS);
  assert_always(std::ssize(stokes_i->lower_freq) == MNUBINS);
  for (int nnu = 0; nnu < std::ssize(stokes_i->lower_freq); nnu++) {
    specpol_file << (stokes_i->lower_freq[nnu] + (stokes_i->delta_freq[nnu] / 2)) << ' ';

    write_specpol_param(specpol_file, emissionpol_file, absorptionpol_file, *stokes_i, nnu, do_emission_absorption);
    write_specpol_param(specpol_file, emissionpol_file, absorptionpol_file, *stokes_q, nnu, do_emission_absorption);
    write_specpol_param(specpol_file, emissionpol_file, absorptionpol_file, *stokes_u, nnu, do_emission_absorption);

    specpol_file << '\n';
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
void init_spectra(Spectra &spectra, const double nu_min, const double nu_max, const bool do_emission_absorption) {
  // setup the time and frequency bins using a logarithmic spacing in both t and nu

  assert_always(MNUBINS > 0);
  const double dlognu = (log(nu_max) - log(nu_min)) / MNUBINS;

  spectra.nu_min = nu_min;
  spectra.nu_max = nu_max;
  spectra.do_emission_absorption = do_emission_absorption;
  const bool print_memusage =
      (spectra.fluxalltimesteps.empty() || (do_emission_absorption && spectra.absorptionalltimesteps.empty()));

  for (ptrdiff_t nnu = 0; nnu < MNUBINS; nnu++) {
    spectra.lower_freq[nnu] = static_cast<float>(std::exp(log(nu_min) + (nnu * (dlognu))));
    spectra.delta_freq[nnu] =
        static_cast<float>(std::exp(log(nu_min) + ((nnu + 1) * (dlognu))) - spectra.lower_freq[nnu]);
  }

  if (spectra.fluxalltimesteps.empty()) {
    assert_always(spectra.win_fluxalltimesteps == MPI_WIN_NULL);
    std::tie(spectra.fluxalltimesteps, spectra.win_fluxalltimesteps) =
        MPI_shared_malloc_span_keepwin<double>(globals::ntimesteps * MNUBINS);
  }
  assert_always(std::ssize(spectra.fluxalltimesteps) == globals::ntimesteps * MNUBINS);
  std::ranges::fill(spectra.fluxalltimesteps, 0.0);
  MPI_Barrier(globals::mpi_comm_node);

  if (do_emission_absorption) {
    if (spectra.absorptionalltimesteps.empty()) {
      assert_always(spectra.win_absorptionalltimesteps == MPI_WIN_NULL);
      std::tie(spectra.absorptionalltimesteps, spectra.win_absorptionalltimesteps) =
          MPI_shared_malloc_span_keepwin<double>(globals::ntimesteps * MNUBINS * get_nelements() * get_max_nions());
    }
    assert_always(std::ssize(spectra.absorptionalltimesteps) ==
                  globals::ntimesteps * MNUBINS * get_nelements() * get_max_nions());
    std::ranges::fill(spectra.absorptionalltimesteps, 0.0);

    if (spectra.emissionalltimesteps.empty()) {
      assert_always(spectra.win_emissionalltimesteps == MPI_WIN_NULL);
      std::tie(spectra.emissionalltimesteps, spectra.win_emissionalltimesteps) =
          MPI_shared_malloc_span_keepwin<double>(globals::ntimesteps * MNUBINS * get_proccount());
    }
    assert_always(std::ssize(spectra.emissionalltimesteps) == globals::ntimesteps * MNUBINS * get_proccount());
    std::ranges::fill(spectra.emissionalltimesteps, 0.0);

    if (spectra.trueemissionalltimesteps.empty()) {
      assert_always(spectra.win_trueemissionalltimesteps == MPI_WIN_NULL);
      std::tie(spectra.trueemissionalltimesteps, spectra.win_trueemissionalltimesteps) =
          MPI_shared_malloc_span_keepwin<double>(globals::ntimesteps * MNUBINS * get_proccount());
    }
    assert_always(std::ssize(spectra.trueemissionalltimesteps) == globals::ntimesteps * MNUBINS * get_proccount());
    std::ranges::fill(spectra.trueemissionalltimesteps, 0.0);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (print_memusage) {
    logprintlnfmt("[info] mem_usage: set of spectra{} occupy {:.3f} MB (node shared memory)",
                  do_emission_absorption ? " (with emission/absorption tracing)" : "",
                  spectra.mem_usage_bytes() / 1024. / 1024.);
  }
}

// Add a packet to the outgoing spectrum.
void add_to_spec_res(const Packet &pkt, const int dirbin, Spectra &spectra, Spectra *stokes_i, Spectra *stokes_q,
                     Spectra *stokes_u) {
  if (dirbin != -1 && get_escapedirectionbin(pkt.dir) != dirbin) {
    return;  // do not add to the spectrum if the direction bin does not match
  }

  // Need to (1) decide which time bin to put it in and (2) which frequency bin.

  // specific angle bins contain fewer packets than the full sphere, so must be normalised to match
  const double nu_min = spectra.nu_min;
  const double nu_max = spectra.nu_max;
  const double t_arrive = get_arrive_time(pkt);
  if (t_arrive > globals::tmin && t_arrive < globals::tmax && pkt.nu_rf > nu_min && pkt.nu_rf < nu_max) {
    const double dlognu = (log(nu_max) - log(nu_min)) / MNUBINS;
    const auto nts = get_timestep(t_arrive);

    // a binary search into freq_lower would probably be faster than this double logarithm
    const auto nnu = static_cast<ptrdiff_t>((log(pkt.nu_rf) - log(nu_min)) / dlognu);

    const double anglefactor = (dirbin >= 0) ? MABINS : 1.;
    const double deltaE = pkt.e_rf / globals::timesteps[nts].width / spectra.delta_freq.at(nnu) / 4.e12 / PI / PARSEC /
                          PARSEC / globals::nprocs_exspec * anglefactor;

    const auto fluxindex = (nnu * static_cast<ptrdiff_t>(globals::ntimesteps)) + nts;
    atomicadd_always(spectra.fluxalltimesteps[fluxindex], deltaE);

    if (stokes_i != nullptr) {
      atomicadd_always(stokes_i->fluxalltimesteps[fluxindex], pkt.stokes[0] * deltaE);
    }
    if (stokes_q != nullptr) {
      atomicadd_always(stokes_q->fluxalltimesteps[fluxindex], pkt.stokes[1] * deltaE);
    }
    if (stokes_u != nullptr) {
      atomicadd_always(stokes_u->fluxalltimesteps[fluxindex], pkt.stokes[2] * deltaE);
    }

    if (spectra.do_emission_absorption) {
      const auto proccount = get_proccount();

      const auto truenproc = columnindex_from_emissiontype(pkt.trueemissiontype);
      assert_always(truenproc < proccount);
      if (truenproc >= 0) {
        const auto emindex = (nnu * globals::ntimesteps * proccount) + (nts * proccount) + truenproc;
        atomicadd_always(spectra.trueemissionalltimesteps[emindex], deltaE);
      }

      const auto nproc = columnindex_from_emissiontype(pkt.emissiontype);
      assert_always(nproc < proccount);
      if (nproc >= 0) {  // -1 means EMTYPE_NOTSET
        const auto emindex = (nnu * globals::ntimesteps * proccount) + (nts * proccount) + nproc;
        atomicadd_always(spectra.emissionalltimesteps[emindex], deltaE);

        if (stokes_i != nullptr && stokes_i->do_emission_absorption) {
          atomicadd_always(stokes_i->emissionalltimesteps[emindex], pkt.stokes[0] * deltaE);
        }
        if (stokes_q != nullptr && stokes_q->do_emission_absorption) {
          atomicadd_always(stokes_q->emissionalltimesteps[emindex], pkt.stokes[1] * deltaE);
        }
        if (stokes_u != nullptr && stokes_u->do_emission_absorption) {
          atomicadd_always(stokes_u->emissionalltimesteps[emindex], pkt.stokes[2] * deltaE);
        }
      }

      if (TRACE_EMISSION_ABSORPTION_REGION_ON && (dirbin == -1)) {
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
          const auto absindex = get_absindex(nts, nnu_abs) + (element * get_max_nions()) + ion;
          atomicadd_always(spectra.absorptionalltimesteps[absindex], deltaE_absorption);

          if (stokes_i != nullptr && stokes_i->do_emission_absorption) {
            atomicadd_always(stokes_i->absorptionalltimesteps[absindex], pkt.stokes[0] * deltaE_absorption);
          }
          if (stokes_q != nullptr && stokes_q->do_emission_absorption) {
            atomicadd_always(stokes_q->absorptionalltimesteps[absindex], pkt.stokes[1] * deltaE_absorption);
          }
          if (stokes_u != nullptr && stokes_u->do_emission_absorption) {
            atomicadd_always(stokes_u->absorptionalltimesteps[absindex], pkt.stokes[2] * deltaE_absorption);
          }

          if (TRACE_EMISSION_ABSORPTION_REGION_ON && t_arrive >= traceemissabs_timemin &&
              t_arrive <= traceemissabs_timemax) {
            if ((dirbin == -1) && (pkt.nu_rf >= traceemissabs_nulower) && (pkt.nu_rf <= traceemissabs_nuupper)) {
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

void write_partial_lightcurve_spectra(const int nts, std::span<const Packet> pkts) {
  // this is called by sn3d (not exspec) when each rank has its own set of packets
  // in memory
  const bool simulation_complete = (nts >= globals::timestep_finish - 1);
  if (!simulation_complete && nts % 5 != 0) {
    // do not write spectra every timestep, only every 5th timestep
    return;
  }

  // the emission resolved spectra are slow to generate, and require a lot of memory
  const bool do_emission_absorption = WRITE_EMISSIONABSORPTION_SPEC_AT_END && simulation_complete;

  const bool multdimensional = grid::get_model_type() != GridType::SPHERICAL1D;
  const int dirbinend = (multdimensional && simulation_complete) ? MABINS : 0;

  const auto time_func_start = std::time(nullptr);

  for (int dirbin = -1; dirbin < dirbinend; dirbin++) {
    write_partial_lightcurve_spectra_dirbin(nts, pkts, do_emission_absorption, dirbin);
  }

  printout("timestep %d: Saving light curves and %sspectra took %lds\n", nts,
           do_emission_absorption ? "emission/absorption " : "", std::time(nullptr) - time_func_start);
}

void write_light_curve(const std::string &lc_filename, const int dirbin, const std::vector<double> &light_curve_lum,
                       const std::vector<double> &light_curve_lumcmf, const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntimesteps);

  auto lc_file = fstream_required(lc_filename, std::ios::out | std::ios::trunc);

  // UVOIR bolometric light curve
  for (int nts = 0; nts < numtimesteps; nts++) {
    lc_file << globals::timesteps[nts].mid / DAY << ' ' << light_curve_lum[nts] / LSUN << ' '
            << light_curve_lumcmf[nts] / LSUN << '\n';
  }

  if (dirbin == -1) {
    // Now print out the gamma ray deposition rate in the same file.
    for (int m = 0; m < numtimesteps; m++) {
      lc_file << globals::timesteps[m].mid / DAY << ' '
              << globals::timesteps[m].gamma_dep / LSUN / globals::timesteps[m].width << ' '
              << globals::timesteps[m].cmf_lum / globals::timesteps[m].width / LSUN << '\n';
    }
  }
}

// add a packet to the outgoing light-curve.
void add_to_lc_res(const Packet &pkt, const int dirbin, std::span<double> light_curve_lum,
                   std::span<double> light_curve_lumcmf) {
  const double anglefactor = (dirbin >= 0) ? MABINS : 1.;
  if (dirbin == -1) {
    // -1 means all full 4π angle average (no angle filtering)

    // Put this into the time grid
    const double arrive_time = get_arrive_time(pkt);
    if (arrive_time > globals::tmin && arrive_time < globals::tmax) {
      const int nts = get_timestep(arrive_time);
      atomicadd_always(light_curve_lum[nts],
                       pkt.e_rf / globals::timesteps[nts].width * anglefactor / globals::nprocs_exspec);
    }

    const double inverse_gamma = std::sqrt(1. - (globals::vmax * globals::vmax / CLIGHTSQUARED));

    // Now do the cmf light curve.
    const double arrive_time_cmf = pkt.escape_time * inverse_gamma;

    if (arrive_time_cmf > globals::tmin && arrive_time_cmf < globals::tmax) {
      const int nts = get_timestep(arrive_time_cmf);
      atomicadd_always(light_curve_lumcmf[nts], pkt.e_cmf / globals::timesteps[nts].width * anglefactor /
                                                    globals::nprocs_exspec / inverse_gamma);
    }

  } else if (get_escapedirectionbin(pkt.dir) == dirbin) {
    // packets that escape in the select angle bin

    const double t_arrive = get_arrive_time(pkt);
    if (t_arrive > globals::tmin && t_arrive < globals::tmax) {
      const int nts = get_timestep(t_arrive);
      atomicadd_always(light_curve_lum[nts],
                       pkt.e_rf / globals::timesteps[nts].width * anglefactor / globals::nprocs_exspec);
    }
  }
}
