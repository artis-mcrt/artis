// Binning of escaped packets into spectra and light curves for each observer direction bin,
// optionally decomposed into the emission and absorption contributions of each atomic process.

#include "spectrum_lightcurve.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <ios>
#include <iterator>
#include <print>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "exspec.h"
#include "globals.h"
#include "grid.h"
#include "mpi_logging.h"
#include "packet.h"
#include "sn3d.h"
#include "vectors.h"

namespace {

// Set to true to print the strongest line emission/absorption contributions in a fixed wavelength and time
// window (see the traceemissabs_* limits below). This only works in exspec, which is the only caller of
// init_spectrum_trace(), i.e. the only place where the traceemissionabsorption list is allocated.
constexpr bool TRACE_EMISSION_ABSORPTION_REGION_ON = false;

constexpr double traceemissabs_lambdamin = 1000.;  // in Angstroms
constexpr double traceemissabs_lambdamax = 25000.;
constexpr double traceemissabs_nulower = (1.e8 * CLIGHT / traceemissabs_lambdamax);
constexpr double traceemissabs_nuupper = (1.e8 * CLIGHT / traceemissabs_lambdamin);
constexpr double traceemissabs_timemin = (320. * DAY);
constexpr double traceemissabs_timemax = (340. * DAY);

struct EmissionAbsorptionContrib {
  double energyemitted;
  double emission_weightedvelocity_sum;
  double energyabsorbed;
  double absorption_weightedvelocity_sum;
  int lineindex;  // this will be important when the list gets sorted
};

std::vector<EmissionAbsorptionContrib> traceemissionabsorption;
double traceemission_totalenergy = 0.;
double traceabsorption_totalenergy = 0.;

Spectra rpkt_spectra_I;

// the other "atomicadd" function is atomic only for multithreaded modes (STDPAR or OpenMP), but here we need it to be
// atomic for node-shared memory between processes even in single-threaded mode
template <typename T, typename U>
constexpr void atomicadd_always(T& var, U&& val) {
  std::atomic_ref<T>(var).fetch_add(std::forward<U>(val), std::memory_order_relaxed);
}

void printout_tracemission_stats() {
  const auto maxlinesprinted = 500Z;

  // mode is 0 for emission and 1 for absorption
  for (int mode = 0; mode < 2; mode++) {
    if (mode == 0) {
      std::ranges::SORT_OR_STABLE_SORT(traceemissionabsorption,
                                       [](const auto& a, const auto& b) { return a.energyemitted > b.energyemitted; });
      printlnlog("lambda [{:5.1f}, {:5.1f}] nu {:g} {:g}", traceemissabs_lambdamin, traceemissabs_lambdamax,
                 traceemissabs_nulower, traceemissabs_nuupper);

      printlnlog(
          "Top line emission contributions in the range lambda [{:5.1f}, {:5.1f}] time [{:5.1f}d, {:5.1f}d] ({:g} erg)",
          traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_timemin / DAY, traceemissabs_timemax / DAY,
          traceemission_totalenergy);
    } else {
      std::ranges::SORT_OR_STABLE_SORT(traceemissionabsorption, std::ranges::greater{},
                                       &EmissionAbsorptionContrib::energyabsorbed);
      printlnlog(
          "Top line absorption contributions in the range lambda [{:5.1f}, {:5.1f}] time [{:5.1f}d, {:5.1f}d] ({:g} "
          "erg)",
          traceemissabs_lambdamin, traceemissabs_lambdamax, traceemissabs_timemin / DAY, traceemissabs_timemax / DAY,
          traceabsorption_totalenergy);
    }

    printlnlog("{:>16} {:>4} {:>9} {:>5} {:>5} {:>8} {:>5} {:>7} {:>7} {:>7} {:>7}", "energy (frac)", "Z", "ionstage",
               "upper", "lower", "coll_str", "forb", "lambda", "<v_rad>", "B_lu", "B_ul");

    // display the top entries of the sorted list
    const auto nlines_limited = std::min(std::ssize(globals::linelist.nu), maxlinesprinted);
    for (auto i = 0Z; i < nlines_limited; i++) {
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
        const int element = globals::linelist.elementindex[lineindex];
        const int ion = globals::linelist.ionindex[lineindex];
        const double linelambda = 1e8 * CLIGHT / globals::linelist.nu[lineindex];
        // flux-weighted average radial velocity of emission in km/s
        double v_rad{NAN};
        if (mode == 0) {
          v_rad =
              traceemissionabsorption[i].emission_weightedvelocity_sum / traceemissionabsorption[i].energyemitted / 1e5;
        } else {
          v_rad = traceemissionabsorption[i].absorption_weightedvelocity_sum /
                  traceemissionabsorption[i].energyabsorbed / 1e5;
        }

        const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
        const auto lower_uniquelevelindex = globals::linelist.uniquelevelindex_lower[lineindex];
        const auto upper_uniquelevelindex = globals::linelist.uniquelevelindex_upper[lineindex];
        const int lower = lower_uniquelevelindex - ionuniquelevelindexstart;
        const int upper = upper_uniquelevelindex - ionuniquelevelindexstart;

        const double B_ul = globals::linelist.B_ul[lineindex];
        const double B_lu = globals::linelist.B_lu[lineindex];

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

        printlnlog("{:7.2e} ({:5.1f}%) {:4} {:9} {:5} {:5} {:8.1f} {:5} {:7.1f} {:7.1f} {:7.1e} {:7.1e}", encontrib,
                   100 * encontrib / totalenergy, get_atomicnumber(element), get_ionstage(element, ion), upper, lower,
                   globals::alltrans.coll_str[downtransid], static_cast<int>(globals::alltrans.forbidden[downtransid]),
                   linelambda, v_rad, B_lu, B_ul);
      } else {
        break;
      }
    }
    printlnlog("");
  }

  traceemissionabsorption.clear();
}

// number of different emission processes (bf and bb for each ion, and free-free)
auto get_proccount() -> int { return (2 * get_nelements() * get_max_nions()) + 1; }

auto columnindex_from_emissiontype(const int et) -> int {
  if (et >= 0) {
    // bb-emission
    const int element = globals::linelist.elementindex[et];
    const int ion = globals::linelist.ionindex[et];
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
    // no bf continua are in use, so a bf emission type should be impossible; count it in the free-free column
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

[[nodiscard]] auto get_absorption_spectrum_index(const ptrdiff_t nts, const ptrdiff_t nnu_abs) -> ptrdiff_t {
  const ptrdiff_t nelements = get_nelements();
  const ptrdiff_t max_nions = get_max_nions();
  return (nnu_abs * globals::ntimesteps * nelements * max_nions) + (nts * nelements * max_nions);
}

[[nodiscard]] inline auto get_timestep(const double time) -> int {
  assert_always(time >= globals::tmin);
  assert_always(time < globals::tmax);
  for (int nts = 0; nts < globals::ntimesteps; nts++) {
    const double tsend = (nts < (globals::ntimesteps - 1)) ? globals::timesteps[nts + 1].start : globals::tmax;
    if (time >= globals::timesteps[nts].start && time < tsend) {
      return nts;
    }
  }
  assert_always(false);  // could not find matching timestep

  return -1;
}

void write_partial_lightcurve_spectra_dirbin(const int nts, std::span<const Packet> packets,
                                             const bool do_emission_absorption, const int dirbin) {
  THREADLOCALONHOST std::vector<double> rpkt_light_curve_lum;
  THREADLOCALONHOST std::vector<double> rpkt_light_curve_lumcmf;
  THREADLOCALONHOST std::vector<double> gamma_light_curve_lum;
  THREADLOCALONHOST std::vector<double> gamma_light_curve_lumcmf;
  reserve_resize(rpkt_light_curve_lum, globals::ntimesteps);
  std::ranges::fill(rpkt_light_curve_lum, 0.);
  reserve_resize(rpkt_light_curve_lumcmf, globals::ntimesteps);
  std::ranges::fill(rpkt_light_curve_lumcmf, 0.);
  if constexpr (KEEP_ESCAPED_GAMMAS) {
    reserve_resize(gamma_light_curve_lum, globals::ntimesteps);
    std::ranges::fill(gamma_light_curve_lum, 0.);
    reserve_resize(gamma_light_curve_lumcmf, globals::ntimesteps);
    std::ranges::fill(gamma_light_curve_lumcmf, 0.);
  }

  init_spectra(rpkt_spectra_I, NU_MIN_R, NU_MAX_R, do_emission_absorption);

  MPI_Barrier_node();
#if defined REPRODUCIBLE && REPRODUCIBLE
  for (int node_rank = 0; node_rank < globals::node_nprocs; node_rank++) {
    // do one rank at a time to keep the results reproducible (instead of simultaneous atomic adds to shared memory)
#else
  {
    // all ranks on the node simultaneously contribute to the light curves and spectra in shared memory using
    // atomic operations
    const int node_rank = globals::rank_in_node;
#endif
    if (node_rank == globals::rank_in_node) {
      for (const auto& pkt : packets) {
        if (pkt.type == TYPE_ESCAPE) {
          if (pkt.escape_type == TYPE_RPKT) {
            add_to_lc_res(pkt, dirbin, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
            add_to_spec_res(pkt, dirbin, rpkt_spectra_I, nullptr, nullptr);
          } else if (KEEP_ESCAPED_GAMMAS && dirbin == -1 && pkt.escape_type == TYPE_GAMMA) {
            add_to_lc_res(pkt, dirbin, gamma_light_curve_lum, gamma_light_curve_lumcmf);
          }
        }
      }
    }
    MPI_Barrier_node();
  }

  const int numtimesteps = nts + 1;  // only produce spectra and light curves up to one past nts
  assert_always(numtimesteps <= globals::ntimesteps);

  MPI_Barrier_allranks();
  if (globals::rank_in_node == 0) {
    MPI_Allreduce_safe(rpkt_spectra_I.fluxalltimesteps, MPI_SUM, globals::mpi_comm_internode);
    if (rpkt_spectra_I.do_emission_absorption) {
      MPI_Allreduce_safe(rpkt_spectra_I.absorptionalltimesteps, MPI_SUM, globals::mpi_comm_internode);
      MPI_Allreduce_safe(rpkt_spectra_I.emissionalltimesteps, MPI_SUM, globals::mpi_comm_internode);
      MPI_Allreduce_safe(rpkt_spectra_I.trueemissionalltimesteps, MPI_SUM, globals::mpi_comm_internode);
    }
  }
  MPI_Allreduce_safe(rpkt_light_curve_lum, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce_safe(rpkt_light_curve_lumcmf, MPI_SUM, MPI_COMM_WORLD);
  if constexpr (KEEP_ESCAPED_GAMMAS) {
    MPI_Allreduce_safe(gamma_light_curve_lum, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce_safe(gamma_light_curve_lumcmf, MPI_SUM, MPI_COMM_WORLD);
  }
  MPI_Barrier_allranks();

  if (dirbin == -1) {
    write_light_curve("light_curve.out", rpkt_light_curve_lum, rpkt_light_curve_lumcmf, numtimesteps);
    if constexpr (KEEP_ESCAPED_GAMMAS) {
      write_light_curve("gamma_light_curve.out", gamma_light_curve_lum, gamma_light_curve_lumcmf, numtimesteps);
    }
    write_spectra("spec.out", "emission.out", "emissiontrue.out", "absorption.out", rpkt_spectra_I, numtimesteps);
  } else {
    if (globals::my_rank == 0 && !std::filesystem::exists(outdir_resfiles)) {
      std::filesystem::create_directory(outdir_resfiles);
    }
    MPI_Barrier_allranks();

    write_light_curve(std::format("{}light_curve_res_{:02d}.out", outdir_resfiles, dirbin), rpkt_light_curve_lum,
                      rpkt_light_curve_lumcmf, numtimesteps);
    write_spectra(std::format("{}spec_res_{:02d}.out", outdir_resfiles, dirbin),
                  std::format("{}emission_res_{:02d}.out", outdir_resfiles, dirbin),
                  std::format("{}emissiontrue_res_{:02d}.out", outdir_resfiles, dirbin),
                  std::format("{}absorption_res_{:02d}.out", outdir_resfiles, dirbin), rpkt_spectra_I, numtimesteps);
  }
  MPI_Barrier_allranks();
}

void write_spectrum_file(const std::string& spec_filename, const Spectra& spectra, const int numtimesteps) {
  auto spec_file = fstream_required(spec_filename, std::ios::out | std::ios::trunc);
  std::print(spec_file, "0 ");
  for (int p = 0; p < numtimesteps; p++) {
    std::print(spec_file, "{:g} ", globals::timesteps[p].mid / DAY);
  }
  std::println(spec_file, "");

  const auto ntimesteps_all = static_cast<ptrdiff_t>(globals::ntimesteps);
  for (auto nubin = 0Z; nubin < MNUBINS; nubin++) {
    std::print(spec_file, "{:g} ", (spectra.lower_freq[nubin] + (spectra.delta_freq[nubin] / 2)));

    for (auto nts = 0Z; nts < numtimesteps; nts++) {
      std::print(spec_file, "{:g} ", spectra.fluxalltimesteps[(nubin * ntimesteps_all) + nts]);
    }
    std::println(spec_file, "");
  }
}

// Write an emission-type spectrum (emission or true emission) with a line for each frequency bin of
// each timestep, holding one column per emission process (see get_proccount).
void write_emission_spectrum_file(const std::string& emission_filename,
                                  const std::span<const double> emission_alltimesteps, const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntimesteps);
  assert_always(!emission_filename.empty());
  auto emission_file = fstream_required(emission_filename, std::ios::out | std::ios::trunc);
  const auto ntimesteps_all = static_cast<ptrdiff_t>(globals::ntimesteps);
  const auto proccount = static_cast<ptrdiff_t>(get_proccount());
  for (auto nubin = 0Z; nubin < MNUBINS; nubin++) {
    for (auto nts = 0Z; nts < numtimesteps; nts++) {
      const auto emindex_nts_nubin = (nubin * ntimesteps_all * proccount) + (nts * proccount);
      for (int nproc = 0; nproc < proccount; nproc++) {
        std::print(emission_file, "{:g} ", emission_alltimesteps[emindex_nts_nubin + nproc]);
      }
      std::println(emission_file, "");
    }
  }
}

void write_absorption_spectrum_file(const std::string& absorption_filename, const Spectra& spectra,
                                    const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntimesteps);
  assert_always(!absorption_filename.empty());
  auto absorption_file = fstream_required(absorption_filename, std::ios::out | std::ios::trunc);
  const int ioncount = get_nelements() * get_max_nions();  // may be higher than the true included ion count
  for (auto nubin = 0Z; nubin < MNUBINS; nubin++) {
    for (auto nts = 0Z; nts < numtimesteps; nts++) {
      for (int i = 0; i < ioncount; i++) {
        std::print(absorption_file, "{:g} ",
                   spectra.absorptionalltimesteps[get_absorption_spectrum_index(nts, nubin) + i]);
      }
      std::println(absorption_file, "");
    }
  }
}

}  // anonymous namespace

void write_spectra(const std::string& spec_filename, const std::string& emission_filename,
                   const std::string& trueemission_filename, const std::string& absorption_filename,
                   const Spectra& spectra, const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntimesteps);

  if (TRACE_EMISSION_ABSORPTION_REGION_ON && spectra.do_emission_absorption && !traceemissionabsorption.empty() &&
      globals::my_rank == 0) {
    printout_tracemission_stats();
  }

  // only one rank should write each file. Try to choose different ranks on different nodes, if available
  const auto this_rank_writes_file = [](const int filenum) {
    return (filenum % globals::node_count == globals::node_id) &&
           (filenum % globals::node_nprocs == globals::rank_in_node);
  };

  if (this_rank_writes_file(1)) {
    write_spectrum_file(spec_filename, spectra, numtimesteps);
  }

  if (spectra.do_emission_absorption) {
    if (this_rank_writes_file(2)) {
      write_emission_spectrum_file(emission_filename, spectra.emissionalltimesteps.span(), numtimesteps);
    }

    if (this_rank_writes_file(3)) {
      write_emission_spectrum_file(trueemission_filename, spectra.trueemissionalltimesteps.span(), numtimesteps);
    }

    if (this_rank_writes_file(4)) {
      write_absorption_spectrum_file(absorption_filename, spectra, numtimesteps);
    }
  }
}

void write_specpol(const std::string& specpol_filename, const std::string& emission_filename,
                   const std::string& absorption_filename, const Spectra& spectra_I, const Spectra& spectra_Q,
                   const Spectra& spectra_U) {
  if (globals::my_rank != 0) {
    return;
  }
  auto specpol_file = fstream_required(specpol_filename, std::ios::out | std::ios::trunc);
  std::fstream emissionpol_file{};
  std::fstream absorptionpol_file{};

  const bool do_emission_absorption = spectra_I.do_emission_absorption;

  if (do_emission_absorption) {
    emissionpol_file = fstream_required(emission_filename, std::ios::out | std::ios::trunc);
    absorptionpol_file = fstream_required(absorption_filename, std::ios::out | std::ios::trunc);
    printlnlog("Writing {}, {}, and {}", specpol_filename, emission_filename, absorption_filename);
  } else {
    printlnlog("Writing {}", specpol_filename);
  }

  const int proccount = get_proccount();
  const int ioncount = get_nelements() * get_max_nions();  // may be higher than the true included ion count
  const auto ntimesteps = static_cast<ptrdiff_t>(globals::ntimesteps);

  std::print(specpol_file, "{:g}", 0.0);

  for (int l = 0; l < 3; l++) {
    for (int p = 0; p < globals::ntimesteps; p++) {
      std::print(specpol_file, " {:g}", globals::timesteps[p].mid / DAY);
    }
  }

  std::println(specpol_file, "");

  assert_always(std::ssize(spectra_I.delta_freq) == MNUBINS);
  assert_always(std::ssize(spectra_I.lower_freq) == MNUBINS);
  for (int nnu = 0; nnu < std::ssize(spectra_I.lower_freq); nnu++) {
    std::print(specpol_file, "{:g}", (spectra_I.lower_freq[nnu] + (spectra_I.delta_freq[nnu] / 2)));

    for (const auto* spec : {&spectra_I, &spectra_Q, &spectra_U}) {
      for (auto nts = 0Z; nts < ntimesteps; nts++) {
        std::print(specpol_file, " {:g}", spec->fluxalltimesteps[(nnu * ntimesteps) + nts]);

        if (do_emission_absorption) {
          for (int nproc = 0; nproc < proccount; nproc++) {
            const auto emindex = (nnu * ntimesteps * proccount) + (nts * proccount) + nproc;
            if (nproc > 0) {
              std::print(emissionpol_file, " ");
            }
            std::print(emissionpol_file, "{:g}", spec->emissionalltimesteps[emindex]);
          }
          std::println(emissionpol_file, "");

          for (int i = 0; i < ioncount; i++) {
            if (i > 0) {
              std::print(absorptionpol_file, " ");
            }
            std::print(absorptionpol_file, "{:g}",
                       spec->absorptionalltimesteps[get_absorption_spectrum_index(nts, nnu) + i]);
          }
          std::println(absorptionpol_file, "");
        }
      }
    }

    std::println(specpol_file, "");
  }
}

void init_spectrum_trace() {
  if (TRACE_EMISSION_ABSORPTION_REGION_ON) {
    traceemission_totalenergy = 0.;
    reserve_resize(traceemissionabsorption, globals::nlines);
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
void init_spectra(Spectra& spectra, const double nu_min, const double nu_max, const bool do_emission_absorption) {
  // setup the time and frequency bins using a logarithmic spacing in both t and nu

  assert_always(MNUBINS > 0);
  const double dlognu = (log(nu_max) - log(nu_min)) / MNUBINS;

  spectra.dlognu = dlognu;
  spectra.nu_min = nu_min;
  spectra.nu_max = nu_max;
  spectra.do_emission_absorption = do_emission_absorption;
  const bool print_memusage =
      (spectra.fluxalltimesteps.empty() || (do_emission_absorption && spectra.absorptionalltimesteps.empty()));

  for (auto nnu = 0Z; nnu < MNUBINS; nnu++) {
    spectra.lower_freq[nnu] = static_cast<float>(std::exp(log(nu_min) + (nnu * dlognu)));
    spectra.delta_freq[nnu] =
        static_cast<float>(std::exp(log(nu_min) + ((nnu + 1) * dlognu)) - spectra.lower_freq[nnu]);
  }

  if (spectra.fluxalltimesteps.empty()) {
    spectra.fluxalltimesteps = MPI_shared_array<double>(globals::ntimesteps * MNUBINS);
  }
  assert_always(std::ssize(spectra.fluxalltimesteps) == globals::ntimesteps * MNUBINS);
  std::ranges::fill(spectra.fluxalltimesteps, 0.0);
  MPI_Barrier_node();

  if (do_emission_absorption) {
    if (spectra.absorptionalltimesteps.empty()) {
      spectra.absorptionalltimesteps =
          MPI_shared_array<double>(globals::ntimesteps * MNUBINS * get_nelements() * get_max_nions());
    }
    assert_always(std::ssize(spectra.absorptionalltimesteps) ==
                  globals::ntimesteps * MNUBINS * get_nelements() * get_max_nions());
    std::ranges::fill(spectra.absorptionalltimesteps, 0.0);

    if (spectra.emissionalltimesteps.empty()) {
      spectra.emissionalltimesteps = MPI_shared_array<double>(globals::ntimesteps * MNUBINS * get_proccount());
    }
    assert_always(std::ssize(spectra.emissionalltimesteps) == globals::ntimesteps * MNUBINS * get_proccount());
    std::ranges::fill(spectra.emissionalltimesteps, 0.0);

    if (spectra.trueemissionalltimesteps.empty()) {
      spectra.trueemissionalltimesteps = MPI_shared_array<double>(globals::ntimesteps * MNUBINS * get_proccount());
    }
    assert_always(std::ssize(spectra.trueemissionalltimesteps) == globals::ntimesteps * MNUBINS * get_proccount());
    std::ranges::fill(spectra.trueemissionalltimesteps, 0.0);
  }
  MPI_Barrier_allranks();

  if (print_memusage) {
    printlnlog("[info] mem_usage: set of spectra{} occupy {:.3f} MB (node shared memory)",
               do_emission_absorption ? " (with emission/absorption tracing)" : "",
               spectra.mem_usage_bytes() / 1024. / 1024.);
  }
}

// Add a packet to the outgoing spectrum.
void add_to_spec_res(const Packet& pkt, const int dirbin, Spectra& spectra_I, Spectra* spectra_Q, Spectra* spectra_U) {
  if (dirbin != -1 && get_escapedirectionbin(pkt.dir) != dirbin) {
    return;  // do not add to the spectrum if the direction bin does not match
  }

  // Need to (1) decide which time bin to put it in and (2) which frequency bin.

  // specific angle bins contain fewer packets than the full sphere, so must be normalised to match
  const double nu_min = spectra_I.nu_min;
  const double nu_max = spectra_I.nu_max;
  const double dlognu = spectra_I.dlognu;
  const double t_arrive = pkt.escape_time - (dot(pkt.pos, pkt.dir) / CLIGHT_PROP);
  if (t_arrive > globals::tmin && t_arrive < globals::tmax && pkt.nu_rf > nu_min && pkt.nu_rf < nu_max) {
    const auto nts = get_timestep(t_arrive);

    // a binary search into freq_lower would probably be faster than this double logarithm
    const auto nnu = get_logbinindex(pkt.nu_rf, nu_min, dlognu, MNUBINS);

    const double solidanglefactor = (dirbin >= 0) ? MABINS : 1.;
    const double deltaE = pkt.e_rf / globals::timesteps[nts].width / spectra_I.delta_freq[nnu] / 4.e12 / PI / PARSEC /
                          PARSEC / globals::nprocs_exspec * solidanglefactor;

    const auto fluxindex = (nnu * static_cast<ptrdiff_t>(globals::ntimesteps)) + nts;
    atomicadd_always(spectra_I.fluxalltimesteps[fluxindex], deltaE);

    if (spectra_Q != nullptr) {
      atomicadd_always(spectra_Q->fluxalltimesteps[fluxindex], pkt.stokes_q * deltaE);
    }
    if (spectra_U != nullptr) {
      atomicadd_always(spectra_U->fluxalltimesteps[fluxindex], pkt.stokes_u * deltaE);
    }

    if (spectra_I.do_emission_absorption) {
      const auto proccount = get_proccount();

      const auto truenproc = columnindex_from_emissiontype(pkt.trueemissiontype);
      assert_always(truenproc < proccount);
      if (truenproc >= 0) {
        const auto emindex = (nnu * globals::ntimesteps * proccount) + (nts * proccount) + truenproc;
        atomicadd_always(spectra_I.trueemissionalltimesteps[emindex], deltaE);
      }

      const auto nproc = columnindex_from_emissiontype(pkt.emissiontype);
      assert_always(nproc < proccount);
      if (nproc >= 0) {  // -1 means EMTYPE_NOTSET
        const auto emindex = (nnu * globals::ntimesteps * proccount) + (nts * proccount) + nproc;
        atomicadd_always(spectra_I.emissionalltimesteps[emindex], deltaE);

        if (spectra_Q != nullptr && spectra_Q->do_emission_absorption) {
          atomicadd_always(spectra_Q->emissionalltimesteps[emindex], pkt.stokes_q * deltaE);
        }
        if (spectra_U != nullptr && spectra_U->do_emission_absorption) {
          atomicadd_always(spectra_U->emissionalltimesteps[emindex], pkt.stokes_u * deltaE);
        }
      }

      // traceemissionabsorption is only allocated by init_spectrum_trace(), which sn3d never calls, so
      // check that it has been set up before indexing into it
      if (TRACE_EMISSION_ABSORPTION_REGION_ON && (dirbin == -1) && !traceemissionabsorption.empty()) {
        const int et = pkt.trueemissiontype;
        if ((et >= 0) && (t_arrive >= traceemissabs_timemin && t_arrive <= traceemissabs_timemax) &&
            (pkt.nu_rf >= traceemissabs_nulower && pkt.nu_rf <= traceemissabs_nuupper))

        {
          traceemissionabsorption[et].energyemitted += deltaE;
          traceemissionabsorption[et].emission_weightedvelocity_sum +=
              vec_len(pkt.trueem_pos) / pkt.trueem_time * deltaE;
          traceemission_totalenergy += deltaE;
        }
      }

      if (pkt.absorptionfreq > nu_min && pkt.absorptionfreq < nu_max) {
        const auto nnu_abs = get_logbinindex(pkt.absorptionfreq, nu_min, dlognu, MNUBINS);
        const double deltaE_absorption = pkt.e_rf / globals::timesteps[nts].width / spectra_I.delta_freq[nnu_abs] /
                                         4.e12 / PI / PARSEC / PARSEC / globals::nprocs_exspec * solidanglefactor;
        const int at = pkt.absorptiontype;
        if (at >= 0) {
          // bb-absorption
          const int element = globals::linelist.elementindex[at];
          const int ion = globals::linelist.ionindex[at];
          const auto absindex = get_absorption_spectrum_index(nts, nnu_abs) + (element * get_max_nions()) + ion;
          atomicadd_always(spectra_I.absorptionalltimesteps[absindex], deltaE_absorption);

          if (spectra_Q != nullptr && spectra_Q->do_emission_absorption) {
            atomicadd_always(spectra_Q->absorptionalltimesteps[absindex], pkt.stokes_q * deltaE_absorption);
          }
          if (spectra_U != nullptr && spectra_U->do_emission_absorption) {
            atomicadd_always(spectra_U->absorptionalltimesteps[absindex], pkt.stokes_u * deltaE_absorption);
          }

          if ((TRACE_EMISSION_ABSORPTION_REGION_ON && !traceemissionabsorption.empty() &&
               t_arrive >= traceemissabs_timemin && t_arrive <= traceemissabs_timemax) &&
              ((dirbin == -1) && (pkt.nu_rf >= traceemissabs_nulower) && (pkt.nu_rf <= traceemissabs_nuupper))) {
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

void write_partial_lightcurve_spectra(const int nts, std::span<const Packet> pkts) {
  // this is called by sn3d (not exspec) when each rank has its own set of packets
  // in memory
  const bool simulation_complete = (nts >= globals::timestep_finish - 1);

  // the emission resolved spectra are slow to generate, and require a lot of memory
  const bool do_emission_absorption = WRITE_EMISSIONABSORPTION_SPEC_AT_END && simulation_complete;

  const bool multdimensional = grid::get_modelgridtype() != GridType::SPHERICAL1D;
  const int dirbinend = (multdimensional && simulation_complete) ? MABINS : 0;

  const auto time_func_start = std::chrono::steady_clock::now();

  for (int dirbin = -1; dirbin < dirbinend; dirbin++) {
    write_partial_lightcurve_spectra_dirbin(nts, pkts, do_emission_absorption, dirbin);
  }

  const auto duration_write_spectra =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - time_func_start).count();
  printlnlog("timestep {}: Saving light curves and {}spectra took {:.1f}s", nts,
             do_emission_absorption ? "emission/absorption " : "", duration_write_spectra);
}

void write_light_curve(const std::string& lc_filename, const std::span<const double> light_curve_lum,
                       const std::span<const double> light_curve_lumcmf, const int numtimesteps) {
  if (globals::node_id != 0 || globals::rank_in_node != 0) {
    return;
  }
  assert_always(numtimesteps <= globals::ntimesteps);

  auto lc_file = fstream_required(lc_filename, std::ios::out | std::ios::trunc);

  // UVOIR bolometric light curve
  for (int nts = 0; nts < numtimesteps; nts++) {
    std::println(lc_file, "{:g} {:g} {:g}", globals::timesteps[nts].mid / DAY, light_curve_lum[nts] / LSUN,
                 light_curve_lumcmf[nts] / LSUN);
  }
}

// add a packet to the outgoing light-curve.
void add_to_lc_res(const Packet& pkt, const int dirbin, std::span<double> light_curve_lum,
                   std::span<double> light_curve_lumcmf) {
  if (dirbin >= 0 && get_escapedirectionbin(pkt.dir) != dirbin) {
    return;
  }
  const double solidanglefactor = (dirbin >= 0) ? MABINS : 1.;
  // dirbin -1 means all full 4π angle average (no angle filtering)

  const double t_arrive = pkt.escape_time - (dot(pkt.pos, pkt.dir) / CLIGHT_PROP);
  if (t_arrive > globals::tmin && t_arrive < globals::tmax) {
    const int nts = get_timestep(t_arrive);
    atomicadd_always(light_curve_lum[nts],
                     pkt.e_rf / globals::timesteps[nts].width * solidanglefactor / globals::nprocs_exspec);
  }

  const double inverse_gamma = std::sqrt(1. - (globals::vmax * globals::vmax / CLIGHTSQUARED));

  // Now do the cmf light curve. Unlike t_arrive above, this has no light-travel-time term:
  // it is the escape time transformed (time-dilated) into the comoving frame.
  const double t_escape_cmf = pkt.escape_time * inverse_gamma;

  if (t_escape_cmf > globals::tmin && t_escape_cmf < globals::tmax) {
    const int nts = get_timestep(t_escape_cmf);
    atomicadd_always(light_curve_lumcmf[nts], pkt.e_cmf / globals::timesteps[nts].width * solidanglefactor /
                                                  globals::nprocs_exspec / inverse_gamma);
  }
}
