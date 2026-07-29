// Main program of the exspec post-processing tool: reads the packet files written by an sn3d
// run and bins the escaped packets into spectra and light curves for each observer direction
// bin, optionally with per-process emission and absorption contributions.

#include "exspec.h"

#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <span>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "mpi_logging.h"
#include "packet.h"
#include "sn3d.h"
#include "spectrum_lightcurve.h"
#include "version.h"

namespace {

void do_direction_bin(const int dirbin, const std::vector<std::vector<Packet>>& packets_by_rank) {
  THREADLOCALONHOST std::vector<double> rpkt_light_curve_lum;
  reserve_resize(rpkt_light_curve_lum, globals::ntimesteps);
  std::ranges::fill(rpkt_light_curve_lum, 0.);
  THREADLOCALONHOST std::vector<double> rpkt_light_curve_lumcmf;
  reserve_resize(rpkt_light_curve_lumcmf, globals::ntimesteps);
  std::ranges::fill(rpkt_light_curve_lumcmf, 0.);
  THREADLOCALONHOST std::vector<double> gamma_light_curve_lum;
  reserve_resize(gamma_light_curve_lum, globals::ntimesteps);
  std::ranges::fill(gamma_light_curve_lum, 0.);
  THREADLOCALONHOST std::vector<double> gamma_light_curve_lumcmf;
  reserve_resize(gamma_light_curve_lumcmf, globals::ntimesteps);
  std::ranges::fill(gamma_light_curve_lumcmf, 0.);

  THREADLOCALONHOST Spectra rpkt_spectra_I;
  // Set up the spectrum grid and initialise the bins to zero.
  init_spectra(rpkt_spectra_I, NU_MIN_R, NU_MAX_R, true);

  THREADLOCALONHOST Spectra rpkt_spectra_Q;
  THREADLOCALONHOST Spectra rpkt_spectra_U;

  if constexpr (POL_ON) {
    init_spectra(rpkt_spectra_Q, NU_MIN_R, NU_MAX_R, true);
    init_spectra(rpkt_spectra_U, NU_MIN_R, NU_MAX_R, true);
  }

  constexpr double nu_min_gamma = 0.05 * MEV / H;
  constexpr double nu_max_gamma = 4. * MEV / H;
  THREADLOCALONHOST Spectra gamma_spectra;
  init_spectra(gamma_spectra, nu_min_gamma, nu_max_gamma, false);
  assert_always(globals::nprocs_exspec > 0);
  for (int p = 0; p < globals::nprocs_exspec; p++) {
    const auto& pkts_thisrank = packets_by_rank[p];

    int nesc_gamma = 0;
    int nesc_rpkt = 0;
    for (const auto& pkt : pkts_thisrank) {
      if (pkt.type != TYPE_ESCAPE) {
        continue;
      }

      if (pkt.escape_type == TYPE_RPKT) {
        nesc_rpkt++;
        add_to_lc_res(pkt, dirbin, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
        add_to_spec_res(pkt, dirbin, rpkt_spectra_I, POL_ON ? &rpkt_spectra_Q : nullptr,
                        POL_ON ? &rpkt_spectra_U : nullptr);
      } else if (pkt.escape_type == TYPE_GAMMA) {
        nesc_gamma++;
        if (dirbin == -1) {
          add_to_lc_res(pkt, dirbin, gamma_light_curve_lum, gamma_light_curve_lumcmf);
          add_to_spec_res(pkt, dirbin, gamma_spectra, nullptr, nullptr);
        }
      }
    }
    if (dirbin == -1) {
      printlnlog("  rank {}: {} escaped r-packets and {} escaped gamma-pkts", p, nesc_rpkt, nesc_gamma);
    }
  }

  if (dirbin == -1) {
    // all directions integrated spectra and light curves
    write_light_curve("light_curve.out", rpkt_light_curve_lum, rpkt_light_curve_lumcmf, globals::ntimesteps);
    write_spectra("spec.out", "emission.out", "emissiontrue.out", "absorption.out", rpkt_spectra_I,
                  globals::ntimesteps);

    if constexpr (POL_ON) {
      write_specpol("specpol.out", "emissionpol.out", "absorptionpol.out", rpkt_spectra_I, rpkt_spectra_Q,
                    rpkt_spectra_U);
    }

    if constexpr (KEEP_ESCAPED_GAMMAS) {
      write_light_curve("gamma_light_curve.out", gamma_light_curve_lum, gamma_light_curve_lumcmf, globals::ntimesteps);
      write_spectra("gamma_spec.out", "", "", "", gamma_spectra, globals::ntimesteps);
    }

    // consistency check (log only): the frequency-integrated spectrum must reproduce the light curve, minus
    // the packets whose frequencies fall outside the spectrum's frequency range
    for (int nts = 0; nts < globals::ntimesteps; nts++) {
      double lum_from_spec = 0.;
      for (ptrdiff_t nnu = 0; nnu < MNUBINS; nnu++) {
        lum_from_spec += rpkt_spectra_I.fluxalltimesteps[(nnu * static_cast<ptrdiff_t>(globals::ntimesteps)) + nts] *
                         rpkt_spectra_I.delta_freq[nnu];
      }
      // undo the flux normalisation applied in add_to_spec_res() to get back to a luminosity
      lum_from_spec *= 4.e12 * PI * PARSEC * PARSEC;
      const double lum_lightcurve = rpkt_light_curve_lum[nts];
      if (lum_lightcurve > 0. && lum_from_spec > (lum_lightcurve * 1.001)) {
        printlnlog(
            "[warning] consistency check failed for timestep {}: frequency-integrated spec.out luminosity {:g} "
            "[erg/s] exceeds the light_curve.out luminosity {:g} [erg/s], but the spectrum's packets should be a "
            "subset of the light curve's packets",
            nts, lum_from_spec, lum_lightcurve);
      }
    }

    printlnlog("wrote the angle-averaged light curves and spectra");
  } else {
    // direction bin a
    // line-of-sight dependent spectra and light curves

    if (!std::filesystem::exists(outdir_resfiles)) {
      std::filesystem::create_directory(outdir_resfiles);
    }
    write_light_curve(std::format("{}light_curve_res_{:02d}.out", outdir_resfiles, dirbin), rpkt_light_curve_lum,
                      rpkt_light_curve_lumcmf, globals::ntimesteps);
    write_spectra(std::format("{}spec_res_{:02d}.out", outdir_resfiles, dirbin),
                  std::format("{}emission_res_{:02d}.out", outdir_resfiles, dirbin),
                  std::format("{}emissiontrue_res_{:02d}.out", outdir_resfiles, dirbin),
                  std::format("{}absorption_res_{:02d}.out", outdir_resfiles, dirbin), rpkt_spectra_I,
                  globals::ntimesteps);

    if constexpr (POL_ON) {
      write_specpol(std::format("{}specpol_res_{:02d}.out", outdir_resfiles, dirbin),
                    std::format("{}emissionpol_res_{:02d}.out", outdir_resfiles, dirbin),
                    std::format("{}absorptionpol_res_{:02d}.out", outdir_resfiles, dirbin), rpkt_spectra_I,
                    rpkt_spectra_Q, rpkt_spectra_U);
    }

    printlnlog("finished direction bin {} (highest bin is {})", dirbin, MABINS - 1);
  }
}

}  // anonymous namespace

auto main(int argc, char* argv[]) -> int {
  const auto sys_time_start = std::chrono::steady_clock::now();

  MPI_Init(&argc, &argv);

  globals::setup_mpi_vars();

  check_already_running();

  if (globals::my_rank == 0) {
    set_log_file("exspec.txt");
  }

  printlnlog("git branch: {}", GIT_BRANCH);

  printlnlog("git version: {}", GIT_VERSION);

  printlnlog("git status: {}", GIT_STATUS);

  printlnlog("exspec compiled at {} on {}", __TIME__, __DATE__);

#if defined TESTMODE && TESTMODE
  printlnlog("TESTMODE is ON");
#endif

  printlnlog("process id (pid): {}", getpid());
  printlnlog("MPI enabled:");
  printlnlog("  rank_global {} of [0..{}] in MPI_COMM_WORLD", globals::my_rank, globals::nprocs - 1);
  printlnlog("  rank_in_node {} of [0..{}] in node {} of [0..{}]", globals::rank_in_node, globals::node_nprocs - 1,
             globals::node_id, globals::node_count - 1);

  // single rank only for now
  assert_always(globals::my_rank == 0);
  assert_always(globals::nprocs == 1);

  // Read in parameters from input.txt
  read_parameterfile({});

  read_atomicdata();

  grid::read_ejecta_model();

  setup_timesteps();

  init_spectrum_trace();  // needed for TRACE_EMISSION_ABSORPTION_REGION_ON

  // nprocs_exspec is the number of rank output files to process with exspec
  // (not the number of ranks used to run exspec, which is always 1 for now)

  std::vector<std::vector<Packet>> packets_by_rank;
  reserve_resize(packets_by_rank, globals::nprocs_exspec);

  for (int p = 0; p < globals::nprocs_exspec; p++) {
    packets_by_rank[p] = read_text_packets(std::format("packets{:02d}_{:04d}.out", 0, p));
  }

  const int dirbinend = (grid::get_modelgridtype() == GridType::SPHERICAL1D) ? 0 : MABINS;
  // a is the escape direction angle bin
  for (int dirbin = -1; dirbin < dirbinend; dirbin++) {
    do_direction_bin(dirbin, packets_by_rank);
  }

  const auto exspec_duration = std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start).count();
  printlnlog("exspec finished (took {:.1f} seconds)", exspec_duration);

  MPI_Finalize();

  std::filesystem::remove("artis.pid");

  return 0;
}
