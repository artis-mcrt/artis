#include "exspec.h"

#include <unistd.h>

#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <new>
#include <span>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "packet.h"
#include "sn3d.h"
#include "spectrum_lightcurve.h"
#include "version.h"

std::fstream output_file;

namespace {

void do_angle_bin(const int a, std::span<Packet> pkts, bool load_allrank_packets, Spectra& rpkt_spectra,
                  Spectra& stokes_i, Spectra& stokes_q, Spectra& stokes_u, Spectra& gamma_spectra) {
  std::vector<double> rpkt_light_curve_lum(globals::ntimesteps, 0.);
  std::vector<double> rpkt_light_curve_lumcmf(globals::ntimesteps, 0.);
  std::vector<double> gamma_light_curve_lum(globals::ntimesteps, 0.);
  std::vector<double> gamma_light_curve_lumcmf(globals::ntimesteps, 0.);

  // Set up the spectrum grid and initialise the bins to zero.
  init_spectra(rpkt_spectra, NU_MIN_R, NU_MAX_R, true);

  if constexpr (POL_ON) {
    init_spectra(stokes_i, NU_MIN_R, NU_MAX_R, true);
    init_spectra(stokes_q, NU_MIN_R, NU_MAX_R, true);
    init_spectra(stokes_u, NU_MIN_R, NU_MAX_R, true);
  }

  const double nu_min_gamma = 0.05 * MEV / H;
  const double nu_max_gamma = 4. * MEV / H;
  init_spectra(gamma_spectra, nu_min_gamma, nu_max_gamma, false);
  assert_always(globals::nprocs_exspec > 0);
  for (int p = 0; p < globals::nprocs_exspec; p++) {
    auto pkts_start = pkts.subspan(load_allrank_packets ? p * globals::npkts : 0, globals::npkts);

    if (a == -1 || !load_allrank_packets) {
      auto pktfilename = std::format("packets{:02d}_{:04d}.out", 0, p);
      logprintlnfmt("Reading {} (file {} of {})", pktfilename, p + 1, globals::nprocs_exspec);

      if (std::filesystem::exists(pktfilename)) {
        read_packets(pktfilename, pkts_start);
      } else {
        logprintlnfmt("   WARNING {} does not exist - trying temp packets file at beginning of timestep {}...",
                      pktfilename, globals::timestep_initial);
        read_temp_packetsfile(globals::timestep_initial, p, pkts_start);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (p % globals::nprocs != globals::my_rank) {
      logprintlnfmt("skipping packets file {} {}", p + 1, globals::nprocs);
      continue;
    }

    int nesc_tot = 0;
    int nesc_gamma = 0;
    int nesc_rpkt = 0;
    for (int ii = 0; ii < globals::npkts; ii++) {
      if (pkts_start[ii].type == TYPE_ESCAPE) {
        nesc_tot++;
        if (pkts_start[ii].escape_type == TYPE_RPKT) {
          nesc_rpkt++;
          add_to_lc_res(pkts_start[ii], a, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
          add_to_spec_res(pkts_start[ii], a, rpkt_spectra, POL_ON ? &stokes_i : nullptr, POL_ON ? &stokes_q : nullptr,
                          POL_ON ? &stokes_u : nullptr);
        } else if (pkts_start[ii].escape_type == TYPE_GAMMA) {
          nesc_gamma++;
          if (a == -1) {
            add_to_lc_res(pkts_start[ii], a, gamma_light_curve_lum, gamma_light_curve_lumcmf);
            add_to_spec_res(pkts_start[ii], a, gamma_spectra, nullptr, nullptr, nullptr);
          }
        }
      }
    }
    if (a == -1 || !load_allrank_packets) {
      logprintlnfmt("  {} of {} packets escaped ({} gamma-pkts and {} r-pkts)", nesc_tot, globals::npkts, nesc_gamma,
                    nesc_rpkt);
    }
  }

  if (a == -1) {
    // angle-averaged spectra and light curves
    write_light_curve("light_curve.out", -1, rpkt_light_curve_lum, rpkt_light_curve_lumcmf, globals::ntimesteps);
    write_light_curve("gamma_light_curve.out", -1, gamma_light_curve_lum, gamma_light_curve_lumcmf,
                      globals::ntimesteps);

    write_spectra("spec.out", "emission.out", "emissiontrue.out", "absorption.out", rpkt_spectra, globals::ntimesteps);

    if constexpr (POL_ON) {
      write_specpol("specpol.out", "emissionpol.out", "absorptionpol.out", &stokes_i, &stokes_q, &stokes_u);
    }

    write_spectra("gamma_spec.out", "", "", "", gamma_spectra, globals::ntimesteps);

    logprintlnfmt("finished angle-averaged stuff");
  } else {
    // direction bin a
    // line-of-sight dependent spectra and light curves

    if (!std::filesystem::exists(outdir_resfiles)) {
      std::filesystem::create_directory(outdir_resfiles);
    }
    write_light_curve(std::format("{}light_curve_res_{:02d}.out", outdir_resfiles, a), a, rpkt_light_curve_lum,
                      rpkt_light_curve_lumcmf, globals::ntimesteps);
    write_spectra(std::format("{}spec_res_{:02d}.out", outdir_resfiles, a),
                  std::format("{}emission_res_{:02d}.out", outdir_resfiles, a),
                  std::format("{}emissiontrue_res_{:02d}.out", outdir_resfiles, a),
                  std::format("{}absorption_res_{:02d}.out", outdir_resfiles, a), rpkt_spectra, globals::ntimesteps);

    if constexpr (POL_ON) {
      write_specpol(std::format("{}specpol_res_{:02d}.out", outdir_resfiles, a),
                    std::format("{}emissionpol_res_{:02d}.out", outdir_resfiles, a),
                    std::format("{}absorptionpol_res_{:02d}.out", outdir_resfiles, a), &stokes_i, &stokes_q, &stokes_u);
    }

    logprintlnfmt("Did {} of {} angle bins.", a + 1, MABINS);
  }
}

}  // anonymous namespace

auto main(int argc, char* argv[]) -> int {
  const auto sys_time_start = std::time(nullptr);

  MPI_Init(&argc, &argv);

  globals::setup_mpi_vars();

  check_already_running();

  if (globals::my_rank == 0) {
    output_file = fstream_required("exspec.txt", std::ios::out | std::ios::trunc);
  }

  logprintlnfmt("git branch {}", GIT_BRANCH);

  logprintlnfmt("git version: {}", GIT_VERSION);

  logprintlnfmt("git status {}", GIT_STATUS);

  logprintlnfmt("exspec compiled at {} on {}", __TIME__, __DATE__);

#if defined TESTMODE && TESTMODE
  logprintlnfmt("TESTMODE is ON");
#endif

  logprintlnfmt("process id (pid): {}", getpid());
  logprintlnfmt("MPI enabled:");
  logprintlnfmt("  rank_global {} of [0..{}] in MPI_COMM_WORLD", globals::my_rank, globals::nprocs - 1);
  logprintlnfmt("  rank_in_node {} of [0..{}] in node {} of [0..{}]", globals::rank_in_node, globals::node_nprocs - 1,
                globals::node_id, globals::node_count - 1);

  // single rank only for now
  assert_always(globals::my_rank == 0);
  assert_always(globals::nprocs == 1);

  // Get input stuff
  input(globals::my_rank);

  // nprocs_exspec is the number of rank output files to process with expec
  // however, we might be running exspec with 1 or just a few ranks

  std::vector<Packet> pkts;
  bool load_allrank_packets = false;
  const ptrdiff_t npkts_allranks = static_cast<ptrdiff_t>(globals::nprocs_exspec) * globals::npkts;
  try {
    // try to allocate memory for all packets from all ranks
    resize_exactly(pkts, npkts_allranks);
    logprintlnfmt(
        "mem_usage: loading {} packets from each {} processes simultaneously (total {} packets, {:.1f} MB memory)",
        globals::npkts, globals::nprocs_exspec, npkts_allranks, npkts_allranks * sizeof(Packet) / 1024. / 1024.);
    load_allrank_packets = true;
  } catch (const std::bad_alloc& e) {
    // if we can't allocate memory for all packets, try to allocate memory for just one rank
    load_allrank_packets = false;
    logprintlnfmt("mem_usage: malloc failed to allocate memory for all packets");
    logprintlnfmt(
        "mem_usage: loading {} packets from each of {} processes sequentially (total {} packets, {:.1f} MB memory)",
        globals::npkts, globals::nprocs_exspec, npkts_allranks, npkts_allranks * sizeof(Packet) / 1024. / 1024.);
    resize_exactly(pkts, globals::npkts);
  }

  init_spectrum_trace();  // needed for TRACE_EMISSION_ABSORPTION_REGION_ON

  Spectra rpkt_spectra;

  Spectra stokes_i;
  Spectra stokes_q;
  Spectra stokes_u;

  Spectra gamma_spectra;

  setup_timesteps();

  const int dirbinend = (grid::get_model_type() == GridType::SPHERICAL1D) ? 0 : MABINS;
  // a is the escape direction angle bin
  for (int dirbin = -1; dirbin < dirbinend; dirbin++) {
    do_angle_bin(dirbin, pkts, load_allrank_packets, rpkt_spectra, stokes_i, stokes_q, stokes_u, gamma_spectra);
  }

  decay::cleanup();
  logprintlnfmt("exspec finished at {} (tstart + {} seconds)", std::time(nullptr), std::time(nullptr) - sys_time_start);

  globals::mpi_finalized = true;
  MPI_Finalize();

  if (std::filesystem::exists("artis.pid")) {
    std::filesystem::remove("artis.pid");
  }

  return 0;
}
