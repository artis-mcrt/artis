// Main program of the exspec post-processing tool: reads the packet files written by an sn3d
// run and bins the escaped packets into spectra and light curves for each observer direction
// bin, optionally with per-process emission and absorption contributions.
//
// sn3d writes the same files at its last requested timestep. exspec makes them again from the packet files,
// e.g. after a change of MNUBINS or of the frequency range.

#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <span>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "globals.h"
#include "grid.h"
#include "input.h"
#include "mpi_logging.h"
#include "packet.h"
#include "sn3d.h"
#include "spectrum_lightcurve.h"
#include "version.h"

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

  // nprocs_exspec is the number of rank output files to process with exspec
  // (not the number of ranks used to run exspec, which is always 1 for now)
  assert_always(globals::nprocs_exspec > 0);

  // one vector for each rank file. A file holds far fewer than MPKTS packets when KEEP_ESCAPED_GAMMAS is false.
  std::vector<std::vector<Packet>> packets_by_rank;
  reserve_resize(packets_by_rank, globals::nprocs_exspec);
  for (int rank = 0; rank < globals::nprocs_exspec; rank++) {
    packets_by_rank[rank] = read_text_packets(std::format("packets{:02d}_{:04d}.out", 0, rank));
    const auto escaped_rpkt_count = std::ranges::count_if(packets_by_rank[rank], [](const Packet& pkt) {
      return pkt.type == TYPE_ESCAPE && pkt.escape_type == TYPE_RPKT;
    });
    const auto escaped_gamma_count = std::ranges::count_if(packets_by_rank[rank], [](const Packet& pkt) {
      return pkt.type == TYPE_ESCAPE && pkt.escape_type == TYPE_GAMMA;
    });
    printlnlog("  rank {}: {} escaped r-packets and {} escaped gamma-pkts", rank, escaped_rpkt_count,
               escaped_gamma_count);
  }

  // the index of the last timestep also selects the emission, absorption, and direction bin files
  const std::vector<std::span<const Packet>> packet_spans_by_rank(packets_by_rank.begin(), packets_by_rank.end());
  write_light_curves_and_spectra(globals::ntimesteps - 1, packet_spans_by_rank);

  const auto exspec_duration = std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start).count();
  printlnlog("exspec finished (took {:.1f} seconds)", exspec_duration);

  MPI_Finalize();

  std::filesystem::remove("artis.pid");

  return 0;
}
