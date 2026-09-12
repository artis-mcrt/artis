// Main program of the exspec post-processing tool: reads the packet files written by an sn3d
// run and bins the escaped packets into spectra and light curves for each observer direction
// bin, optionally with per-process emission and absorption contributions.
//
// sn3d writes the same files at its last requested timestep. exspec makes them again from the packet files,
// e.g. after a change of MNUBINS or of the frequency range.

#include "exspec.h"

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

#include "artisoptions.h"
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

  // the packets of all ranks in the order of the rank files, so that the sums are reproducible
  std::vector<Packet> packets;
  packets.reserve(static_cast<size_t>(globals::nprocs_exspec) * MPKTS);
  for (int p = 0; p < globals::nprocs_exspec; p++) {
    const auto packets_thisrank = read_text_packets(std::format("packets{:02d}_{:04d}.out", 0, p));
    const auto nesc_rpkt = std::ranges::count_if(
        packets_thisrank, [](const Packet& pkt) { return pkt.type == TYPE_ESCAPE && pkt.escape_type == TYPE_RPKT; });
    const auto nesc_gamma = std::ranges::count_if(
        packets_thisrank, [](const Packet& pkt) { return pkt.type == TYPE_ESCAPE && pkt.escape_type == TYPE_GAMMA; });
    printlnlog("  rank {}: {} escaped r-packets and {} escaped gamma-pkts", p, nesc_rpkt, nesc_gamma);
    packets.insert(packets.end(), packets_thisrank.begin(), packets_thisrank.end());
  }

  // the index of the last timestep also selects the emission, absorption, and direction bin files
  write_partial_lightcurve_spectra(globals::ntimesteps - 1, packets);

  const auto exspec_duration = std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start).count();
  printlnlog("exspec finished (took {:.1f} seconds)", exspec_duration);

  MPI_Finalize();

  std::filesystem::remove("artis.pid");

  return 0;
}
