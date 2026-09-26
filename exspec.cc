// Main program of the exspec post-processing tool: reads the packet files written by an sn3d
// run and bins the escaped packets into spectra and light curves for each observer direction
// bin, optionally with per-process emission and absorption contributions.
//
// sn3d writes the same files at its last requested timestep. exspec makes them again from the packet files,
// e.g. after a change of MNUBINS or of the frequency range.

#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <functional>
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

  // the log lines of the other ranks go nowhere
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

  // Read in parameters from input.txt
  read_parameterfile({});

  // nprocs_exspec is the number of packet files, one for each sn3d rank
  assert_always(globals::nprocs_exspec > 0);
  if (globals::nprocs > globals::nprocs_exspec) {
    fatal_crash("exspec runs with {} ranks but there are only {} packet files", globals::nprocs,
                globals::nprocs_exspec);
  }

  read_atomicdata();

  grid::read_ejecta_model();

  setup_timesteps();

  // each exspec rank reads a contiguous block of the packet files
  const auto [firstfile, nfiles] = get_range_chunk(globals::nprocs_exspec, globals::nprocs, globals::my_rank);
  printlnlog("{} packet files, read by {} exspec ranks", globals::nprocs_exspec, globals::nprocs);

  // one vector for each packet file
  std::vector<std::vector<Packet>> packets_by_file;
  packets_by_file.reserve(nfiles);
  // the counts of each packet file, summed over the exspec ranks so that rank 0 can log every file
  std::vector<std::int64_t> packet_count(globals::nprocs_exspec);
  std::vector<std::int64_t> escaped_rpkt_count(globals::nprocs_exspec);
  std::vector<std::int64_t> escaped_gamma_count(globals::nprocs_exspec);
  for (auto sn3d_rank = firstfile; sn3d_rank < firstfile + nfiles; sn3d_rank++) {
    packets_by_file.push_back(read_text_packets(std::format("packets/packets{:02d}_{:04d}.out", 0, sn3d_rank)));
    const auto& packets = packets_by_file.back();
    packet_count[sn3d_rank] = std::ssize(packets);
    escaped_rpkt_count[sn3d_rank] = std::ranges::count_if(
        packets, [](const Packet& pkt) { return pkt.type == TYPE_ESCAPE && pkt.escape_type == TYPE_RPKT; });
    escaped_gamma_count[sn3d_rank] = std::ranges::count_if(
        packets, [](const Packet& pkt) { return pkt.type == TYPE_ESCAPE && pkt.escape_type == TYPE_GAMMA; });
  }
  MPI_Allreduce_safe(packet_count, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce_safe(escaped_rpkt_count, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce_safe(escaped_gamma_count, MPI_SUM, MPI_COMM_WORLD);
  for (auto sn3d_rank = 0Z; sn3d_rank < globals::nprocs_exspec; sn3d_rank++) {
    printlnlog("  packets{:02d}_{:04d}.out: {} packets, {} escaped r-packets and {} escaped gamma-pkts", 0, sn3d_rank,
               packet_count[sn3d_rank], escaped_rpkt_count[sn3d_rank], escaped_gamma_count[sn3d_rank]);
  }
  printlnlog("total: {} packets, {} escaped r-packets and {} escaped gamma-pkts",
             std::ranges::fold_left(packet_count, 0Z, std::plus{}),
             std::ranges::fold_left(escaped_rpkt_count, 0Z, std::plus{}),
             std::ranges::fold_left(escaped_gamma_count, 0Z, std::plus{}));

  // the index of the last timestep also selects the emission, absorption, and direction bin files
  const std::vector<std::span<const Packet>> packet_spans_by_file(packets_by_file.begin(), packets_by_file.end());
  write_light_curves_and_spectra(globals::ntimesteps - 1, packet_spans_by_file);

  const auto exspec_duration = std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start).count();
  printlnlog("exspec finished (took {:.1f} seconds)", exspec_duration);

  MPI_Finalize();

  std::filesystem::remove("artis.pid");

  return 0;
}
