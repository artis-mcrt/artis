// Main program of the exspec post-processing tool: reads the packet files written by an sn3d
// run and bins the escaped packets into spectra and light curves for each observer direction
// bin, optionally with per-process emission and absorption contributions.
//
// sn3d writes the same files at its last requested timestep. exspec makes them again from the packet files,
// e.g. after a change of MNUBINS or of the frequency range.

#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <span>
#include <string>
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

  // rank 0 keeps the name exspec.txt, which the scripts and the workflows read
  set_log_file(globals::my_rank == 0 ? std::string("exspec.txt") : std::format("exspec_{}.txt", globals::my_rank));

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

  read_atomicdata();

  grid::read_ejecta_model();

  setup_timesteps();

  // nprocs_exspec is the number of packet files, one for each sn3d rank
  assert_always(globals::nprocs_exspec > 0);
  if (globals::nprocs > globals::nprocs_exspec) {
    fatal_crash("exspec runs with {} ranks but there are only {} packet files. Use at most {} ranks.", globals::nprocs,
                globals::nprocs_exspec, globals::nprocs_exspec);
  }

  // Each exspec rank reads a contiguous block of the packet files. With one node, the blocks keep the order of
  // the packets equal to a run with a single rank.
  const auto [firstfile, nfiles] = get_range_chunk(globals::nprocs_exspec, globals::nprocs, globals::my_rank);
  printlnlog("this rank reads {} of the {} packet files, from packets{:02d}_{:04d}.out", nfiles, globals::nprocs_exspec,
             0, firstfile);

  // one vector for each packet file. A file holds far fewer than MPKTS packets when KEEP_ESCAPED_GAMMAS is false.
  std::vector<std::vector<Packet>> packets_by_rank;
  reserve_resize(packets_by_rank, nfiles);
  for (ptrdiff_t i = 0; i < nfiles; i++) {
    const auto rank = firstfile + i;
    packets_by_rank[i] = read_text_packets(std::format("packets{:02d}_{:04d}.out", 0, rank));
    const auto escaped_rpkt_count = std::ranges::count_if(
        packets_by_rank[i], [](const Packet& pkt) { return pkt.type == TYPE_ESCAPE && pkt.escape_type == TYPE_RPKT; });
    const auto escaped_gamma_count = std::ranges::count_if(
        packets_by_rank[i], [](const Packet& pkt) { return pkt.type == TYPE_ESCAPE && pkt.escape_type == TYPE_GAMMA; });
    printlnlog("  rank {}: {} escaped r-packets and {} escaped gamma-pkts", rank, escaped_rpkt_count,
               escaped_gamma_count);
  }

  // the index of the last timestep also selects the emission, absorption, and direction bin files
  const std::vector<std::span<const Packet>> packet_spans_by_rank(packets_by_rank.begin(), packets_by_rank.end());
  write_light_curves_and_spectra(globals::ntimesteps - 1, packet_spans_by_rank);

  const auto exspec_duration = std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start).count();
  printlnlog("exspec finished (took {:.1f} seconds)", exspec_duration);

  MPI_Finalize();

  if (globals::my_rank == 0) {
    std::filesystem::remove("artis.pid");
  }

  return 0;
}
