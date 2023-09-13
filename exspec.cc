#include "exspec.h"

#include <unistd.h>

#include <cstdio>
#include <filesystem>
#include <memory>

#include "decay.h"
#include "grid.h"
#include "input.h"
#include "light_curve.h"
#include "sn3d.h"
#include "spectrum.h"
#include "version.h"

// threadprivate variables
FILE *output_file = nullptr;
int tid = 0;
bool use_cellhist = false;
std::mt19937 stdrng;
gsl_integration_workspace *gslworkspace = nullptr;

auto main(int argc, char *argv[]) -> int {
  const time_t sys_time_start = time(nullptr);

#ifdef MPI_ON
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &globals::rank_global);
  MPI_Comm_size(MPI_COMM_WORLD, &globals::nprocs);

  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globals::rank_global, MPI_INFO_NULL,
                      &globals::mpi_comm_node);
  // get the local rank within this node
  MPI_Comm_rank(globals::mpi_comm_node, &globals::rank_in_node);
  // get the number of ranks on the node
  MPI_Comm_size(globals::mpi_comm_node, &globals::node_nprocs);
  MPI_Barrier(MPI_COMM_WORLD);

  // make an inter-node communicator (using local rank as the key for group membership)
  MPI_Comm_split(MPI_COMM_WORLD, globals::rank_in_node, globals::rank_global, &globals::mpi_comm_internode);

  // take the node id from the local rank 0 (node master) and broadcast it
  if (globals::rank_in_node == 0) {
    MPI_Comm_rank(globals::mpi_comm_internode, &globals::node_id);
    MPI_Comm_size(globals::mpi_comm_internode, &globals::node_count);
  }

  MPI_Bcast(&globals::node_id, 1, MPI_INT, 0, globals::mpi_comm_node);
  MPI_Bcast(&globals::node_count, 1, MPI_INT, 0, globals::mpi_comm_node);
  MPI_Barrier(MPI_COMM_WORLD);
#else
  globals::rank_global = 0;
  globals::nprocs = 1;
  globals::rank_in_node = 0;
  globals::node_nprocs = 1;
  globals::node_id = 0;
  globals::node_count = 0;
#endif
  char filename[MAXFILENAMELENGTH];

  globals::startofline = std::make_unique<bool[]>(get_max_threads());
  if (globals::rank_global == 0) {
    check_already_running();
  }

  // make sure rank 0 checked for a pid file before we proceed
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (globals::rank_global == 0) {
    snprintf(filename, MAXFILENAMELENGTH, "exspec.txt");
    output_file = fopen_required(filename, "w");
    setvbuf(output_file, nullptr, _IOLBF, 1);
  }

  printout("git branch %s\n", GIT_BRANCH);

  printout("git version: %s\n", GIT_VERSION);

  printout("git status %s\n", GIT_STATUS);

  // printout("Hash of most recent commit: %s\n",GIT_HASH);
  printout("exspec compiled at %s on %s\n", __TIME__, __DATE__);

#if defined TESTMODE && TESTMODE
  printout("TESTMODE is ON\n");
#endif

#ifdef MPI_ON
  printout("process id (pid): %d\n", getpid());
  printout("MPI enabled:\n");
  printout("  rank %d of [0..%d] in MPI_COMM_WORLD\n", globals::rank_global, globals::nprocs - 1);
  printout("  node %d of [0..%d]\n", globals::node_id, globals::node_count - 1);
  printout("  rank %d of [0..%d] within this node (MPI_COMM_WORLD_SHARED)\n", globals::rank_in_node,
           globals::node_nprocs - 1);
#else
  printout("MPI is disabled in this build\n");
#endif

  // single rank only for now
  assert_always(globals::rank_global == 0);
  assert_always(globals::nprocs == 1);

  printout("Begining exspec.\n");

  /// Get input stuff
  printout("time before input %ld\n", time(nullptr));
  input(globals::rank_global);
  printout("time after input %ld\n", time(nullptr));
  // nprocs_exspec is the number of rank output files to process with expec
  // however, we might be running exspec with 1 or just a few ranks
  globals::nprocs = globals::nprocs_exspec;

  constexpr double maxpktmem_mb = 16000;
  bool load_allrank_packets = false;

  if ((globals::nprocs_exspec * globals::npkts * sizeof(struct packet) / 1024. / 1024.) < maxpktmem_mb) {
    printout(
        "mem_usage: loading packets from all %d processes simultaneously (total %d packets, %.1f MB memory is within "
        "limit of %.1f MB)\n",
        globals::nprocs_exspec, globals::nprocs_exspec * globals::npkts,
        globals::nprocs_exspec * globals::npkts * sizeof(struct packet) / 1024. / 1024., maxpktmem_mb);
    load_allrank_packets = true;
  } else {
    printout(
        "mem_usage: loading packets from each of %d processes sequentially (total %d packets, %.1f MB memory would be "
        "above limit of %.1f MB)\n",
        globals::nprocs_exspec, globals::nprocs_exspec * globals::npkts,
        globals::nprocs_exspec * globals::npkts * sizeof(struct packet) / 1024. / 1024., maxpktmem_mb);
    load_allrank_packets = false;
  }

  const int npkts_loaded = load_allrank_packets ? globals::nprocs_exspec * globals::npkts : globals::npkts;
  auto *pkts = static_cast<struct packet *>(malloc(npkts_loaded * sizeof(struct packet)));

  init_spectrum_trace();  // needed for TRACE_EMISSION_ABSORPTION_REGION_ON

  auto rpkt_spectra = alloc_spectra(globals::do_emission_res);

  std::unique_ptr<struct spec> stokes_i = nullptr;
  std::unique_ptr<struct spec> stokes_q = nullptr;
  std::unique_ptr<struct spec> stokes_u = nullptr;

  if constexpr (POL_ON) {
    stokes_i = alloc_spectra(globals::do_emission_res);
    stokes_q = alloc_spectra(globals::do_emission_res);
    stokes_u = alloc_spectra(globals::do_emission_res);
  }

  auto gamma_spectra = alloc_spectra(false);

  /// Initialise the grid. Call routine that sets up the initial positions
  /// and sizes of the grid cells.
  // grid_init();
  time_init();

  const int amax = ((grid::get_model_type() == grid::RHO_1D_READ)) ? 0 : MABINS;
  // a is the escape direction angle bin
  for (int a = -1; a < amax; a++) {
    /// Set up the light curve grid and initialise the bins to zero.
    std::vector<double> rpkt_light_curve_lum(globals::ntstep, 0.);
    std::vector<double> rpkt_light_curve_lumcmf(globals::ntstep, 0.);
    std::vector<double> gamma_light_curve_lum(globals::ntstep, 0.);
    std::vector<double> gamma_light_curve_lumcmf(globals::ntstep, 0.);
    /// Set up the spectrum grid and initialise the bins to zero.

    init_spectra(*rpkt_spectra, NU_MIN_R, NU_MAX_R, globals::do_emission_res);

    if constexpr (POL_ON) {
      init_spectra(*stokes_i, NU_MIN_R, NU_MAX_R, globals::do_emission_res);
      init_spectra(*stokes_q, NU_MIN_R, NU_MAX_R, globals::do_emission_res);
      init_spectra(*stokes_u, NU_MIN_R, NU_MAX_R, globals::do_emission_res);
    }

    const double nu_min_gamma = 0.05 * MEV / H;
    const double nu_max_gamma = 4. * MEV / H;
    init_spectra(*gamma_spectra, nu_min_gamma, nu_max_gamma, false);

    for (int p = 0; p < globals::nprocs_exspec; p++) {
      struct packet *pkts_start = load_allrank_packets ? &pkts[p * globals::npkts] : pkts;

      if (a == -1 || !load_allrank_packets) {
        char pktfilename[MAXFILENAMELENGTH];
        snprintf(pktfilename, MAXFILENAMELENGTH, "packets%.2d_%.4d.out", 0, p);
        printout("reading %s (file %d of %d)\n", pktfilename, p + 1, globals::nprocs_exspec);

        if (access(pktfilename, F_OK) == 0) {
          read_packets(pktfilename, pkts_start);
        } else {
          printout("   WARNING %s does not exist - trying temp packets file at beginning of timestep %d...\n   ",
                   pktfilename, globals::itstep);
          read_temp_packetsfile(globals::itstep, p, pkts_start);
        }
      }

      int nesc_tot = 0;
      int nesc_gamma = 0;
      int nesc_rpkt = 0;
      for (int ii = 0; ii < globals::npkts; ii++) {
        // printout("packet %d escape_type %d type %d", ii, pkts[ii].escape_type, pkts[ii].type);
        if (pkts_start[ii].type == TYPE_ESCAPE) {
          nesc_tot++;
          if (pkts_start[ii].escape_type == TYPE_RPKT) {
            nesc_rpkt++;
            add_to_lc_res(&pkts_start[ii], a, rpkt_light_curve_lum.data(), rpkt_light_curve_lumcmf.data());
            add_to_spec_res(&pkts_start[ii], a, *rpkt_spectra, stokes_i.get(), stokes_q.get(), stokes_u.get());
          } else if (pkts_start[ii].escape_type == TYPE_GAMMA) {
            nesc_gamma++;
            if (a == -1) {
              add_to_lc_res(&pkts_start[ii], a, gamma_light_curve_lum.data(), gamma_light_curve_lumcmf.data());
              add_to_spec_res(&pkts_start[ii], a, *gamma_spectra, nullptr, nullptr, nullptr);
            }
          }
        }
      }
      if (a == -1 || !load_allrank_packets) {
        printout("  %d of %d packets escaped (%d gamma-pkts and %d r-pkts)\n", nesc_tot, globals::npkts, nesc_gamma,
                 nesc_rpkt);
      }
    }

    if (a == -1) {
      // angle-averaged spectra and light curves
      write_light_curve("light_curve.out", -1, rpkt_light_curve_lum.data(), rpkt_light_curve_lumcmf.data(),
                        globals::ntstep);
      write_light_curve("gamma_light_curve.out", -1, gamma_light_curve_lum.data(), gamma_light_curve_lumcmf.data(),
                        globals::ntstep);

      write_spectrum("spec.out", "emission.out", "emissiontrue.out", "absorption.out", *rpkt_spectra, globals::ntstep);

      if constexpr (POL_ON) {
        write_specpol("specpol.out", "emissionpol.out", "absorptionpol.out", stokes_i.get(), stokes_q.get(),
                      stokes_u.get());
      }

      write_spectrum("gamma_spec.out", nullptr, nullptr, nullptr, *gamma_spectra, globals::ntstep);

      printout("finished angle-averaged stuff\n");
    } else {
      // direction bin a
      // line-of-sight dependent spectra and light curves

      char lc_filename[MAXFILENAMELENGTH] = "";
      snprintf(lc_filename, MAXFILENAMELENGTH, "light_curve_res_%.2d.out", a);

      char spec_filename[MAXFILENAMELENGTH] = "";
      snprintf(spec_filename, MAXFILENAMELENGTH, "spec_res_%.2d.out", a);

      char emission_filename[MAXFILENAMELENGTH] = "";
      snprintf(emission_filename, MAXFILENAMELENGTH, "emission_res_%.2d.out", a);

      char trueemission_filename[MAXFILENAMELENGTH] = "";
      snprintf(trueemission_filename, MAXFILENAMELENGTH, "emissiontrue_res_%.2d.out", a);

      char absorption_filename[MAXFILENAMELENGTH] = "";
      snprintf(absorption_filename, MAXFILENAMELENGTH, "absorption_res_%.2d.out", a);

      write_light_curve(lc_filename, a, rpkt_light_curve_lum.data(), rpkt_light_curve_lumcmf.data(), globals::ntstep);
      write_spectrum(spec_filename, emission_filename, trueemission_filename, absorption_filename, *rpkt_spectra,
                     globals::ntstep);

      if constexpr (POL_ON) {
        char specpol_filename[MAXFILENAMELENGTH] = "";
        snprintf(specpol_filename, MAXFILENAMELENGTH, "specpol_res_%.2d.out", a);

        char emissionpol_filename[MAXFILENAMELENGTH] = "";
        snprintf(emissionpol_filename, MAXFILENAMELENGTH, "emissionpol_res_%.2d.out", a);

        char absorptionpol_filename[MAXFILENAMELENGTH] = "";
        snprintf(absorptionpol_filename, MAXFILENAMELENGTH, "absorptionpol_res_%.2d.out", a);

        write_specpol(specpol_filename, emissionpol_filename, absorptionpol_filename, stokes_i.get(), stokes_q.get(),
                      stokes_u.get());
      }

      printout("Did %d of %d angle bins.\n", a + 1, MABINS);
    }
  }

  rpkt_spectra.reset();
  stokes_i.reset();
  stokes_q.reset();
  stokes_u.reset();
  gamma_spectra.reset();

  free(pkts);
  decay::cleanup();

  printout("exspec finished at %ld (tstart + %ld seconds)\n", time(nullptr), time(nullptr) - sys_time_start);
  fclose(output_file);

#ifdef MPI_ON
  MPI_Finalize();
#endif

  if (std::filesystem::exists("artis.pid")) {
    std::filesystem::remove("artis.pid");
  }

  return 0;
}
