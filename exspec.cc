/* 2007-10-30 -- MK
   Non-grey treatment of UVOIR opacity as opacity_case 4 added.
   Still not fully commented.
   Comments are marked by ///  Deactivated code by // */
/* 2007-01-17 -- MK
   Several minor modifications (some marked in the code with //MK), these include
     - global printout() routine (located in sn3d.c)
     - opacity_cases 2 and 3 added (changes in grid_init.c and update_grid.c,
       original opacity stuff was moved there from input.c) */
/* This is a code copied from Lucy 2004 paper on t-dependent supernova
   explosions. */

#include "exspec.h"

#include <unistd.h>

#include <cstdio>

#include "decay.h"
#include "grid.h"
#include "input.h"
#include "light_curve.h"
#include "sn3d.h"
#include "spectrum.h"

const bool do_exspec = true;

// threadprivate variables
FILE *output_file = nullptr;
int tid = 0;
bool use_cellhist = false;
bool neutral_flag = false;
gsl_rng *rng = nullptr;
std::mt19937_64 *stdrng = nullptr;
gsl_integration_workspace *gslworkspace = nullptr;

int main(int argc, char **argv) {
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
  char filename[128];

  globals::startofline = std::make_unique<bool[]>(get_max_threads());
  if (globals::rank_global == 0) {
    check_already_running();

    snprintf(filename, 128, "exspec.txt");
    output_file = fopen_required(filename, "w");
    setvbuf(output_file, nullptr, _IOLBF, 1);
  }

  // single rank only for now
  assert_always(globals::rank_global == 0);
  assert_always(globals::nprocs == 1);

  const time_t sys_time_start = time(nullptr);

  printout("Begining do_exspec.\n");

  /// Get input stuff
  printout("time before input %ld\n", time(nullptr));
  input(globals::rank_global);
  printout("time after input %ld\n", time(nullptr));
  // nprocs_exspec is the number of rank output files to process with expec
  // however, we might be running exspec with 1 or just a few ranks
  globals::nprocs = globals::nprocs_exspec;

  constexpr double maxpktmem_mb = 6000;
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
  struct packet *pkts = static_cast<struct packet *>(malloc(npkts_loaded * sizeof(struct packet)));

  init_spectrum_trace();  // needed for TRACE_EMISSION_ABSORPTION_REGION_ON

  struct spec *rpkt_spectra = alloc_spectra(globals::do_emission_res);

  struct spec *stokes_i = nullptr;
  struct spec *stokes_q = nullptr;
  struct spec *stokes_u = nullptr;

  if constexpr (POL_ON) {
    stokes_i = alloc_spectra(globals::do_emission_res);
    stokes_q = alloc_spectra(globals::do_emission_res);
    stokes_u = alloc_spectra(globals::do_emission_res);
  }

  struct spec *gamma_spectra = alloc_spectra(false);

  /// Initialise the grid. Call routine that sets up the initial positions
  /// and sizes of the grid cells.
  // grid_init();
  time_init();

  const int amax = ((grid::get_model_type() == grid::RHO_1D_READ)) ? 0 : MABINS;
  // a is the escape direction angle bin
  for (int a = -1; a < amax; a++) {
    /// Set up the light curve grid and initialise the bins to zero.
    double *rpkt_light_curve_lum = static_cast<double *>(calloc(globals::ntstep, sizeof(double)));
    double *rpkt_light_curve_lumcmf = static_cast<double *>(calloc(globals::ntstep, sizeof(double)));
    double *gamma_light_curve_lum = static_cast<double *>(calloc(globals::ntstep, sizeof(double)));
    double *gamma_light_curve_lumcmf = static_cast<double *>(calloc(globals::ntstep, sizeof(double)));
    /// Set up the spectrum grid and initialise the bins to zero.

    init_spectra(rpkt_spectra, NU_MIN_R, NU_MAX_R, globals::do_emission_res);

    if constexpr (POL_ON) {
      init_spectra(stokes_i, NU_MIN_R, NU_MAX_R, globals::do_emission_res);
      init_spectra(stokes_q, NU_MIN_R, NU_MAX_R, globals::do_emission_res);
      init_spectra(stokes_u, NU_MIN_R, NU_MAX_R, globals::do_emission_res);
    }

    const double nu_min_gamma = 0.05 * MEV / H;
    const double nu_max_gamma = 4. * MEV / H;
    init_spectra(gamma_spectra, nu_min_gamma, nu_max_gamma, false);

    for (int p = 0; p < globals::nprocs_exspec; p++) {
      struct packet *pkts_start = load_allrank_packets ? &pkts[p * globals::npkts] : pkts;

      if (a == -1 || !load_allrank_packets) {
        char pktfilename[128];
        snprintf(pktfilename, 128, "packets%.2d_%.4d.out", 0, p);
        printout("reading %s (file %d of %d)\n", pktfilename, p + 1, globals::nprocs_exspec);

        if (!access(pktfilename, F_OK)) {
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
            add_to_lc_res(&pkts_start[ii], a, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
            add_to_spec_res(&pkts_start[ii], a, rpkt_spectra, stokes_i, stokes_q, stokes_u);
          } else if (pkts_start[ii].escape_type == TYPE_GAMMA) {
            nesc_gamma++;
            if (a == -1) {
              add_to_lc_res(&pkts_start[ii], a, gamma_light_curve_lum, gamma_light_curve_lumcmf);
              add_to_spec_res(&pkts_start[ii], a, gamma_spectra, nullptr, nullptr, nullptr);
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
      write_light_curve("light_curve.out", -1, rpkt_light_curve_lum, rpkt_light_curve_lumcmf, globals::ntstep);
      write_light_curve("gamma_light_curve.out", -1, gamma_light_curve_lum, gamma_light_curve_lumcmf, globals::ntstep);

      write_spectrum("spec.out", "emission.out", "emissiontrue.out", "absorption.out", rpkt_spectra, globals::ntstep);

      if constexpr (POL_ON) {
        write_specpol("specpol.out", "emissionpol.out", "absorptionpol.out", stokes_i, stokes_q, stokes_u);
      }

      write_spectrum("gamma_spec.out", nullptr, nullptr, nullptr, gamma_spectra, globals::ntstep);

      printout("finished angle-averaged stuff\n");
    } else {
      // direction bin a
      // line-of-sight dependent spectra and light curves

      char lc_filename[128] = "";
      snprintf(lc_filename, 128, "light_curve_res_%.2d.out", a);

      char spec_filename[128] = "";
      snprintf(spec_filename, 128, "spec_res_%.2d.out", a);

      char emission_filename[128] = "";
      snprintf(emission_filename, 128, "emission_res_%.2d.out", a);

      char trueemission_filename[128] = "";
      snprintf(trueemission_filename, 128, "emissiontrue_res_%.2d.out", a);

      char absorption_filename[128] = "";
      snprintf(absorption_filename, 128, "absorption_res_%.2d.out", a);

      write_light_curve(lc_filename, a, rpkt_light_curve_lum, rpkt_light_curve_lumcmf, globals::ntstep);
      write_spectrum(spec_filename, emission_filename, trueemission_filename, absorption_filename, rpkt_spectra,
                     globals::ntstep);

      if constexpr (POL_ON) {
        char specpol_filename[128] = "";
        snprintf(specpol_filename, 128, "specpol_res_%.2d.out", a);

        char emissionpol_filename[128] = "";
        snprintf(emissionpol_filename, 128, "emissionpol_res_%.2d.out", a);

        char absorptionpol_filename[128] = "";
        snprintf(absorptionpol_filename, 128, "absorptionpol_res_%.2d.out", a);

        write_specpol(specpol_filename, emissionpol_filename, absorptionpol_filename, stokes_i, stokes_q, stokes_u);
      }

      printout("Did %d of %d angle bins.\n", a + 1, MABINS);
    }

    free(rpkt_light_curve_lum);
    free(rpkt_light_curve_lumcmf);

    free(gamma_light_curve_lum);
    free(gamma_light_curve_lumcmf);
  }

  free_spectra(rpkt_spectra);

  if (stokes_i != nullptr) {
    free_spectra(stokes_i);
  }

  if (stokes_q != nullptr) {
    free_spectra(stokes_q);
  }

  if (stokes_u != nullptr) {
    free_spectra(stokes_u);
  }

  free_spectra(gamma_spectra);

  free(pkts);
  decay::cleanup();

  printout("exspec finished at %ld (tstart + %ld seconds)\n", time(nullptr), time(nullptr) - sys_time_start);
  fclose(output_file);

#ifdef MPI_ON
  MPI_Finalize();
#endif

  return 0;
}
