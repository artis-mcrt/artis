#include "exspec.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#ifdef MPI_ON
#include <mpi.h>
#endif
#ifndef GPU_ON
#include <random>
#endif
#include <vector>

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

#ifndef GPU_ON
std::mt19937 stdrng{std::random_device{}()};
#endif

std::ofstream output_file;

namespace {

void do_angle_bin(const int a, Packet *pkts, bool load_allrank_packets, Spectra &rpkt_spectra, Spectra &stokes_i,
                  Spectra &stokes_q, Spectra &stokes_u, Spectra &gamma_spectra) {
  std::vector<double> rpkt_light_curve_lum(globals::ntimesteps, 0.);
  std::vector<double> rpkt_light_curve_lumcmf(globals::ntimesteps, 0.);
  std::vector<double> gamma_light_curve_lum(globals::ntimesteps, 0.);
  std::vector<double> gamma_light_curve_lumcmf(globals::ntimesteps, 0.);

  // Set up the spectrum grid and initialise the bins to zero.
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
    Packet *pkts_start = load_allrank_packets ? &pkts[p * globals::npkts] : pkts;

    if (a == -1 || !load_allrank_packets) {
      char pktfilename[MAXFILENAMELENGTH];
      snprintf(pktfilename, MAXFILENAMELENGTH, "packets%.2d_%.4d.out", 0, p);
      printout("reading %s (file %d of %d)\n", pktfilename, p + 1, globals::nprocs_exspec);

      if (std::filesystem::exists(pktfilename)) {
        read_packets(pktfilename, pkts_start);
      } else {
        printout("   WARNING %s does not exist - trying temp packets file at beginning of timestep %d...\n",
                 pktfilename, globals::timestep_initial);
        read_temp_packetsfile(globals::timestep_initial, p, pkts_start);
      }
    }

#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (p % globals::nprocs != globals::my_rank) {
      printout("skipping packets file %d %d\n", p + 1, globals::nprocs);
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
      printout("  %d of %d packets escaped (%d gamma-pkts and %d r-pkts)\n", nesc_tot, globals::npkts, nesc_gamma,
               nesc_rpkt);
    }
  }

  if (a == -1) {
    // angle-averaged spectra and light curves
    write_light_curve("light_curve.out", -1, rpkt_light_curve_lum, rpkt_light_curve_lumcmf, globals::ntimesteps);
    write_light_curve("gamma_light_curve.out", -1, gamma_light_curve_lum, gamma_light_curve_lumcmf,
                      globals::ntimesteps);

    write_spectrum("spec.out", "emission.out", "emissiontrue.out", "absorption.out", rpkt_spectra, globals::ntimesteps);

    if constexpr (POL_ON) {
      write_specpol("specpol.out", "emissionpol.out", "absorptionpol.out", &stokes_i, &stokes_q, &stokes_u);
    }

    write_spectrum("gamma_spec.out", "", "", "", gamma_spectra, globals::ntimesteps);

    printout("finished angle-averaged stuff\n");
  } else {
    // direction bin a
    // line-of-sight dependent spectra and light curves

    char lc_filename[MAXFILENAMELENGTH] = "";
    snprintf(lc_filename, sizeof(lc_filename), "light_curve_res_%.2d.out", a);

    char spec_filename[MAXFILENAMELENGTH] = "";
    snprintf(spec_filename, sizeof(spec_filename), "spec_res_%.2d.out", a);

    char emission_filename[MAXFILENAMELENGTH] = "";
    snprintf(emission_filename, sizeof(emission_filename), "emission_res_%.2d.out", a);

    char trueemission_filename[MAXFILENAMELENGTH] = "";
    snprintf(trueemission_filename, sizeof(trueemission_filename), "emissiontrue_res_%.2d.out", a);

    char absorption_filename[MAXFILENAMELENGTH] = "";
    snprintf(absorption_filename, sizeof(absorption_filename), "absorption_res_%.2d.out", a);

    write_light_curve(lc_filename, a, rpkt_light_curve_lum, rpkt_light_curve_lumcmf, globals::ntimesteps);
    write_spectrum(spec_filename, emission_filename, trueemission_filename, absorption_filename, rpkt_spectra,
                   globals::ntimesteps);

    if constexpr (POL_ON) {
      char specpol_filename[MAXFILENAMELENGTH] = "";
      snprintf(specpol_filename, sizeof(specpol_filename), "specpol_res_%.2d.out", a);

      char emissionpol_filename[MAXFILENAMELENGTH] = "";
      snprintf(emissionpol_filename, sizeof(emissionpol_filename), "emissionpol_res_%.2d.out", a);

      char absorptionpol_filename[MAXFILENAMELENGTH] = "";
      snprintf(absorptionpol_filename, sizeof(absorptionpol_filename), "absorptionpol_res_%.2d.out", a);

      write_specpol(specpol_filename, emissionpol_filename, absorptionpol_filename, &stokes_i, &stokes_q, &stokes_u);
    }

    printout("Did %d of %d angle bins.\n", a + 1, MABINS);
  }
}

}  // anonymous namespace

auto main(int argc, char *argv[]) -> int {  // NOLINT(misc-unused-parameters)
  const auto sys_time_start = std::time(nullptr);

#ifdef MPI_ON
  MPI_Init(&argc, &argv);
#endif

  globals::setup_mpi_vars();

  check_already_running();

  char filename[MAXFILENAMELENGTH];
  if (globals::my_rank == 0) {
    snprintf(filename, MAXFILENAMELENGTH, "exspec.txt");
    output_file = std::ofstream(filename);
    assert_always(output_file.is_open());
  }

  printout("git branch %s\n", GIT_BRANCH);

  printout("git version: %s\n", GIT_VERSION);

  printout("git status %s\n", GIT_STATUS);

  printout("exspec compiled at %s on %s\n", __TIME__, __DATE__);

#if defined TESTMODE && TESTMODE
  printout("TESTMODE is ON\n");
#endif

#ifdef MPI_ON
  printout("process id (pid): %d\n", getpid());
  printout("MPI enabled:\n");
  printout("  rank_global %d of [0..%d] in MPI_COMM_WORLD\n", globals::my_rank, globals::nprocs - 1);
  printout("  rank_in_node %d of [0..%d] in node %d of [0..%d]\n", globals::rank_in_node, globals::node_nprocs - 1,
           globals::node_id, globals::node_count - 1);
#else
  printout("MPI is disabled in this build\n");
#endif

  // single rank only for now
  assert_always(globals::my_rank == 0);
  assert_always(globals::nprocs == 1);

  printout("Beginning exspec.\n");

  // Get input stuff
  printout("time before input %ld\n", std::time(nullptr));
  input(globals::my_rank);
  printout("time after input %ld\n", std::time(nullptr));

  // nprocs_exspec is the number of rank output files to process with expec
  // however, we might be running exspec with 1 or just a few ranks

  auto *pkts = static_cast<Packet *>(malloc(globals::nprocs_exspec * globals::npkts * sizeof(Packet)));
  const bool load_allrank_packets = (pkts != nullptr);
  if (load_allrank_packets) {
    printout("mem_usage: loading %d packets from each %d processes simultaneously (total %d packets, %.1f MB memory)\n",
             globals::npkts, globals::nprocs_exspec, globals::nprocs_exspec * globals::npkts,
             globals::nprocs_exspec * globals::npkts * sizeof(Packet) / 1024. / 1024.);
  } else {
    printout("mem_usage: malloc failed to allocate memory for all packets\n");
    printout(
        "mem_usage: loading %d packets from each of %d processes sequentially (total %d packets, %.1f MB memory)\n",
        globals::npkts, globals::nprocs_exspec, globals::nprocs_exspec * globals::npkts,
        globals::nprocs_exspec * globals::npkts * sizeof(Packet) / 1024. / 1024.);
    pkts = static_cast<Packet *>(malloc(globals::npkts * sizeof(Packet)));
    assert_always(pkts != nullptr);
  }

  init_spectrum_trace();  // needed for TRACE_EMISSION_ABSORPTION_REGION_ON

  Spectra rpkt_spectra;

  Spectra stokes_i;
  Spectra stokes_q;
  Spectra stokes_u;

  Spectra gamma_spectra;

  time_init();

  const int amax = ((grid::get_model_type() == GridType::SPHERICAL1D)) ? 0 : MABINS;
  // a is the escape direction angle bin
  for (int a = -1; a < amax; a++) {
    do_angle_bin(a, pkts, load_allrank_packets, rpkt_spectra, stokes_i, stokes_q, stokes_u, gamma_spectra);
  }

  free(pkts);
  decay::cleanup();
  printout("exspec finished at %ld (tstart + %ld seconds)\n", std::time(nullptr), std::time(nullptr) - sys_time_start);

#ifdef MPI_ON
  MPI_Finalize();
#endif

  if (std::filesystem::exists("artis.pid")) {
    std::filesystem::remove("artis.pid");
  }

  return 0;
}
