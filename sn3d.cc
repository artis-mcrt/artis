/* 2007-10-30 -- MK
   Non-grey treatment of UVOIR opacity as opacity_case 4 added.
   Still not fully commented.
   Comments are marked by //  Deactivated code by // */
/* 2007-01-17 -- MK
   Several minor modifications (some marked in the code with //MK), these include
     - global printout() routine (located in sn3d.c)
     - opacity_cases 2 and 3 added (changes in grid_init.c and update_grid.c,
       original opacity stuff was moved there from input.c) */
/* This is a code copied from Lucy 2004 paper on t-dependent supernova
   explosions. */

#include "sn3d.h"

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <getopt.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#ifndef GPU_ON
#include <random>
#endif
#include <span>
#ifdef STDPAR_ON
#include <thread>
#endif

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "macroatom.h"
#ifdef MPI_ON
#include "rpkt.h"
#endif
#include "nltepop.h"
#include "nonthermal.h"
#include "packet.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "spectrum_lightcurve.h"
#include "stats.h"
#include "update_grid.h"
#include "update_packets.h"
#include "version.h"
#include "vpkt.h"

#ifndef GPU_ON
std::mt19937 stdrng{std::random_device{}()};
#endif

std::ofstream output_file;

namespace {
constexpr bool VERIFY_WRITTEN_PACKETS_FILES = false;

FILE *linestat_file{};
auto real_time_start = -1;
auto time_timestep_start = -1;  // this will be set after the first update of the grid and before packet prop
FILE *estimators_file{};

void initialise_linestat_file() {
  if (globals::simulation_continued_from_saved && !RECORD_LINESTAT) {
    // only write linestat.out on the first run, unless it contains statistics for each timestep
    return;
  }

  linestat_file = fopen_required("linestat.out", "w");

  for (int i = 0; i < globals::nlines; i++) {
    fprintf(linestat_file, "%g ", CLIGHT / globals::linelist[i].nu);
  }
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++) {
    fprintf(linestat_file, "%d ", get_atomicnumber(globals::linelist[i].elementindex));
  }
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++) {
    fprintf(linestat_file, "%d ", get_ionstage(globals::linelist[i].elementindex, globals::linelist[i].ionindex));
  }
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++) {
    fprintf(linestat_file, "%d ", globals::linelist[i].upperlevelindex + 1);
  }
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++) {
    fprintf(linestat_file, "%d ", globals::linelist[i].lowerlevelindex + 1);
  }
  fprintf(linestat_file, "\n");

  fflush(linestat_file);
}

void write_deposition_file() {
  const int my_rank = globals::my_rank;
  const int nts = globals::timestep;
  printout("Calculating deposition rates...\n");
  auto const time_write_deposition_file_start = std::time(nullptr);
  double mtot = 0.;
  const int nstart_nonempty = grid::get_nstart_nonempty(my_rank);
  const int ndo_nonempty = grid::get_ndo_nonempty(my_rank);

  // calculate analytical decay rates
  const double t_mid_nts = globals::timesteps[nts].mid;

  // power in [erg/s]
  globals::timesteps[nts].eps_positron_ana_power = 0.;
  globals::timesteps[nts].eps_electron_ana_power = 0.;
  globals::timesteps[nts].eps_alpha_ana_power = 0.;
  globals::timesteps[nts].qdot_betaminus = 0.;
  globals::timesteps[nts].qdot_alpha = 0.;
  globals::timesteps[nts].qdot_total = 0.;

  for (int nonemptymgi = nstart_nonempty; nonemptymgi < (nstart_nonempty + ndo_nonempty); nonemptymgi++) {
    const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    const double cellmass = grid::get_rho_tmin(mgi) * grid::get_modelcell_assocvolume_tmin(mgi);

    globals::timesteps[nts].eps_positron_ana_power +=
        cellmass * decay::get_particle_injection_rate(nonemptymgi, t_mid_nts, decay::DECAYTYPE_BETAPLUS);
    globals::timesteps[nts].eps_electron_ana_power +=
        cellmass * decay::get_particle_injection_rate(nonemptymgi, t_mid_nts, decay::DECAYTYPE_BETAMINUS);
    globals::timesteps[nts].eps_alpha_ana_power +=
        cellmass * decay::get_particle_injection_rate(nonemptymgi, t_mid_nts, decay::DECAYTYPE_ALPHA);

    mtot += cellmass;

    for (const auto decaytype : decay::all_decaytypes) {
      // Qdot here has been multiplied by mass, so it is in units of [erg/s]
      const double qdot_cell = decay::get_qdot_modelcell(nonemptymgi, t_mid_nts, decaytype) * cellmass;
      globals::timesteps[nts].qdot_total += qdot_cell;
      if (decaytype == decay::DECAYTYPE_BETAMINUS) {
        globals::timesteps[nts].qdot_betaminus += qdot_cell;
      } else if (decaytype == decay::DECAYTYPE_ALPHA) {
        globals::timesteps[nts].qdot_alpha += qdot_cell;
      }
    }
  }

#ifdef MPI_ON
  // in MPI mode, each process only calculated the contribution of a subset of cells
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].eps_positron_ana_power, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].eps_electron_ana_power, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].eps_alpha_ana_power, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].qdot_betaminus, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].qdot_alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].qdot_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

#ifdef MPI_ON
  MPI_Allreduce(MPI_IN_PLACE, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (my_rank == 0) {
    FILE *dep_file = fopen_required("deposition.out.tmp", "w");
    fprintf(dep_file,
            "#ts tmid_days tmid_s total_dep_Lsun gammadep_discrete_Lsun gammadep_Lsun positrondep_Lsun "
            "eps_positron_ana_Lsun elecdep_Lsun eps_elec_Lsun eps_elec_ana_Lsun alphadep_Lsun eps_alpha_Lsun "
            "eps_alpha_ana_Lsun eps_gamma_Lsun Qdot_betaminus_ana_erg/s/g Qdotalpha_ana_erg/s/g eps_erg/s/g "
            "Qdot_ana_erg/s/g positrondep_discrete_Lsun elecdep_discrete_Lsun alphadep_discrete_Lsun\n");

    for (int i = 0; i <= nts; i++) {
      const double t_mid = globals::timesteps[i].mid;
      const double t_width = globals::timesteps[i].width;
      const double total_dep = (globals::timesteps[i].gamma_dep + globals::timesteps[i].positron_dep +
                                globals::timesteps[i].electron_dep + globals::timesteps[i].alpha_dep);

      const double epsilon_tot = (globals::timesteps[i].gamma_emission + globals::timesteps[i].positron_emission +
                                  globals::timesteps[i].electron_emission + globals::timesteps[i].alpha_emission) /
                                 mtot / t_width;

      fprintf(dep_file, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", i, t_mid / DAY, t_mid,
              total_dep / t_width / LSUN, globals::timesteps[i].gamma_dep_discrete / t_width / LSUN,
              globals::timesteps[i].gamma_dep / t_width / LSUN, globals::timesteps[i].positron_dep / t_width / LSUN,
              globals::timesteps[i].eps_positron_ana_power / LSUN, globals::timesteps[i].electron_dep / t_width / LSUN,
              globals::timesteps[i].electron_emission / t_width / LSUN,
              globals::timesteps[i].eps_electron_ana_power / LSUN, globals::timesteps[i].alpha_dep / t_width / LSUN,
              globals::timesteps[i].alpha_emission / t_width / LSUN, globals::timesteps[i].eps_alpha_ana_power / LSUN,
              globals::timesteps[i].gamma_emission / t_width / LSUN, globals::timesteps[i].qdot_betaminus / mtot,
              globals::timesteps[i].qdot_alpha / mtot, epsilon_tot, globals::timesteps[i].qdot_total / mtot,
              globals::timesteps[i].positron_dep_discrete / t_width / LSUN,
              globals::timesteps[i].electron_dep_discrete / t_width / LSUN,
              globals::timesteps[i].alpha_dep_discrete / t_width / LSUN);
    }
    fclose(dep_file);

    std::remove("deposition.out");
    std::rename("deposition.out.tmp", "deposition.out");
  }

  printout("calculating and writing deposition.out took %ld seconds\n",
           std::time(nullptr) - time_write_deposition_file_start);
}

#ifdef MPI_ON
void mpi_communicate_grid_properties() {
  const auto nincludedions = get_includedions();
  const auto nelements = get_nelements();
  for (int root = 0; root < globals::nprocs; root++) {
    MPI_Barrier(MPI_COMM_WORLD);

    int root_node_id = globals::node_id;
    MPI_Bcast(&root_node_id, 1, MPI_INT, root, MPI_COMM_WORLD);

    const int root_nstart_nonempty = grid::get_nstart_nonempty(root);
    const int root_ndo_nonempty = grid::get_ndo_nonempty(root);

    for (int nonemptymgi = root_nstart_nonempty; nonemptymgi < (root_nstart_nonempty + root_ndo_nonempty);
         nonemptymgi++) {
      assert_always(root_ndo_nonempty > 0);

      radfield::do_MPI_Bcast(nonemptymgi, root, root_node_id);

      nonthermal::nt_MPI_Bcast(nonemptymgi, root, root_node_id);

      if (globals::total_nlte_levels > 0 && globals::rank_in_node == 0) {
        MPI_Bcast(&grid::nltepops_allcells[nonemptymgi * globals::total_nlte_levels], globals::total_nlte_levels,
                  MPI_DOUBLE, root_node_id, globals::mpi_comm_internode);
      }

      if (USE_LUT_PHOTOION && globals::nbfcontinua_ground > 0) {
        assert_always(globals::corrphotoionrenorm != nullptr);
        if (globals::rank_in_node == 0) {
          MPI_Bcast(&globals::corrphotoionrenorm[nonemptymgi * globals::nbfcontinua_ground],
                    globals::nbfcontinua_ground, MPI_DOUBLE, root_node_id, globals::mpi_comm_internode);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&globals::gammaestimator[nonemptymgi * globals::nbfcontinua_ground], globals::nbfcontinua_ground,
                  MPI_DOUBLE, root, MPI_COMM_WORLD);
      }
      if (globals::rank_in_node == 0) {
        MPI_Bcast(&grid::elem_meanweight_allcells[nonemptymgi * nelements], nelements, MPI_FLOAT, root_node_id,
                  globals::mpi_comm_internode);
        MPI_Bcast(&grid::elem_massfracs_allcells[nonemptymgi * nelements], nelements, MPI_FLOAT, root_node_id,
                  globals::mpi_comm_internode);
        MPI_Bcast(&grid::ion_groundlevelpops_allcells[nonemptymgi * nincludedions], nincludedions, MPI_FLOAT,
                  root_node_id, globals::mpi_comm_internode);
        MPI_Bcast(&grid::ion_partfuncts_allcells[nonemptymgi * nincludedions], nincludedions, MPI_FLOAT, root_node_id,
                  globals::mpi_comm_internode);
        MPI_Bcast(&grid::ion_cooling_contribs_allcells[nonemptymgi * nincludedions], nincludedions, MPI_DOUBLE,
                  root_node_id, globals::mpi_comm_internode);
        MPI_Bcast(&grid::modelgrid[grid::get_mgi_of_nonemptymgi(nonemptymgi)], sizeof(grid::ModelGridCell), MPI_BYTE,
                  root_node_id, globals::mpi_comm_internode);
      }

      MPI_Bcast_binned_opacities(nonemptymgi, root_node_id);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

void mpi_reduce_estimators(const int nts) {
  const int nonempty_npts_model = grid::get_nonempty_npts_model();
  radfield::reduce_estimators();
  MPI_Barrier(MPI_COMM_WORLD);
  assert_always(!globals::ffheatingestimator.empty());
  MPI_Allreduce(MPI_IN_PLACE, globals::ffheatingestimator.data(), globals::ffheatingestimator.size(), MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  assert_always(!globals::colheatingestimator.empty());
  MPI_Allreduce(MPI_IN_PLACE, globals::colheatingestimator.data(), globals::colheatingestimator.size(), MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (globals::nbfcontinua_ground > 0) {
    if constexpr (USE_LUT_PHOTOION) {
      MPI_Barrier(MPI_COMM_WORLD);
      assert_always(!globals::gammaestimator.empty());
      MPI_Allreduce(MPI_IN_PLACE, globals::gammaestimator.data(), globals::gammaestimator.size(), MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
    }

    if constexpr (USE_LUT_BFHEATING) {
      MPI_Barrier(MPI_COMM_WORLD);
      assert_always(!globals::bfheatingestimator.empty());
      MPI_Allreduce(MPI_IN_PLACE, globals::bfheatingestimator.data(), nonempty_npts_model * globals::nbfcontinua_ground,
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
  }

  if constexpr (RECORD_LINESTAT) {
    MPI_Barrier(MPI_COMM_WORLD);
    assert_always(!globals::ecounter.empty());
    MPI_Allreduce(MPI_IN_PLACE, globals::ecounter.data(), globals::ecounter.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    assert_always(!globals::acounter.empty());
    MPI_Allreduce(MPI_IN_PLACE, globals::acounter.data(), globals::acounter.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }

  assert_always(static_cast<int>(globals::dep_estimator_gamma.size()) == nonempty_npts_model);
  MPI_Allreduce(MPI_IN_PLACE, globals::dep_estimator_gamma.data(), nonempty_npts_model, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  assert_always(static_cast<int>(globals::dep_estimator_positron.size()) == nonempty_npts_model);
  MPI_Allreduce(MPI_IN_PLACE, globals::dep_estimator_positron.data(), nonempty_npts_model, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  assert_always(static_cast<int>(globals::dep_estimator_electron.size()) == nonempty_npts_model);
  MPI_Allreduce(MPI_IN_PLACE, globals::dep_estimator_electron.data(), nonempty_npts_model, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  assert_always(static_cast<int>(globals::dep_estimator_alpha.size()) == nonempty_npts_model);
  MPI_Allreduce(MPI_IN_PLACE, globals::dep_estimator_alpha.data(), nonempty_npts_model, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].cmf_lum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].cmf_lum /= globals::nprocs;

  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].gamma_dep_discrete, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].gamma_dep_discrete /= globals::nprocs;

  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].positron_dep_discrete, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].positron_dep_discrete /= globals::nprocs;

  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].positron_emission, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].positron_emission /= globals::nprocs;

  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].electron_dep_discrete, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].electron_dep_discrete /= globals::nprocs;

  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].electron_emission, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].electron_emission /= globals::nprocs;

  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].alpha_dep_discrete, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].alpha_dep_discrete /= globals::nprocs;

  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].alpha_emission, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].alpha_emission /= globals::nprocs;

  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].gamma_emission, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].gamma_emission /= globals::nprocs;

  if constexpr (TRACK_ION_STATS) {
    stats::reduce_estimators();
  }

  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

void write_temp_packetsfile(const int timestep, const int my_rank, const Packet *pkt) {
  // write packets binary file (and retry if the write fails)
  char filename[MAXFILENAMELENGTH];
  snprintf(filename, MAXFILENAMELENGTH, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  bool write_success = false;
  while (!write_success) {
    printout("Writing %s...", filename);
    FILE *packets_file = fopen(filename, "wb");
    if (packets_file == nullptr) {
      printout("ERROR: Could not open file '%s' for mode 'wb'.\n", filename);
      write_success = false;
    } else {
      write_success =
          (std::fwrite(pkt, sizeof(Packet), globals::npkts, packets_file) == static_cast<size_t>(globals::npkts));
      if (!write_success) {
        printout("fwrite() FAILED! will retry...\n");
      }

      fclose(packets_file);
    }

    if (write_success) {
      printout("done\n");
    }
  }
}

void remove_temp_packetsfile(const int timestep, const int my_rank) {
  char filename[MAXFILENAMELENGTH];
  snprintf(filename, MAXFILENAMELENGTH, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  if (std::filesystem::exists(filename)) {
    std::remove(filename);
    printout("Deleted %s\n", filename);
  }
}

void remove_grid_restart_data(const int timestep) {
  char prevfilename[MAXFILENAMELENGTH];
  snprintf(prevfilename, MAXFILENAMELENGTH, "gridsave_ts%d.tmp", timestep);

  if (std::filesystem::exists(prevfilename)) {
    std::remove(prevfilename);
    printout("Deleted %s\n", prevfilename);
  }
}

auto walltime_sufficient_to_continue(const int nts, const int nts_prev, const int walltimelimitseconds) -> bool {
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // time is measured from just before packet propagation from one timestep to the next
  const int estimated_time_per_timestep = std::time(nullptr) - time_timestep_start;
  printout("TIME: time between timesteps is %d seconds (measured packet prop of ts %d and update grid of ts %d)\n",
           estimated_time_per_timestep, nts_prev, nts);

  bool do_this_full_loop = true;
  if (walltimelimitseconds > 0) {
    const int wallclock_used_seconds = std::time(nullptr) - real_time_start;
    const int wallclock_remaining_seconds = walltimelimitseconds - wallclock_used_seconds;
    printout("TIMED_RESTARTS: Used %d of %d seconds of wall time.\n", wallclock_used_seconds, walltimelimitseconds);

    // This flag being false will make it update_grid, and then exit
    do_this_full_loop = (wallclock_remaining_seconds >= (1.5 * estimated_time_per_timestep));

#ifdef MPI_ON
    // communicate whatever decision the rank 0 process decided, just in case they differ
    MPI_Bcast(&do_this_full_loop, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif
    if (do_this_full_loop) {
      printout("TIMED_RESTARTS: Going to continue since remaining time %d s >= 1.5 * time_per_timestep\n",
               wallclock_remaining_seconds);
    } else {
      printout("TIMED_RESTARTS: Going to terminate since remaining time %d s < 1.5 * time_per_timestep\n",
               wallclock_remaining_seconds);
    }
  }
  return do_this_full_loop;
}

void save_grid_and_packets(const int nts, const Packet *packets) {
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  const auto my_rank = globals::my_rank;

  bool write_successful = false;
  while (!write_successful) {
    const auto time_write_packets_file_start = std::time(nullptr);
    printout("time before write temporary packets file %ld\n", time_write_packets_file_start);

    // save packet state at start of current timestep (before propagation)
    write_temp_packetsfile(nts, globals::my_rank, packets);

    vpkt_write_timestep(nts, globals::my_rank, false);

    const auto time_write_packets_finished_thisrank = std::time(nullptr);

#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const auto timenow = std::time(nullptr);

    printout("time after write temporary packets file %ld (took %lds, waited %lds, total %lds)\n", timenow,
             time_write_packets_finished_thisrank - time_write_packets_file_start,
             timenow - time_write_packets_finished_thisrank, timenow - time_write_packets_file_start);

    if constexpr (VERIFY_WRITTEN_PACKETS_FILES) {
#ifdef MPI_ON
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      const auto time_readback_packets_start = std::time(nullptr);

      printout("reading back temporary packets file to check validity...\n");

      // read packets file back to check that the disk write didn't fail
      write_successful = verify_temp_packetsfile(nts, my_rank, packets);

#ifdef MPI_ON
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      printout("Verifying packets files for all ranks took %ld seconds.\n",
               std::time(nullptr) - time_readback_packets_start);
    } else {
      write_successful = true;
    }
  }

  if (my_rank == 0) {
    grid::write_grid_restart_data(nts);
    update_parameterfile(nts);
  }

  if (!KEEP_ALL_RESTART_FILES) {
    // ensure new packets files have been written by all processes before we remove the old set
#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (my_rank == 0) {
      remove_grid_restart_data(nts - 1);
    }

    // delete temp packets files from previous timestep now that all restart data for the new timestep is available
    remove_temp_packetsfile(nts - 1, my_rank);
    vpkt_remove_temp_file(nts - 1, my_rank);
  }
}

void zero_estimators() {
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  radfield::zero_estimators();
  if constexpr (TRACK_ION_STATS) {
    for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
      stats::reset_ion_stats(modelgridindex);
    }
  }

  std::ranges::fill(globals::ffheatingestimator, 0.);
  std::ranges::fill(globals::colheatingestimator, 0.);
  std::ranges::fill(globals::dep_estimator_gamma, 0.);
  std::ranges::fill(globals::dep_estimator_positron, 0.);
  std::ranges::fill(globals::dep_estimator_electron, 0.);
  std::ranges::fill(globals::dep_estimator_alpha, 0.);

  if constexpr (USE_LUT_PHOTOION) {
    if (globals::nbfcontinua_ground > 0) {
      std::ranges::fill(globals::gammaestimator, 0.);
    }
  }

  if constexpr (USE_LUT_BFHEATING) {
    if (globals::nbfcontinua_ground > 0) {
      std::ranges::fill(globals::bfheatingestimator, 0.);
    }
  }

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void normalise_deposition_estimators(int nts) {
  const double dt = globals::timesteps[nts].width;
  const auto nprocs = globals::nprocs;

  globals::timesteps[nts].gamma_dep = 0.;
  globals::timesteps[nts].positron_dep = 0.;
  globals::timesteps[nts].electron_dep = 0.;
  globals::timesteps[nts].alpha_dep = 0.;

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);

    const double dV =
        grid::get_modelcell_assocvolume_tmin(mgi) * std::pow(globals::timesteps[nts].mid / globals::tmin, 3);

    // contribute the energy deposited (in erg) by each process in this cell to the timestep total
    globals::timesteps[nts].gamma_dep += globals::dep_estimator_gamma[nonemptymgi] / nprocs;
    globals::timesteps[nts].positron_dep += globals::dep_estimator_positron[nonemptymgi] / nprocs;
    globals::timesteps[nts].electron_dep += globals::dep_estimator_electron[nonemptymgi] / nprocs;
    globals::timesteps[nts].alpha_dep += globals::dep_estimator_alpha[nonemptymgi] / nprocs;

    // normalise the estimators to units of erg/s/cm^3
    const double estimator_normfactor = 1 / dV / dt / nprocs;

    globals::dep_estimator_gamma[nonemptymgi] *= estimator_normfactor;
    globals::dep_estimator_positron[nonemptymgi] *= estimator_normfactor;
    globals::dep_estimator_electron[nonemptymgi] *= estimator_normfactor;
    globals::dep_estimator_alpha[nonemptymgi] *= estimator_normfactor;

    assert_testmodeonly(globals::dep_estimator_gamma[nonemptymgi] >= 0.);
    assert_testmodeonly(std::isfinite(globals::dep_estimator_gamma[nonemptymgi]));
  }
}

auto do_timestep(const int nts, const int titer, Packet *packets, const int walltimelimitseconds) -> bool {
  const auto my_rank = globals::my_rank;
  bool do_this_full_loop = true;
  const int nts_prev = (titer != 0 || nts == 0) ? nts : nts - 1;
  if ((titer > 0) || (globals::simulation_continued_from_saved && (nts == globals::timestep_initial))) {
    // Read the packets file to reset before each additional iteration on the timestep
    read_temp_packetsfile(nts, globals::my_rank, packets);
  }

  // Some counters on pkt-actions need to be reset to do statistics
  stats::pkt_action_counters_reset();

  if (nts == 0) {
    radfield::initialise_prev_titer_photoionestimators();
  }

  if constexpr (RECORD_LINESTAT) {
    // The same for absorption/emission of r-pkts in lines
    for (int i = 0; i < globals::nlines; i++) {
      globals::acounter[i] = 0;
      globals::ecounter[i] = 0;
    }
  }

  // Update the matter quantities in the grid for the new timestep.

  update_grid(estimators_file, nts, nts_prev, titer, real_time_start);

  const auto sys_time_start_communicate_grid = std::time(nullptr);

  // Each process has now updated its own set of cells. The results now need to be communicated between processes.
#ifdef MPI_ON
  mpi_communicate_grid_properties();
#endif

  printout("timestep %d: time after grid properties have been communicated %ld (took %ld seconds)\n", nts,
           std::time(nullptr), std::time(nullptr) - sys_time_start_communicate_grid);

  // If this is not the 0th time step of the current job step,
  // write out a snapshot of the grid properties for further restarts
  // and update input.txt accordingly
  if (((nts - globals::timestep_initial) != 0)) {
    save_grid_and_packets(nts, packets);
    do_this_full_loop = walltime_sufficient_to_continue(nts, nts_prev, walltimelimitseconds);
  }
  time_timestep_start = std::time(nullptr);

  // set all the estimators to zero before moving packets. This is now done
  // after update_grid so that, if requires, the gamma-ray heating estimator is known there
  // and also the photoion and stimrecomb estimators
  zero_estimators();

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if ((nts < globals::timestep_finish) && do_this_full_loop) {
    // Now process the packets.

    update_packets(nts, std::span{packets, static_cast<size_t>(globals::npkts)});

#ifdef MPI_ON
    // All the processes have their own versions of the estimators for this time step now.
    // Since these are going to be needed in the next time step, we will gather all the
    // estimators together now, sum them, and distribute the results

    const auto time_communicate_estimators_start = std::time(nullptr);
    mpi_reduce_estimators(nts);
#endif

#ifdef MPI_ON
    printout("timestep %d: time after estimators have been communicated %ld (took %ld seconds)\n", nts,
             std::time(nullptr), std::time(nullptr) - time_communicate_estimators_start);
#endif

    // The estimators have been summed across all processes and distributed.
    // They will now be normalised independently on all processes.

    normalise_deposition_estimators(nts);

    write_deposition_file();

    write_partial_lightcurve_spectra(my_rank, nts, packets);

    printout("During timestep %d on MPI process %d, %d pellets decayed and %d packets escaped. (t=%gd)\n", nts, my_rank,
             globals::timesteps[nts].pellet_decays, globals::nesc, globals::timesteps[nts].mid / DAY);

    if (VPKT_ON) {
      printout("During timestep %d on MPI process %d, %d virtual packets were generated and %d escaped.\n", nts,
               my_rank, nvpkt, nvpkt_esc1 + nvpkt_esc2 + nvpkt_esc3);
      printout(
          "%d virtual packets came from an electron scattering event, %d from a kpkt deactivation and %d from a "
          "macroatom deactivation.\n",
          nvpkt_esc1, nvpkt_esc2, nvpkt_esc3);

      nvpkt = 0;
      nvpkt_esc1 = 0;
      nvpkt_esc2 = 0;
      nvpkt_esc3 = 0;
    }

    if constexpr (RECORD_LINESTAT) {
      if (my_rank == 0) {
        // Print net absorption/emission in lines to the linestat_file
        // Currently linestat information is only properly implemented for MPI only runs
        // For hybrid runs only data from thread 0 is recorded
        for (int i = 0; i < globals::nlines; i++) {
          fprintf(linestat_file, "%d ", globals::ecounter[i]);
        }
        fprintf(linestat_file, "\n");
        for (int i = 0; i < globals::nlines; i++) {
          fprintf(linestat_file, "%d ", globals::acounter[i]);
        }
        fprintf(linestat_file, "\n");
        fflush(linestat_file);
      }
    }

    if (nts == globals::timestep_finish - 1) {
      char filename[MAXFILENAMELENGTH];
      snprintf(filename, MAXFILENAMELENGTH, "packets%.2d_%.4d.out", 0, my_rank);
      // snprintf(filename, MAXFILENAMELENGTH, "packets%.2d_%.4d.out", middle_iteration, my_rank);
      write_packets(filename, packets);

      vpkt_write_timestep(nts, my_rank, true);

      printout("time after write final packets file %ld\n", std::time(nullptr));

      // final packets*.out have been written, so remove the temporary packets files
      // commented out because you might still want to resume the simulation
      // snprintf(filename, MAXFILENAMELENGTH, "packets%d_%d_odd.tmp", 0, my_rank);
      // std::remove(filename);
      // snprintf(filename, MAXFILENAMELENGTH, "packets%d_%d_even.tmp", 0, my_rank);
      // std::remove(filename);
    }
  }
  return !do_this_full_loop;
}

}  // anonymous namespace

auto main(int argc, char *argv[]) -> int {
  real_time_start = std::time(nullptr);
  char filename[MAXFILENAMELENGTH];

  // if DETAILED_BF_ESTIMATORS_ON is true, USE_LUT_PHOTOION must be false
  assert_always(!DETAILED_BF_ESTIMATORS_ON || !USE_LUT_PHOTOION);

  if constexpr (VPKT_ON) {
    nvpkt = 0;
    nvpkt_esc1 = 0;
    nvpkt_esc2 = 0;
    nvpkt_esc3 = 0;
  }

#ifdef MPI_ON
  MPI_Init(&argc, &argv);
#endif

  globals::setup_mpi_vars();

  check_already_running();

  const int my_rank = globals::my_rank;

#if defined(_OPENMP) && !defined(GPU_ON)
  // Explicitly turn off dynamic threads because we use the threadprivate directive!!!
  omp_set_dynamic(0);
#pragma omp parallel private(filename)
#endif
  {
    // initialise the thread and rank specific output file
    snprintf(filename, MAXFILENAMELENGTH, "output_%d-%d.txt", my_rank, get_thread_num());
    output_file = std::ofstream(filename);
    assert_always(output_file.is_open());

#ifdef _OPENMP
    printout("OpenMP parallelisation is active with %d threads (max %d)\n", omp_get_num_threads(), get_max_threads());
#else
    printout("OpenMP parallelisation is not enabled in this build (this is normal)\n");
#endif
  }

#ifdef STDPAR_ON
  printout("C++ standard parallelism (stdpar) is enabled with %d hardware threads\n",
           std::thread::hardware_concurrency());
#endif

#ifdef GPU_ON
  printout("GPU_ON is enabled\n");
#endif

  printout("time at start %d\n", real_time_start);

  printout("integration method is %s\n", USE_SIMPSON_INTEGRATOR ? "Simpson rule" : "GSL qag");

#ifdef WALLTIMELIMITSECONDS
  int walltimelimitseconds = WALLTIMELIMITSECONDS;
#else
  int walltimelimitseconds = -1;
#endif

  int opt = 0;
  while ((opt = getopt(argc, argv, "w:")) != -1) {  // NOLINT(concurrency-mt-unsafe)
    if (opt == 'w') {
      printout("Command line argument specifies wall time hours '%s', setting ", optarg);
      const float walltimehours = strtof(optarg, nullptr);
      walltimelimitseconds = static_cast<int>(walltimehours * 3600);
      printout("walltimelimitseconds = %d\n", walltimelimitseconds);
    } else {
      fprintf(stderr, "Usage: %s [-w WALLTIMELIMITHOURS]\n", argv[0]);
      std::abort();
    }
  }

  auto *const packets = static_cast<Packet *>(malloc(globals::npkts * sizeof(Packet)));

  assert_always(packets != nullptr);

  printout("git branch %s\n", GIT_BRANCH);

  printout("git version: %s\n", GIT_VERSION);

  printout("git status %s\n", GIT_STATUS);

  printout("sn3d compiled at %s on %s\n", __TIME__, __DATE__);

#if defined TESTMODE && TESTMODE
  printout("TESTMODE is ON\n");
#endif

#ifdef MPI_ON
  printout("process id (pid): %d\n", getpid());
  printout("MPI enabled:\n");
  printout("  rank %d of [0..%d] in MPI_COMM_WORLD\n", globals::my_rank, globals::nprocs - 1);
  printout("  rank %d of [0..%d] in node %d of [0..%d]\n", globals::rank_in_node, globals::node_nprocs - 1,
           globals::node_id, globals::node_count - 1);
#ifdef MAX_NODE_SIZE
  printout(
      "WARNING: Compiled with MAX_NODE_SIZE %d, which may mean mean that there are more nodes reported than physically "
      "present\n",
      MAX_NODE_SIZE);
#endif
#else
  printout("MPI is disabled in this build\n");
#endif

  input(my_rank);
  if (globals::simulation_continued_from_saved) {
    assert_always(globals::nprocs_exspec == globals::nprocs);
  } else {
    globals::nprocs_exspec = globals::nprocs;
  }

  if (my_rank == 0) {
    initialise_linestat_file();
  }

  printout("time after input %ld\n", std::time(nullptr));
  printout("timesteps %d\n", globals::ntimesteps);

  // Precalculate the rate coefficients for spontaneous and stimulated recombination
  // and for photoionisation. With the nebular approximation they only depend on T_e
  // T_R and W. W is easily factored out. For stimulated recombination we must assume
  // T_e = T_R for this precalculation.

  ratecoefficients_init();

#ifdef MPI_ON
  printout("barrier after tabulation of rate coefficients: time before barrier %ld, ", std::time(nullptr));
  MPI_Barrier(MPI_COMM_WORLD);
  printout("time after barrier %ld\n", std::time(nullptr));
#endif

  stats::init();

  // Record the chosen syn_dir
  FILE *syn_file = fopen_required("syn_dir.txt", "w");
  fprintf(syn_file, "%g %g %g", globals::syn_dir[0], globals::syn_dir[1], globals::syn_dir[2]);
  fclose(syn_file);
  printout("time write syn_dir.txt file %ld\n", std::time(nullptr));

  bool terminate_early = false;

  time_init();

  if (my_rank == 0) {
    write_timestep_file();
  }

  printout("time grid_init %ld\n", std::time(nullptr));
  grid::grid_init(my_rank);

  printout("Simulation propagates %g packets per process (total %g with nprocs %d)\n", 1. * globals::npkts,
           1. * globals::npkts * globals::nprocs, globals::nprocs);

  printout("[info] mem_usage: packets occupy %.3f MB\n", MPKTS * sizeof(Packet) / 1024. / 1024.);

  if (!globals::simulation_continued_from_saved) {
    std::remove("deposition.out");
    packet_init(packets);
    zero_estimators();
  }

  // For the parallelisation of update_grid, the process needs to be told which cells belong to it.
  // The next loop is over all grid cells. For parallelisation, we want to split this loop between
  // processes. This is done by assigning each MPI process nblock cells. The residual n_leftover
  // cells are sent to processes 0 ... process n_leftover -1.
  const int nstart = grid::get_nstart(my_rank);
  const int ndo = grid::get_ndo(my_rank);
  const int ndo_nonempty = grid::get_ndo_nonempty(my_rank);
  printout("process rank %d (global max rank %d) assigned %d modelgrid cells (%d nonempty)", my_rank,
           globals::nprocs - 1, ndo, ndo_nonempty);
  if (ndo > 0) {
    printout(": cells [%d..%d] (model has max mgi %d)\n", nstart, nstart + ndo - 1, grid::get_npts_model() - 1);
  } else {
    printout("\n");
  }

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  globals::timestep = globals::timestep_initial;

  macroatom_open_file(my_rank);
  if (ndo > 0) {
    assert_always(estimators_file == nullptr);
    snprintf(filename, MAXFILENAMELENGTH, "estimators_%.4d.out", my_rank);
    estimators_file = fopen_required(filename, "w");

    if (globals::total_nlte_levels > 0 && ndo_nonempty > 0) {
      nltepop_open_file(my_rank);
    }
  }

  // initialise or read in virtual packet spectra
  vpkt_init(globals::timestep, my_rank, globals::simulation_continued_from_saved);

  while (globals::timestep < globals::timestep_finish && !terminate_early) {
#ifdef MPI_ON
    //        const auto time_before_barrier = std::time(nullptr);
    MPI_Barrier(MPI_COMM_WORLD);
    //        const auto time_after_barrier = std::time(nullptr);
    //        printout("timestep %d: time before barrier %d, time after barrier %d\n", nts, time_before_barrier,
    //        time_after_barrier);
#endif

    // titer example: Do 3 iterations on timestep 0-6
    // globals::n_titer = (nts < 6) ? 3: 1;

    globals::n_titer = 1;
    globals::lte_iteration = (globals::timestep < globals::num_lte_timesteps);
    assert_always(globals::num_lte_timesteps > 0);  // The first time step must solve the ionisation balance in LTE

    for (int titer = 0; titer < globals::n_titer; titer++) {
      terminate_early = do_timestep(globals::timestep, titer, packets, walltimelimitseconds);
#ifdef DO_TITER
      // No iterations over the zeroth timestep, set titer > n_titer
      if (globals::timestep == 0) titer = globals::n_titer + 1;
#endif
    }

    globals::timestep++;
  }

  // The main calculation is now over. The packets now have all stored the time, place and direction
  // at which they left the grid. Also their rest frame energies and frequencies.
  // Spectra and light curves are now extracted using exspec which is another make target of this
  // code.

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (linestat_file != nullptr) {
    fclose(linestat_file);
  }

  if ((globals::ntimesteps != globals::timestep_finish) || (terminate_early)) {
    printout("RESTART_NEEDED to continue model\n");
  } else {
    printout("No need for restart\n");
  }

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  const auto real_time_end = std::time(nullptr);
  printout("sn3d finished at %ld (this job wallclock hours %.2f * %d processes * %d threads = %.1f core hours)\n",
           real_time_end, (real_time_end - real_time_start) / 3600., globals::nprocs, get_max_threads(),
           (real_time_end - real_time_start) / 3600. * globals::nprocs * get_max_threads());

  if (estimators_file != nullptr) {
    fclose(estimators_file);
  }

  macroatom_close_file();
  nltepop_close_file();

  radfield::close_file();
  nonthermal::close_file();

  free(packets);

  decay::cleanup();

#ifdef MPI_ON
  MPI_Finalize();
#endif

  const std::filesystem::path pid_file_path("artis.pid");
  if (std::filesystem::exists(pid_file_path)) {
    std::filesystem::remove(pid_file_path);
  }

  return 0;
}
