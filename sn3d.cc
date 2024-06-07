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

#include "sn3d.h"

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <getopt.h>
#include <sys/unistd.h>
#include <unistd.h>

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <random>
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

std::mt19937 stdrng{std::random_device{}()};

std::ofstream output_file;

namespace {

FILE *linestat_file{};
auto real_time_start = -1;
auto time_timestep_start = -1;  // this will be set after the first update of the grid and before packet prop
FILE *estimators_file{};

#ifdef MPI_ON
size_t mpi_grid_buffer_size = 0;
char *mpi_grid_buffer{};
#endif

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

void write_deposition_file(const int nts, const int my_rank, const int nstart, const int ndo) {
  printout("Calculating deposition rates...\n");
  auto const time_write_deposition_file_start = std::time(nullptr);
  double mtot = 0.;

  // calculate analytical decay rates
  // for (int i = 0; i <= nts; i++)
  {
    const int i = nts;
    const double t_mid = globals::timesteps[i].mid;

    // power in [erg/s]
    globals::timesteps[i].eps_positron_ana_power = 0.;
    globals::timesteps[i].eps_electron_ana_power = 0.;
    globals::timesteps[i].eps_alpha_ana_power = 0.;
    globals::timesteps[i].qdot_betaminus = 0.;
    globals::timesteps[i].qdot_alpha = 0.;
    globals::timesteps[i].qdot_total = 0.;

    for (int mgi = nstart; mgi < (nstart + ndo); mgi++)
    // for (int mgi = 0; mgi < grid::get_npts_model(); mgi++)
    {
      if (grid::get_numassociatedcells(mgi) > 0) {
        const double cellmass = grid::get_rho_tmin(mgi) * grid::get_modelcell_assocvolume_tmin(mgi);

        globals::timesteps[i].eps_positron_ana_power +=
            cellmass * decay::get_particle_injection_rate(mgi, t_mid, decay::DECAYTYPE_BETAPLUS);
        globals::timesteps[i].eps_electron_ana_power +=
            cellmass * decay::get_particle_injection_rate(mgi, t_mid, decay::DECAYTYPE_BETAMINUS);
        globals::timesteps[i].eps_alpha_ana_power +=
            cellmass * decay::get_particle_injection_rate(mgi, t_mid, decay::DECAYTYPE_ALPHA);

        if (i == nts) {
          mtot += cellmass;
        }

        for (const auto decaytype : decay::all_decaytypes) {
          // Qdot here has been multiplied by mass, so it is in units of [erg/s]
          const double qdot_cell = decay::get_qdot_modelcell(mgi, t_mid, decaytype) * cellmass;
          globals::timesteps[i].qdot_total += qdot_cell;
          if (decaytype == decay::DECAYTYPE_BETAMINUS) {
            globals::timesteps[i].qdot_betaminus += qdot_cell;
          } else if (decaytype == decay::DECAYTYPE_ALPHA) {
            globals::timesteps[i].qdot_alpha += qdot_cell;
          }
        }
      }
    }

#ifdef MPI_ON
    // in MPI mode, each process only did some fraction of the cells
    MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[i].eps_positron_ana_power, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[i].eps_electron_ana_power, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[i].eps_alpha_ana_power, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[i].qdot_betaminus, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[i].qdot_alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[i].qdot_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  }

#ifdef MPI_ON
  MPI_Allreduce(MPI_IN_PLACE, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (my_rank == 0) {
    FILE *dep_file = fopen_required("deposition.out.tmp", "w");
    fprintf(
        dep_file,
        "#ts tmid_days tmid_s total_dep_Lsun gammadep_Lsun gammadeppathint_Lsun positrondep_Lsun eps_positron_ana_Lsun "
        "elecdep_Lsun eps_elec_Lsun eps_elec_ana_Lsun alphadep_Lsun eps_alpha_Lsun eps_alpha_ana_Lsun eps_gamma_Lsun "
        "Qdot_betaminus_ana_erg/s/g Qdotalpha_ana_erg/s/g eps_erg/s/g Qdot_ana_erg/s/g\n");

    for (int i = 0; i <= nts; i++) {
      const double t_mid = globals::timesteps[i].mid;
      const double t_width = globals::timesteps[i].width;
      const double total_dep = (globals::timesteps[i].gamma_dep + globals::timesteps[i].positron_dep +
                                globals::timesteps[i].electron_dep + globals::timesteps[i].alpha_dep);

      // dep is used here for positrons and alphas because it is the same as the emission rate
      const double epsilon_mc = (globals::timesteps[i].gamma_emission + globals::timesteps[i].positron_dep +
                                 globals::timesteps[i].electron_emission + globals::timesteps[i].alpha_emission) /
                                mtot / t_width;

      fprintf(dep_file, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", i, t_mid / DAY, t_mid,
              total_dep / t_width / LSUN, globals::timesteps[i].gamma_dep / t_width / LSUN,
              globals::timesteps[i].gamma_dep_pathint / t_width / LSUN,
              globals::timesteps[i].positron_dep / t_width / LSUN, globals::timesteps[i].eps_positron_ana_power / LSUN,
              globals::timesteps[i].electron_dep / t_width / LSUN,
              globals::timesteps[i].electron_emission / t_width / LSUN,
              globals::timesteps[i].eps_electron_ana_power / LSUN, globals::timesteps[i].alpha_dep / t_width / LSUN,
              globals::timesteps[i].alpha_emission / t_width / LSUN, globals::timesteps[i].eps_alpha_ana_power / LSUN,
              globals::timesteps[i].gamma_emission / t_width / LSUN, globals::timesteps[i].qdot_betaminus / mtot,
              globals::timesteps[i].qdot_alpha / mtot, epsilon_mc, globals::timesteps[i].qdot_total / mtot);
    }
    fclose(dep_file);

    std::remove("deposition.out");
    std::rename("deposition.out.tmp", "deposition.out");
  }

  printout("calculating and writing deposition.out took %ld seconds\n",
           std::time(nullptr) - time_write_deposition_file_start);
}

#ifdef MPI_ON
void mpi_communicate_grid_properties(const int my_rank, const int nprocs, const int nstart, const int ndo,
                                     char *mpi_grid_buffer, const size_t mpi_grid_buffer_size) {
  int position = 0;
  for (int root = 0; root < nprocs; root++) {
    MPI_Barrier(MPI_COMM_WORLD);
    int root_nstart = nstart;
    MPI_Bcast(&root_nstart, 1, MPI_INT, root, MPI_COMM_WORLD);
    int root_ndo = ndo;
    MPI_Bcast(&root_ndo, 1, MPI_INT, root, MPI_COMM_WORLD);
    int root_node_id = globals::node_id;
    MPI_Bcast(&root_node_id, 1, MPI_INT, root, MPI_COMM_WORLD);

    for (int modelgridindex = root_nstart; modelgridindex < (root_nstart + root_ndo); modelgridindex++) {
      if (grid::get_numassociatedcells(modelgridindex) < 1) {
        continue;
      }

      radfield::do_MPI_Bcast(modelgridindex, root, root_node_id);

      nonthermal::nt_MPI_Bcast(modelgridindex, root, my_rank);

      if (globals::total_nlte_levels > 0 && globals::rank_in_node == 0) {
        MPI_Bcast(grid::modelgrid[modelgridindex].nlte_pops, globals::total_nlte_levels, MPI_DOUBLE, root_node_id,
                  globals::mpi_comm_internode);
      }

      if constexpr (USE_LUT_PHOTOION) {
        const auto nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
        assert_always(globals::corrphotoionrenorm != nullptr);
        if (globals::rank_in_node == 0) {
          MPI_Bcast(globals::corrphotoionrenorm + (nonemptymgi * globals::nbfcontinua_ground),
                    globals::nbfcontinua_ground, MPI_DOUBLE, root_node_id, globals::mpi_comm_internode);
        }

        assert_always(globals::gammaestimator != nullptr);
        MPI_Bcast(globals::gammaestimator + (nonemptymgi * globals::nbfcontinua_ground), globals::nbfcontinua_ground,
                  MPI_DOUBLE, root, MPI_COMM_WORLD);
      }

      assert_always(grid::modelgrid[modelgridindex].elem_meanweight != nullptr);
      if (globals::rank_in_node == 0) {
        MPI_Bcast(grid::modelgrid[modelgridindex].elem_meanweight, get_nelements(), MPI_FLOAT, root_node_id,
                  globals::mpi_comm_internode);
      }

      MPI_Bcast_binned_opacities(modelgridindex, root_node_id);
    }

    if (root == my_rank) {
      position = 0;
      MPI_Pack(&ndo, 1, MPI_INT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
      for (int mgi = nstart; mgi < (nstart + ndo); mgi++) {
        MPI_Pack(&mgi, 1, MPI_INT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);

        if (grid::get_numassociatedcells(mgi) > 0) {
          MPI_Pack(&grid::modelgrid[mgi].Te, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position,
                   MPI_COMM_WORLD);
          MPI_Pack(&grid::modelgrid[mgi].TR, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position,
                   MPI_COMM_WORLD);
          MPI_Pack(&grid::modelgrid[mgi].TJ, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position,
                   MPI_COMM_WORLD);
          MPI_Pack(&grid::modelgrid[mgi].W, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position,
                   MPI_COMM_WORLD);
          MPI_Pack(&grid::modelgrid[mgi].rho, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position,
                   MPI_COMM_WORLD);
          MPI_Pack(&grid::modelgrid[mgi].nne, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position,
                   MPI_COMM_WORLD);
          MPI_Pack(&grid::modelgrid[mgi].nnetot, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position,
                   MPI_COMM_WORLD);
          MPI_Pack(&grid::modelgrid[mgi].kappagrey, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position,
                   MPI_COMM_WORLD);
          MPI_Pack(&grid::modelgrid[mgi].totalcooling, 1, MPI_DOUBLE, mpi_grid_buffer, mpi_grid_buffer_size, &position,
                   MPI_COMM_WORLD);
          MPI_Pack(&grid::modelgrid[mgi].thick, 1, MPI_SHORT, mpi_grid_buffer, mpi_grid_buffer_size, &position,
                   MPI_COMM_WORLD);

          for (int element = 0; element < get_nelements(); element++) {
            if (get_nions(element) > 0) {
              MPI_Pack(grid::modelgrid[mgi].composition[element].groundlevelpop, get_nions(element), MPI_FLOAT,
                       mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
              MPI_Pack(grid::modelgrid[mgi].composition[element].partfunct, get_nions(element), MPI_FLOAT,
                       mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
              MPI_Pack(grid::modelgrid[mgi].cooling_contrib_ion[element], get_nions(element), MPI_DOUBLE,
                       mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
            }
          }
        }
      }
      printout("[info] mem_usage: MPI_BUFFER: used %d of %zu bytes allocated to mpi_grid_buffer\n", position,
               mpi_grid_buffer_size);
      assert_always(static_cast<size_t>(position) <= mpi_grid_buffer_size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(mpi_grid_buffer, mpi_grid_buffer_size, MPI_PACKED, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    position = 0;
    int nlp = 0;
    MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &nlp, 1, MPI_INT, MPI_COMM_WORLD);
    for (int nn = 0; nn < nlp; nn++) {
      int mgi = 0;
      MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &mgi, 1, MPI_INT, MPI_COMM_WORLD);

      if (grid::get_numassociatedcells(mgi) > 0) {
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &grid::modelgrid[mgi].Te, 1, MPI_FLOAT,
                   MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &grid::modelgrid[mgi].TR, 1, MPI_FLOAT,
                   MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &grid::modelgrid[mgi].TJ, 1, MPI_FLOAT,
                   MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &grid::modelgrid[mgi].W, 1, MPI_FLOAT,
                   MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &grid::modelgrid[mgi].rho, 1, MPI_FLOAT,
                   MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &grid::modelgrid[mgi].nne, 1, MPI_FLOAT,
                   MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &grid::modelgrid[mgi].nnetot, 1, MPI_FLOAT,
                   MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &grid::modelgrid[mgi].kappagrey, 1, MPI_FLOAT,
                   MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &grid::modelgrid[mgi].totalcooling, 1, MPI_DOUBLE,
                   MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &grid::modelgrid[mgi].thick, 1, MPI_SHORT,
                   MPI_COMM_WORLD);

        for (int element = 0; element < get_nelements(); element++) {
          if (get_nions(element) > 0) {
            MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position,
                       grid::modelgrid[mgi].composition[element].groundlevelpop, get_nions(element), MPI_FLOAT,
                       MPI_COMM_WORLD);
            MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position,
                       grid::modelgrid[mgi].composition[element].partfunct, get_nions(element), MPI_FLOAT,
                       MPI_COMM_WORLD);
            MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position,
                       grid::modelgrid[mgi].cooling_contrib_ion[element], get_nions(element), MPI_DOUBLE,
                       MPI_COMM_WORLD);
          }
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void mpi_reduce_estimators(int nts) {
  const int nonempty_npts_model = grid::get_nonempty_npts_model();
  radfield::reduce_estimators();
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, globals::ffheatingestimator, nonempty_npts_model, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, globals::colheatingestimator, nonempty_npts_model, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if constexpr (USE_LUT_PHOTOION) {
    MPI_Barrier(MPI_COMM_WORLD);
    assert_always(globals::gammaestimator != nullptr);
    MPI_Allreduce(MPI_IN_PLACE, globals::gammaestimator, nonempty_npts_model * globals::nbfcontinua_ground, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
  }

  if constexpr (USE_LUT_BFHEATING) {
    MPI_Barrier(MPI_COMM_WORLD);
    assert_always(globals::bfheatingestimator != nullptr);
    MPI_Allreduce(MPI_IN_PLACE, globals::bfheatingestimator, nonempty_npts_model * globals::nbfcontinua_ground,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  if constexpr (RECORD_LINESTAT) {
    MPI_Barrier(MPI_COMM_WORLD);
    assert_always(globals::ecounter != nullptr);
    MPI_Allreduce(MPI_IN_PLACE, globals::ecounter, globals::nlines, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    assert_always(globals::acounter != nullptr);
    MPI_Allreduce(MPI_IN_PLACE, globals::acounter, globals::nlines, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }

  assert_always(static_cast<int>(globals::dep_estimator_gamma.size()) == nonempty_npts_model);
  MPI_Allreduce(MPI_IN_PLACE, globals::dep_estimator_gamma.data(), nonempty_npts_model, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  /// Communicate gamma and positron deposition and write to file
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].cmf_lum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].gamma_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].positron_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].electron_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].electron_emission, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].alpha_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].alpha_emission, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::timesteps[nts].gamma_emission, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  globals::timesteps[nts].cmf_lum /= globals::nprocs;
  globals::timesteps[nts].gamma_dep /= globals::nprocs;
  globals::timesteps[nts].positron_dep /= globals::nprocs;
  globals::timesteps[nts].electron_dep /= globals::nprocs;
  globals::timesteps[nts].electron_emission /= globals::nprocs;
  globals::timesteps[nts].alpha_dep /= globals::nprocs;
  globals::timesteps[nts].alpha_emission /= globals::nprocs;
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
      printout("ERROR: Could not open file '%s' for mode 'wb'. \n", filename);
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

  if (access(filename, F_OK) == 0) {
    std::remove(filename);
    printout("Deleted %s\n", filename);
  }
}

void remove_grid_restart_data(const int timestep) {
  char prevfilename[MAXFILENAMELENGTH];
  snprintf(prevfilename, MAXFILENAMELENGTH, "gridsave_ts%d.tmp", timestep);

  if (access(prevfilename, F_OK) == 0) {
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

void save_grid_and_packets(const int nts, const int my_rank, const Packet *packets) {
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  bool write_verified_sucess = false;
  while (!write_verified_sucess) {
    const auto time_write_packets_file_start = std::time(nullptr);
    printout("time before write temporary packets file %ld\n", time_write_packets_file_start);

    // save packet state at start of current timestep (before propagation)
    write_temp_packetsfile(nts, my_rank, packets);

    vpkt_write_timestep(nts, my_rank, false);

    const auto time_write_packets_file_finished = std::time(nullptr);

#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    printout("time after write temporary packets file %ld (took %ld seconds, waited %ld s for other ranks)\n",
             std::time(nullptr), time_write_packets_file_finished - time_write_packets_file_start,
             std::time(nullptr) - time_write_packets_file_finished);

#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    const auto time_readback_packets_start = std::time(nullptr);

    printout("reading back temporary packets file to check validity...\n");

    // read packets file back to check that the disk write didn't fail
    write_verified_sucess = verify_temp_packetsfile(nts, my_rank, packets);

#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    printout("Verifying packets files for all ranks took %ld seconds.\n",
             std::time(nullptr) - time_readback_packets_start);
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

  std::fill_n(globals::ffheatingestimator, grid::get_nonempty_npts_model(), 0.);
  std::fill_n(globals::colheatingestimator, grid::get_nonempty_npts_model(), 0.);
  std::ranges::fill(globals::dep_estimator_gamma, 0.);

  if constexpr (USE_LUT_PHOTOION) {
    std::fill_n(globals::gammaestimator, grid::get_nonempty_npts_model() * globals::nbfcontinua_ground, 0.);
  }

  if constexpr (USE_LUT_BFHEATING) {
    std::fill_n(globals::bfheatingestimator, grid::get_nonempty_npts_model() * globals::nbfcontinua_ground, 0.);
  }

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void normalise_deposition_estimators(int nts) {
  const double dt = globals::timesteps[nts].width;
  globals::timesteps[nts].gamma_dep_pathint = 0.;
  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);

    const double dV = grid::get_modelcell_assocvolume_tmin(mgi) * pow(globals::timesteps[nts].mid / globals::tmin, 3);

    const double estimator_normfactor = 1 / dV / dt / globals::nprocs;

    globals::timesteps[nts].gamma_dep_pathint += globals::dep_estimator_gamma[nonemptymgi] / globals::nprocs;

    globals::dep_estimator_gamma[nonemptymgi] *= ONEOVER4PI / dV / dt / globals::nprocs;

    assert_testmodeonly(globals::dep_estimator_gamma[nonemptymgi] >= 0.);
    assert_testmodeonly(std::isfinite(globals::dep_estimator_gamma[nonemptymgi]));

    globals::dep_estimator_positron[nonemptymgi] *= estimator_normfactor;

    globals::dep_estimator_electron[nonemptymgi] *= estimator_normfactor;

    globals::dep_estimator_alpha[nonemptymgi] *= estimator_normfactor;
  }
}

auto do_timestep(const int nts, const int titer, const int my_rank, const int nstart, const int ndo, Packet *packets,
                 const int walltimelimitseconds) -> bool {
  bool do_this_full_loop = true;

  const int nts_prev = (titer != 0 || nts == 0) ? nts : nts - 1;
  if ((titer > 0) || (globals::simulation_continued_from_saved && (nts == globals::timestep_initial))) {
    /// Read the packets file to reset before each additional iteration on the timestep
    read_temp_packetsfile(nts, my_rank, packets);
  }

  /// Some counters on pkt-actions need to be reset to do statistics
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

  update_grid(estimators_file, nts, nts_prev, my_rank, nstart, ndo, titer, real_time_start);

  const auto sys_time_start_communicate_grid = std::time(nullptr);

/// Each process has now updated its own set of cells. The results now need to be communicated between processes.
#ifdef MPI_ON
  mpi_communicate_grid_properties(my_rank, globals::nprocs, nstart, ndo, mpi_grid_buffer, mpi_grid_buffer_size);
#endif

  printout("timestep %d: time after grid properties have been communicated %ld (took %ld seconds)\n", nts,
           std::time(nullptr), std::time(nullptr) - sys_time_start_communicate_grid);

  /// If this is not the 0th time step of the current job step,
  /// write out a snapshot of the grid properties for further restarts
  /// and update input.txt accordingly
  if (((nts - globals::timestep_initial) != 0)) {
    save_grid_and_packets(nts, my_rank, packets);
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
    /// Now process the packets.

    update_packets(my_rank, nts, std::span{packets, static_cast<size_t>(globals::npkts)});

#ifdef MPI_ON
    // All the processes have their own versions of the estimators for this time step now.
    // Since these are going to be needed in the next time step, we will gather all the
    // estimators together now, sum them, and distribute the results

    const auto time_communicate_estimators_start = std::time(nullptr);
    mpi_reduce_estimators(nts);
#endif

    // The estimators have been summed across all proceses and distributed.
    // They will now be normalised independently on all processes

    normalise_deposition_estimators(nts);

    write_deposition_file(nts, my_rank, nstart, ndo);

    write_partial_lightcurve_spectra(my_rank, nts, packets);

#ifdef MPI_ON
    printout("timestep %d: time after estimators have been communicated %ld (took %ld seconds)\n", nts,
             std::time(nullptr), std::time(nullptr) - time_communicate_estimators_start);
#endif

    printout("During timestep %d on MPI process %d, %d pellets decayed and %d packets escaped. (t=%gd)\n", nts, my_rank,
             globals::timesteps[nts].pellet_decays, globals::nesc, globals::timesteps[nts].mid / DAY);

    if (VPKT_ON) {
      printout("During timestep %d on MPI process %d, %d virtual packets were generated and %d escaped. \n", nts,
               my_rank, nvpkt, nvpkt_esc1 + nvpkt_esc2 + nvpkt_esc3);
      printout(
          "%d virtual packets came from an electron scattering event, %d from a kpkt deactivation and %d from a "
          "macroatom deactivation. \n",
          nvpkt_esc1, nvpkt_esc2, nvpkt_esc3);

      nvpkt = 0;
      nvpkt_esc1 = 0;
      nvpkt_esc2 = 0;
      nvpkt_esc3 = 0;
    }

    if constexpr (RECORD_LINESTAT) {
      if (my_rank == 0) {
        /// Print net absorption/emission in lines to the linestat_file
        /// Currently linestat information is only properly implemented for MPI only runs
        /// For hybrid runs only data from thread 0 is recorded
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

  const int my_rank = globals::rank_global;

#if defined(_OPENMP) && !defined(GPU_ON)
  // Explicitly turn off dynamic threads because we use the threadprivate directive!!!
  omp_set_dynamic(0);
#pragma omp parallel private(filename)
#endif
  {
    /// initialise the thread and rank specific output file
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

  auto *const packets = static_cast<Packet *>(malloc(MPKTS * sizeof(Packet)));

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
  printout("  rank %d of [0..%d] in MPI_COMM_WORLD\n", globals::rank_global, globals::nprocs - 1);
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

  /// Precalculate the rate coefficients for spontaneous and stimulated recombination
  /// and for photoionisation. With the nebular approximation they only depend on T_e
  /// T_R and W. W is easily factored out. For stimulated recombination we must assume
  /// T_e = T_R for this precalculation.

  ratecoefficients_init();

#ifdef MPI_ON
  printout("barrier after tabulation of rate coefficients: time before barrier %ld, ", std::time(nullptr));
  MPI_Barrier(MPI_COMM_WORLD);
  printout("time after barrier %ld\n", std::time(nullptr));
#endif

  stats::init();

  /// Record the chosen syn_dir
  FILE *syn_file = fopen_required("syn_dir.txt", "w");
  fprintf(syn_file, "%g %g %g", globals::syn_dir[0], globals::syn_dir[1], globals::syn_dir[2]);
  fclose(syn_file);
  printout("time write syn_dir.txt file %ld\n", std::time(nullptr));

  bool terminate_early = false;

  time_init();

  if (my_rank == 0) {
    write_timestep_file();
  }

  /// Initialise the grid. Set up the initial positions and sizes of the grid cells.
  printout("time grid_init %ld\n", std::time(nullptr));
  grid::grid_init(my_rank);

  printout("Simulation propagates %g packets per process (total %g with nprocs %d)\n", 1. * globals::npkts,
           1. * globals::npkts * globals::nprocs, globals::nprocs);

  printout("[info] mem_usage: packets occupy %.3f MB\n", MPKTS * sizeof(Packet) / 1024. / 1024.);

  if (!globals::simulation_continued_from_saved) {
    std::remove("deposition.out");
    /// Next we want to initialise the packets.
    /// Create a bunch of npkts packets
    /// and write them to a binary file for later readin.
    packet_init(packets);
    zero_estimators();
  }

  /// For the parallelisation of update_grid, the process needs to be told which cells belong to it.
  /// The next loop is over all grid cells. For parallelisation, we want to split this loop between
  /// processes. This is done by assigning each MPI process nblock cells. The residual n_leftover
  /// cells are sent to processes 0 ... process n_leftover -1.
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
  const size_t maxndo = grid::get_maxndo();
  /// Initialise the exchange buffer
  /// The factor 4 comes from the fact that our buffer should contain elements of 4 byte
  /// instead of 1 byte chars. But the MPI routines don't care about the buffers datatype
  mpi_grid_buffer_size = 4 * ((12 + 4 * get_includedions()) * (maxndo) + 1);
  printout("reserve mpi_grid_buffer_size %zu space for MPI communication buffer\n", mpi_grid_buffer_size);
  mpi_grid_buffer = static_cast<char *>(malloc(mpi_grid_buffer_size * sizeof(char)));
  assert_always(mpi_grid_buffer != nullptr);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  int nts = globals::timestep_initial;

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
  vpkt_init(nts, my_rank, globals::simulation_continued_from_saved);

  while (nts < globals::timestep_finish && !terminate_early) {
    globals::timestep = nts;
#ifdef MPI_ON
    //        const auto time_before_barrier = std::time(nullptr);
    MPI_Barrier(MPI_COMM_WORLD);
    //        const auto time_after_barrier = std::time(nullptr);
    //        printout("timestep %d: time before barrier %d, time after barrier %d\n", nts, time_before_barrier,
    //        time_after_barrier);
#endif

    /// titer example: Do 3 iterations on timestep 0-6
    // globals::n_titer = (nts < 6) ? 3: 1;

    globals::n_titer = 1;
    globals::lte_iteration = (nts < globals::num_lte_timesteps);
    assert_always(globals::num_lte_timesteps > 0);  // The first time step must solve the ionisation balance in LTE

    for (int titer = 0; titer < globals::n_titer; titer++) {
      terminate_early = do_timestep(nts, titer, my_rank, nstart, ndo, packets, walltimelimitseconds);
#ifdef DO_TITER
      /// No iterations over the zeroth timestep, set titer > n_titer
      if (nts == 0) titer = globals::n_titer + 1;
#endif
    }

    nts++;
  }

  /// The main calculation is now over. The packets now have all stored the time, place and direction
  /// at which they left the grid. Also their rest frame energies and frequencies.
  /// Spectra and light curves are now extracted using exspec which is another make target of this
  /// code.

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
  free(mpi_grid_buffer);
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
