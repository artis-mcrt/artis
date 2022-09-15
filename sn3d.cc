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

#include <getopt.h>
#include <unistd.h>

#include "atomic.h"
#include "decay.h"
#include "emissivities.h"
#include "globals.h"
#include "grey_emissivities.h"
#include "grid.h"
#include "input.h"
// #include "ltepop.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "spectrum.h"
#include "stats.h"
#include "update_grid.h"
#include "update_packets.h"
#include "version.h"
#include "vpkt.h"

const bool KEEP_ALL_RESTART_FILES =
    false;  // once a new gridsave and packets*.tmp have been written, don't delete the previous set

const bool do_exspec = false;

// threadprivate variables
int tid;
__managed__ int myGpuId = 0;
__managed__ bool use_cellhist;
__managed__ bool neutral_flag;
#ifndef __CUDA_ARCH__
gsl_rng *rng = NULL;
#else
__device__ void *rng = NULL;
#endif
gsl_integration_workspace *gslworkspace = NULL;
FILE *output_file = NULL;
static FILE *linestat_file = NULL;
static time_t real_time_start = -1;
static time_t time_timestep_start = -1;  // this will be set after the first update of the grid and before packet prop
static FILE *estimators_file = NULL;

int mpi_grid_buffer_size = 0;
char *mpi_grid_buffer = NULL;

static void initialise_linestat_file(void) {
  linestat_file = fopen_required("linestat.out", "w");

  for (int i = 0; i < globals::nlines; i++) fprintf(linestat_file, "%g ", CLIGHT / globals::linelist[i].nu);
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++)
    fprintf(linestat_file, "%d ", get_element(globals::linelist[i].elementindex));
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++)
    fprintf(linestat_file, "%d ", get_ionstage(globals::linelist[i].elementindex, globals::linelist[i].ionindex));
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++) fprintf(linestat_file, "%d ", globals::linelist[i].upperlevelindex + 1);
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++) fprintf(linestat_file, "%d ", globals::linelist[i].lowerlevelindex + 1);
  fprintf(linestat_file, "\n");

  fflush(linestat_file);
  // setvbuf(linestat_file, NULL, _IOLBF, 1); // flush after every line makes it slow!
}

static void write_deposition_file(const int nts, const int my_rank, const int nstart, const int ndo) {
  printout("Calculating deposition rates...\n");
  time_t time_write_deposition_file_start = time(NULL);
  double mtot = 0.;

  // calculate analytical decay rates
  // for (int i = 0; i <= nts; i++)
  const int i = nts;
  {
    const double t_mid = globals::time_step[i].mid;

    // power in [erg/s]
    globals::time_step[i].eps_positron_ana_power = 0.;
    globals::time_step[i].eps_electron_ana_power = 0.;
    globals::time_step[i].eps_alpha_ana_power = 0.;
    globals::time_step[i].qdot_betaminus = 0.;
    globals::time_step[i].qdot_alpha = 0.;
    globals::time_step[i].qdot_total = 0.;

    for (int mgi = nstart; mgi < (nstart + ndo); mgi++)
    // for (int mgi = 0; mgi < grid::get_npts_model(); mgi++)
    {
      if (grid::get_numassociatedcells(mgi) > 0) {
        const double cellmass = grid::get_rhoinit(mgi) * grid::vol_init_modelcell(mgi);

        globals::time_step[i].eps_positron_ana_power +=
            cellmass * decay::get_particle_injection_rate(mgi, t_mid, decay::DECAYTYPE_BETAPLUS);
        globals::time_step[i].eps_electron_ana_power +=
            cellmass * decay::get_particle_injection_rate(mgi, t_mid, decay::DECAYTYPE_BETAMINUS);
        globals::time_step[i].eps_alpha_ana_power +=
            cellmass * decay::get_particle_injection_rate(mgi, t_mid, decay::DECAYTYPE_ALPHA);

        if (i == nts) {
          mtot += cellmass;
        }

        for (int dectypeindex = 0; dectypeindex < decay::DECAYTYPE_COUNT; dectypeindex++) {
          // Qdot here has been multiplied by mass, so it is in units of [erg/s]
          const double qdot_cell = decay::get_qdot_modelcell(mgi, t_mid, dectypeindex) * cellmass;
          globals::time_step[i].qdot_total += qdot_cell;
          if (dectypeindex == decay::DECAYTYPE_BETAMINUS) {
            globals::time_step[i].qdot_betaminus += qdot_cell;
          } else if (dectypeindex == decay::DECAYTYPE_ALPHA) {
            globals::time_step[i].qdot_alpha += qdot_cell;
          }
        }
      }
    }

#ifdef MPI_ON
    // in MPI mode, each process only did some fraction of the cells
    MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[i].eps_positron_ana_power, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[i].eps_electron_ana_power, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[i].eps_alpha_ana_power, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[i].qdot_betaminus, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[i].qdot_alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[i].qdot_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
      const double t_mid = globals::time_step[i].mid;
      const double t_width = globals::time_step[i].width;
      const double total_dep = (globals::time_step[i].gamma_dep + globals::time_step[i].positron_dep +
                                globals::time_step[i].electron_dep + globals::time_step[i].alpha_dep);

      // dep is used here for positrons and alphas because it is the same as the emission rate
      const double epsilon_mc = (globals::time_step[i].gamma_emission + globals::time_step[i].positron_dep +
                                 globals::time_step[i].electron_emission + globals::time_step[i].alpha_emission) /
                                mtot / t_width;

      fprintf(dep_file, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", i, t_mid / DAY, t_mid,
              total_dep / t_width / LSUN, globals::time_step[i].gamma_dep / t_width / LSUN,
              globals::time_step[i].gamma_dep_pathint / t_width / LSUN,
              globals::time_step[i].positron_dep / t_width / LSUN, globals::time_step[i].eps_positron_ana_power / LSUN,
              globals::time_step[i].electron_dep / t_width / LSUN,
              globals::time_step[i].electron_emission / t_width / LSUN,
              globals::time_step[i].eps_electron_ana_power / LSUN, globals::time_step[i].alpha_dep / t_width / LSUN,
              globals::time_step[i].alpha_emission / t_width / LSUN, globals::time_step[i].eps_alpha_ana_power / LSUN,
              globals::time_step[i].gamma_emission / t_width / LSUN, globals::time_step[i].qdot_betaminus / mtot,
              globals::time_step[i].qdot_alpha / mtot, epsilon_mc, globals::time_step[i].qdot_total / mtot);
    }
    fclose(dep_file);

    std::remove("deposition.out");
    std::rename("deposition.out.tmp", "deposition.out");
  }

  printout("calculating and writing deposition.out took %ld seconds\n", time(NULL) - time_write_deposition_file_start);
}

#ifdef MPI_ON
static void mpi_communicate_grid_properties(const int my_rank, const int nprocs, const int nstart, const int ndo,
                                            const int nts, const int titer, char *mpi_grid_buffer,
                                            const int mpi_grid_buffer_size) {
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
      radfield::do_MPI_Bcast(modelgridindex, root, root_node_id);

      if (grid::get_numassociatedcells(modelgridindex) > 0) {
        nonthermal::nt_MPI_Bcast(modelgridindex, root);
        if (NLTE_POPS_ON && globals::rank_in_node == 0) {
          MPI_Bcast(grid::modelgrid[modelgridindex].nlte_pops, globals::total_nlte_levels, MPI_DOUBLE, root_node_id,
                    globals::mpi_comm_internode);
        }

#if (!NO_LUT_PHOTOION)
        assert_always(globals::corrphotoionrenorm != NULL);
        MPI_Bcast(&globals::corrphotoionrenorm[modelgridindex * get_nelements() * get_max_nions()],
                  get_nelements() * get_max_nions(), MPI_DOUBLE, root, MPI_COMM_WORLD);
        assert_always(globals::gammaestimator != NULL);
        MPI_Bcast(&globals::gammaestimator[modelgridindex * get_nelements() * get_max_nions()],
                  get_nelements() * get_max_nions(), MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif
      }
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
            MPI_Pack(grid::modelgrid[mgi].composition[element].groundlevelpop, get_nions(element), MPI_FLOAT,
                     mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
            MPI_Pack(grid::modelgrid[mgi].composition[element].partfunct, get_nions(element), MPI_FLOAT,
                     mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
            MPI_Pack(grid::modelgrid[mgi].cooling_contrib_ion[element], get_nions(element), MPI_DOUBLE, mpi_grid_buffer,
                     mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          }
        }
      }
      printout("[info] mem_usage: MPI_BUFFER: used %d of %d bytes allocated to mpi_grid_buffer\n", position,
               mpi_grid_buffer_size);
      assert_always(position <= mpi_grid_buffer_size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(mpi_grid_buffer, mpi_grid_buffer_size, MPI_PACKED, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    position = 0;
    int nlp;
    MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &nlp, 1, MPI_INT, MPI_COMM_WORLD);
    for (int nn = 0; nn < nlp; nn++) {
      int mgi;
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
          MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position,
                     grid::modelgrid[mgi].composition[element].groundlevelpop, get_nions(element), MPI_FLOAT,
                     MPI_COMM_WORLD);
          MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position,
                     grid::modelgrid[mgi].composition[element].partfunct, get_nions(element), MPI_FLOAT,
                     MPI_COMM_WORLD);
          MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position,
                     grid::modelgrid[mgi].cooling_contrib_ion[element], get_nions(element), MPI_DOUBLE, MPI_COMM_WORLD);
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

static void mpi_reduce_estimators(int my_rank, int nts) {
  radfield::reduce_estimators();
#ifndef FORCE_LTE
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, globals::ffheatingestimator, grid::get_npts_model(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, globals::colheatingestimator, grid::get_npts_model(), MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#if (!NO_LUT_PHOTOION) || (!NO_LUT_BFHEATING)
  const int arraylen = grid::get_npts_model() * get_nelements() * get_max_nions();
#endif
#if (!NO_LUT_PHOTOION)
  MPI_Barrier(MPI_COMM_WORLD);
  assert_always(globals::gammaestimator != NULL);
  MPI_Allreduce(MPI_IN_PLACE, globals::gammaestimator, arraylen, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#if (!NO_LUT_BFHEATING)
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, globals::bfheatingestimator, arraylen, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#endif

#ifdef RECORD_LINESTAT
  MPI_Barrier(MPI_COMM_WORLD);
  assert_always(globals::ecounter != NULL);
  MPI_Allreduce(MPI_IN_PLACE, globals::ecounter, globals::nlines, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  assert_always(globals::acounter != NULL);
  MPI_Allreduce(MPI_IN_PLACE, globals::acounter, globals::nlines, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  // double deltaV = pow(grid::wid_init * globals::time_step[nts].mid/globals::tmin, 3.0);
  // double deltat = globals::time_step[nts].width;
  if (globals::do_rlc_est != 0) {
    assert_always(globals::rpkt_emiss != NULL);
    MPI_Allreduce(MPI_IN_PLACE, globals::rpkt_emiss, grid::get_npts_model(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (globals::do_comp_est) {
    assert_always(globals::compton_emiss != NULL);
    MPI_Allreduce(MPI_IN_PLACE, globals::compton_emiss, grid::get_npts_model() * EMISS_MAX, MPI_FLOAT, MPI_SUM,
                  MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  /// Communicate gamma and positron deposition and write to file
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].cmf_lum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].gamma_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].positron_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].electron_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].electron_emission, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].alpha_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].alpha_emission, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].gamma_emission, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  globals::time_step[nts].cmf_lum /= globals::nprocs;
  globals::time_step[nts].gamma_dep /= globals::nprocs;
  globals::time_step[nts].positron_dep /= globals::nprocs;
  globals::time_step[nts].electron_dep /= globals::nprocs;
  globals::time_step[nts].electron_emission /= globals::nprocs;
  globals::time_step[nts].alpha_dep /= globals::nprocs;
  globals::time_step[nts].alpha_emission /= globals::nprocs;
  globals::time_step[nts].gamma_emission /= globals::nprocs;

#if TRACK_ION_STATS
  stats::reduce_estimators();
#endif

  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

static void write_temp_packetsfile(const int timestep, const int my_rank, const struct packet *const pkt) {
  // write packets binary file
  char filename[128];
  snprintf(filename, 128, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  printout("Writing %s...", filename);
  FILE *packets_file = fopen_required(filename, "wb");

  fwrite(pkt, sizeof(struct packet), globals::npkts, packets_file);
  fclose(packets_file);
  printout("done\n");
}

static void remove_temp_packetsfile(const int timestep, const int my_rank) {
  char filename[128];
  snprintf(filename, 128, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  if (!access(filename, F_OK)) {
    remove(filename);
    printout("Deleted %s\n", filename);
  }
}

static void remove_grid_restart_data(const int timestep) {
  char prevfilename[128];
  snprintf(prevfilename, 128, "gridsave_ts%d.tmp", timestep);

  if (!access(prevfilename, F_OK)) {
    remove(prevfilename);
    printout("Deleted %s\n", prevfilename);
  }
}

static bool walltime_sufficient_to_continue(const int nts, const int nts_prev, const int walltimelimitseconds) {
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // time is measured from just before packet propagation from one timestep to the next
  const int estimated_time_per_timestep = time(NULL) - time_timestep_start;
  printout("TIME: time between timesteps is %d seconds (measured packet prop of ts %d and update grid of ts %d)\n",
           estimated_time_per_timestep, nts_prev, nts);

  bool do_this_full_loop = true;
  if (walltimelimitseconds > 0) {
    const int wallclock_used_seconds = time(NULL) - real_time_start;
    const int wallclock_remaining_seconds = walltimelimitseconds - wallclock_used_seconds;
    printout("TIMED_RESTARTS: Used %d of %d seconds of wall time.\n", wallclock_used_seconds, walltimelimitseconds);

    // This flag being false will make it update_grid, and then exit
    do_this_full_loop = (wallclock_remaining_seconds >= (1.5 * estimated_time_per_timestep));

#ifdef MPI_ON
    // communicate whatever decision the rank 0 process decided, just in case they differ
    MPI_Bcast(&do_this_full_loop, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif
    if (do_this_full_loop)
      printout("TIMED_RESTARTS: Going to continue since remaining time %d s >= 1.5 * time_per_timestep\n",
               wallclock_remaining_seconds);
    else
      printout("TIMED_RESTARTS: Going to terminate since remaining time %d s < 1.5 * time_per_timestep\n",
               wallclock_remaining_seconds);
  }
  return do_this_full_loop;
}

#if CUDA_ENABLED
void *makemanaged(void *ptr, size_t curSize) {
  if (ptr == NULL) {
    return NULL;
  }
  void *newptr;
  cudaMallocManaged(&newptr, curSize);
  memcpy(newptr, ptr, curSize);
  free(ptr);
  return newptr;
}
#endif

static void save_grid_and_packets(const int nts, const int my_rank, struct packet *packets) {
  const time_t time_write_packets_file_start = time(NULL);
  printout("time before write temporary packets file %ld\n", time_write_packets_file_start);

  // save packet state at start of current timestep (before propagation)
  write_temp_packetsfile(nts, my_rank, packets);

#ifdef VPKT_ON
  char filename[128];
  snprintf(filename, 128, "vspecpol_%d_%d_%s.tmp", 0, my_rank, (nts % 2 == 0) ? "even" : "odd");

  FILE *vspecpol_file = fopen_required(filename, "wb");

  write_vspecpol(vspecpol_file);
  fclose(vspecpol_file);

  // Write temporary files for vpkt_grid
  if (vgrid_flag == 1) {
    snprintf(filename, 128, "vpkt_grid_%d_%d_%s.tmp", 0, my_rank, (nts % 2 == 0) ? "even" : "odd");

    FILE *vpkt_grid_file = fopen_required(filename, "wb");

    write_vpkt_grid(vpkt_grid_file);
    fclose(vpkt_grid_file);
  }
#endif

  printout("time after write temporary packets file %ld (took %ld seconds)\n", time(NULL),
           time(NULL) - time_write_packets_file_start);

  if (my_rank == 0) {
    grid::write_grid_restart_data(nts);
    update_parameterfile(nts);
  }

  if (!KEEP_ALL_RESTART_FILES) {
// ensure new packets files have been written by all processes before we remove the old set
#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (my_rank == 0) remove_grid_restart_data(nts - 1);

    // delete temp packets files from previous timestep now that all restart data for the new timestep is available
    remove_temp_packetsfile(nts - 1, my_rank);
  }
}

static bool do_timestep(const int nts, const int titer, const int my_rank, const int nstart, const int ndo,
                        struct packet *packets, const int walltimelimitseconds) {
  bool do_this_full_loop = true;

  const int nts_prev = (titer != 0 || nts == 0) ? nts : nts - 1;
  if ((titer > 0) || (globals::simulation_continued_from_saved && (nts == globals::itstep))) {
    /// Read the packets file to reset before each additional iteration on the timestep
    read_temp_packetsfile(nts, my_rank, packets);
  }

  /// Some counters on pkt-actions need to be reset to do statistics
  stats::pkt_action_counters_reset();

  if (nts == 0) {
    radfield::initialise_prev_titer_photoionestimators();
  }

#ifdef RECORD_LINESTAT
  // The same for absorption/emission of r-pkts in lines
  for (int i = 0; i < globals::nlines; i++) {
    globals::acounter[i] = 0;
    globals::ecounter[i] = 0;
  }
#endif

  globals::do_comp_est = globals::do_r_lc ? false : estim_switch(nts);

  // Update the matter quantities in the grid for the new timestep.

  update_grid(estimators_file, nts, nts_prev, my_rank, nstart, ndo, titer, real_time_start);

  const time_t sys_time_start_communicate_grid = time(NULL);

/// Each process has now updated its own set of cells. The results now need to be communicated between processes.
#ifdef MPI_ON
  mpi_communicate_grid_properties(my_rank, globals::nprocs, nstart, ndo, nts, titer, mpi_grid_buffer,
                                  mpi_grid_buffer_size);
#endif

  printout("timestep %d: time after grid properties have been communicated %ld (took %ld seconds)\n", nts, time(NULL),
           time(NULL) - sys_time_start_communicate_grid);

  /// If this is not the 0th time step of the current job step,
  /// write out a snapshot of the grid properties for further restarts
  /// and update input.txt accordingly
  if (((nts - globals::itstep) != 0)) {
    save_grid_and_packets(nts, my_rank, packets);
    do_this_full_loop = walltime_sufficient_to_continue(nts, nts_prev, walltimelimitseconds);
  }
  time_timestep_start = time(NULL);

  // set all the estimators to zero before moving packets. This is now done
  // after update_grid so that, if requires, the gamma-ray heating estimator is known there
  // and also the photoion and stimrecomb estimators
  zero_estimators();

  // MPI_Barrier(MPI_COMM_WORLD);
  if ((nts < globals::ftstep) && do_this_full_loop) {
    /// Now process the packets.

    update_packets(my_rank, nts, packets);

#ifdef MPI_ON
    // All the processes have their own versions of the estimators for this time step now.
    // Since these are going to be needed in the next time step, we will gather all the
    // estimators together now, sum them, and distribute the results

    const time_t time_communicate_estimators_start = time(NULL);
    mpi_reduce_estimators(my_rank, nts);
#endif

    // The estimators have been summed across all proceses and distributed.
    // They will now be normalised independently on all processes
    if (globals::do_comp_est) {
      normalise_compton_estimators(nts);
      if (my_rank == 0) {
        write_compton_estimators(nts);
      }
    }

    if (globals::do_rlc_est != 0) {
      normalise_grey(nts);
      if ((globals::do_rlc_est != 3) && (my_rank == 0)) {
        write_grey(nts);
      }
    }

    write_deposition_file(nts, my_rank, nstart, ndo);

    write_partial_lightcurve_spectra(my_rank, nts, packets);

#ifdef MPI_ON
    printout("timestep %d: time after estimators have been communicated %ld (took %ld seconds)\n", nts, time(NULL),
             time(NULL) - time_communicate_estimators_start);
#endif

    printout("During timestep %d on MPI process %d, %d pellets decayed and %d packets escaped. (t=%gd)\n", nts, my_rank,
             globals::time_step[nts].pellet_decays, globals::nesc, globals::time_step[nts].mid / DAY);

#ifdef VPKT_ON
    printout("During timestep %d on MPI process %d, %d virtual packets were generated and %d escaped. \n", nts, my_rank,
             nvpkt, nvpkt_esc1 + nvpkt_esc2 + nvpkt_esc3);
    printout(
        "%d virtual packets came from an electron scattering event, %d from a kpkt deactivation and %d from a "
        "macroatom deactivation. \n",
        nvpkt_esc1, nvpkt_esc2, nvpkt_esc3);

    nvpkt = 0;
    nvpkt_esc1 = 0;
    nvpkt_esc2 = 0;
    nvpkt_esc3 = 0;
#endif

#ifdef RECORD_LINESTAT
    if (my_rank == 0) {
      /// Print net absorption/emission in lines to the linestat_file
      /// Currently linestat information is only properly implemented for MPI only runs
      /// For hybrid runs only data from thread 0 is recorded
      for (int i = 0; i < globals::nlines; i++) fprintf(linestat_file, "%d ", globals::ecounter[i]);
      fprintf(linestat_file, "\n");
      for (int i = 0; i < globals::nlines; i++) fprintf(linestat_file, "%d ", globals::acounter[i]);
      fprintf(linestat_file, "\n");
      fflush(linestat_file);
    }
#endif

    if (nts == globals::ftstep - 1) {
      char filename[128];
      snprintf(filename, 128, "packets%.2d_%.4d.out", 0, my_rank);
      // snprintf(filename, 128, "packets%.2d_%.4d.out", middle_iteration, my_rank);
      write_packets(filename, packets);

// write specpol of the virtual packets
#ifdef VPKT_ON
      snprintf(filename, 128, "vspecpol_%d-%d.out", my_rank, tid);
      FILE *vspecpol_file = fopen_required(filename, "w");

      write_vspecpol(vspecpol_file);
      fclose(vspecpol_file);

      if (vgrid_flag == 1) {
        snprintf(filename, 128, "vpkt_grid_%d-%d.out", my_rank, tid);
        FILE *vpkt_grid_file = fopen_required(filename, "w");
        write_vpkt_grid(vpkt_grid_file);
        fclose(vpkt_grid_file);
      }
#endif

      printout("time after write final packets file %ld\n", time(NULL));

      // final packets*.out have been written, so remove the temporary packets files
      // commented out because you might still want to resume the simulation
      // snprintf(filename, 128, "packets%d_%d_odd.tmp", 0, my_rank);
      // remove(filename);
      // snprintf(filename, 128, "packets%d_%d_even.tmp", 0, my_rank);
      // remove(filename);
    }
  }
  return !do_this_full_loop;
}

int main(int argc, char **argv)
// Main - top level routine.
{
  // FILE *temperature_file;
  char filename[128];

#ifdef VPKT_ON
  nvpkt = 0;
  nvpkt_esc1 = 0;
  nvpkt_esc2 = 0;
  nvpkt_esc3 = 0;
#endif

#if CUDA_ENABLED
  printf("CUDA ENABLED\n");
  int deviceCount = -1;
  checkCudaErrors(cudaGetDeviceCount(&deviceCount));
  printf("deviceCount %d\n", deviceCount);
  assert_always(deviceCount > 0);
  {
    myGpuId = 0;
#ifdef _OPENMP
    const int omp_tid = get_thread_num();
    myGpuId = omp_tid % deviceCount;
    printf("OMP thread id %d\n", omp_tid);
#else
#if MPI_ON
    const int local_rank = atoi(getenv("OMPI_COMM_WORLD_LOCAL_RANK"));
    myGpuId = local_rank % deviceCount;
    printf("local_rank %d\n", local_rank);
#endif
#endif
    printf("myGpuId %d\n", myGpuId);
  }
  cudaSetDevice(myGpuId);
#endif

#ifdef MPI_ON
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &globals::rank_global);
  MPI_Comm_size(MPI_COMM_WORLD, &globals::nprocs);

  // make an intra-node communicator (group ranks that can share memory)
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

#else
  globals::rank_global = 0;
  globals::nprocs = 1;
  globals::rank_in_node = 0;
  globals::node_nprocs = 1;
  globals::node_id = 0;
  globals::node_count = 0;
#endif

  const int my_rank = globals::rank_global;

#ifdef _OPENMP
  /// Explicitly turn off dynamic threads because we use the threadprivate directive!!!
  omp_set_dynamic(0);

#pragma omp parallel private(filename)
#endif
  {
    /// Get the current threads ID, copy it to a threadprivate variable
    tid = get_thread_num();
    /// and initialise the threads outputfile
    snprintf(filename, 128, "output_%d-%d.txt", my_rank, tid);
    output_file = fopen_required(filename, "w");
    /// Makes sure that the output_file is written line-by-line
    setvbuf(output_file, NULL, _IOLBF, 1);

#ifdef _OPENMP
    printout("OpenMP parallelisation active with %d threads (max %d)\n", get_num_threads(), get_max_threads());
#else
    printout("OpenMP is not available in this build\n");
#endif

    gslworkspace = gsl_integration_workspace_alloc(GSLWSIZE);
  }

  real_time_start = time(NULL);
  printout("time at start %ld\n", real_time_start);

#ifdef WALLTIMELIMITSECONDS
  int walltimelimitseconds = WALLTIMELIMITSECONDS;
#else
  int walltimelimitseconds = -1;
#endif

  int opt;
  while ((opt = getopt(argc, argv, "w:")) != -1) {
    switch (opt) {
      case 'w': {
        printout("Command line argument specifies wall time hours '%s', setting ", optarg);
        const float walltimehours = strtof(optarg, NULL);
        walltimelimitseconds = walltimehours * 3600;
        printout("walltimelimitseconds = %d\n", walltimelimitseconds);
        break;
      }

        // case 'p':
        //   printout("Command line argument specifies number of packets '%s', setting ", optarg);
        //   npkts = (int) strtol(optarg, NULL, 10);
        //   printout("npkts = %d\n", npkts);
        //   break;

      default: {
        fprintf(stderr, "Usage: %s [-w WALLTIMELIMITHOURS]\n", argv[0]);
        exit(EXIT_FAILURE);
      }
    }
  }

#if CUDA_ENABLED
  struct packet *packets;
  cudaMallocManaged(&packets, MPKTS * sizeof(struct packet));
#if USECUDA_UPDATEPACKETS
  cudaMemAdvise(packets, MPKTS * sizeof(struct packet), cudaMemAdviseSetPreferredLocation, myGpuId);
#endif
#else
  struct packet *const packets = (struct packet *)calloc(MPKTS, sizeof(struct packet));
#endif

  assert_always(packets != NULL);

#ifndef GIT_BRANCH
#define GIT_BRANCH "UNKNOWN"
#endif
  printout("ARTIS git branch %s\n", GIT_BRANCH);

#ifndef GIT_VERSION
#define GIT_VERSION "UNKNOWN"
#endif
  printout("Current version: %s\n", GIT_VERSION);

  // printout("Hash of most recent commit: %s\n",GIT_HASH);
  printout("sn3d.cc compiled at %s on %s\n", __TIME__, __DATE__);

#if defined TESTMODE && TESTMODE
  printout("TESTMODE is ON");
#endif

#ifdef MPI_ON
  printout("MPI enabled:\n");
  printout("  rank %d of [0..%d] in MPI_COMM_WORLD\n", globals::rank_global, globals::nprocs - 1);
  printout("  node %d of [0..%d]\n", globals::node_id, globals::node_count - 1);
  printout("  rank %d of [0..%d] within this node (MPI_COMM_WORLD_SHARED)\n", globals::rank_in_node,
           globals::node_nprocs - 1);
#else
  printout("MPI is disabled in this build\n");
#endif

#ifdef __CUDACC__
  printout("[CUDA] NVIDIA CUDA is available in this build\n");
#endif

#if CUDA_ENABLED
  printout("[CUDA] NVIDIA CUDA accelerated routines are enabled\n");
  printout("[CUDA] Detected %d CUDA capable device(s)\n", deviceCount);
  printout("[CUDA] This rank will use cudaSetDevice(%d)\n", myGpuId);
  struct cudaDeviceProp deviceProp;
  checkCudaErrors(cudaGetDeviceProperties(&deviceProp, myGpuId));
  printout("[CUDA] totalGlobalMem %7.1f GiB\n", deviceProp.totalGlobalMem / 1024. / 1024. / 1024.);
  printout("[CUDA] multiProcessorCount %d\n", deviceProp.multiProcessorCount);
  printout("[CUDA] maxThreadsPerMultiProcessor %d\n", deviceProp.maxThreadsPerMultiProcessor);
  printout("[CUDA] maxThreadsPerBlock %d\n", deviceProp.maxThreadsPerBlock);
  printout("[CUDA] warpSize %d\n", deviceProp.warpSize);
#endif

  globals::kappa_rpkt_cont = (struct rpkt_cont_opacity *)calloc(get_max_threads(), sizeof(struct rpkt_cont_opacity));
  assert_always(globals::kappa_rpkt_cont != NULL);

  /// Using this and the global variable output_file opens and closes the output_file
  /// only once, which speeds up the simulation with a lots of output switched on (debugging).
  /// The downside is that while the simulation runs, its output is only readable on that
  /// machine where the simulation is running.
  /// NB: printout also needs some changes to get this working!
  /*
  output_file = fopen_required("output.txt", "w");
  /// Makes sure that the output_file is written line-by-line
  setvbuf(output_file, NULL, _IOLBF, 1);
  */

  /*
  ldist_file = fopen_required("ldist.out", "w");
  */

  // snprintf(filename, 128, "tb%.4d.txt", my_rank);
  //(tb_file = fopen_required(filename, "w");
  // setvbuf(tb_file, NULL, _IOLBF, 1);

  // printout("CELLHISTORYSIZE %d\n",CELLHISTORYSIZE);

  input(my_rank);

#ifdef RECORD_LINESTAT
  if (my_rank == 0) {
    initialise_linestat_file();
  }
#endif

  printout("time after input %ld\n", time(NULL));
  printout("timesteps %d\n", globals::ntstep);

  /// Precalculate the rate coefficients for spontaneous and stimulated recombination
  /// and for photoionisation. With the nebular approximation they only depend on T_e
  /// T_R and W. W is easily factored out. For stimulated recombination we must assume
  /// T_e = T_R for this precalculation.
  /// Make this parallel ?
  printout("time before tabulation of rate coefficients %ld\n", time(NULL));
  ratecoefficients_init();
  printout("time after tabulation of rate coefficients %ld\n", time(NULL));
  //  abort();
#ifdef MPI_ON
  printout("barrier after tabulation of rate coefficients: time before barrier %ld, ", time(NULL));
  MPI_Barrier(MPI_COMM_WORLD);
  printout("time after barrier %ld\n", time(NULL));
#endif

  stats::init();

  /// Record the chosen syn_dir
  FILE *syn_file = fopen_required("syn_dir.txt", "w");
  fprintf(syn_file, "%g %g %g", globals::syn_dir[0], globals::syn_dir[1], globals::syn_dir[2]);
  fclose(syn_file);
  printout("time read syn file %ld\n", time(NULL));

  bool terminate_early = false;
  globals::file_set = false;  // LJS: deprecate this switch?
  globals::debuglevel = 4;    /// Selects detail level of debug output, needs still some work.

  /// Initialise the grid. Call routine that sets up the initial positions
  /// and sizes of the grid cells.
  time_init();
  if (my_rank == 0) {
    write_timestep_file();
  }
  printout("time grid_init %ld\n", time(NULL));
  grid::grid_init(my_rank);

  printout("Simulation propagates %g packets per process (total %g with nprocs %d)\n", 1. * globals::npkts,
           1. * globals::npkts * globals::nprocs, globals::nprocs);

  printout("[info] mem_usage: packets occupy %.3f MB\n", MPKTS * sizeof(struct packet) / 1024. / 1024.);

  if (!globals::simulation_continued_from_saved) {
    std::remove("deposition.out");
    /// Next we want to initialise the packets.
    /// Create a bunch of npkts packets
    /// and write them to a binary file for later readin.
    packet_init(my_rank, packets);
    zero_estimators();
  }

  /// For the parallelisation of update_grid, the process needs to be told which cells belong to it.
  /// The next loop is over all grid cells. For parallelisation, we want to split this loop between
  /// processes. This is done by assigning each MPI process nblock cells. The residual n_leftover
  /// cells are sent to processes 0 ... process n_leftover -1.
  int nstart = 0;
  int ndo = 0;
  int ndo_nonempty = 0;
  int maxndo = 0;
  grid::get_nstart_ndo(my_rank, globals::nprocs, &nstart, &ndo, &ndo_nonempty, &maxndo);
  printout("process rank %d (global max rank %d) assigned %d modelgrid cells (%d nonempty)", my_rank,
           globals::nprocs - 1, ndo, ndo_nonempty);
  if (ndo > 0) {
    printout(": cells [%d..%d] (model has max mgi %d)\n", nstart, nstart + ndo - 1, grid::get_npts_model() - 1);
  } else {
    printout("\n");
  }

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
  /// Initialise the exchange buffer
  /// The factor 4 comes from the fact that our buffer should contain elements of 4 byte
  /// instead of 1 byte chars. But the MPI routines don't care about the buffers datatype
  mpi_grid_buffer_size = 4 * ((12 + 4 * get_includedions()) * (maxndo) + 1);
  printout("reserve mpi_grid_buffer_size %d space for MPI communication buffer\n", mpi_grid_buffer_size);
  // char buffer[mpi_grid_buffer_size];
  mpi_grid_buffer = (char *)malloc(mpi_grid_buffer_size * sizeof(char));
  if (mpi_grid_buffer == NULL) {
    printout("[fatal] input: not enough memory to initialize MPI grid buffer ... abort.\n");
    abort();
  }
#endif

  /** That's the end of the initialisation. */
  /** *******************************************/

  /// Now comes the loop over timesteps. Before doing this we need to
  /// define the time steps.
  // time_init();

  /// Standard approach: for loop over given number of time steps
  // for (nts = globals::itstep; nts < ftstep; nts++)
  //{

  /// Now use while loop to allow for timed restarts
  const int last_loop = globals::ftstep;
  int nts = globals::itstep;

  macroatom_open_file(my_rank);
  if (ndo > 0) {
    assert_always(estimators_file == NULL);
    snprintf(filename, 128, "estimators_%.4d.out", my_rank);
    estimators_file = fopen_required(filename, "w");

    if (NLTE_POPS_ON && ndo_nonempty > 0) {
      nltepop_open_file(my_rank);
    }
  }

// Initialise virtual packets file and vspecpol
#ifdef VPKT_ON
  init_vspecpol();
  if (vgrid_flag == 1) {
    init_vpkt_grid();
  }

  if (globals::simulation_continued_from_saved) {
    // Continue simulation: read into temporary files

    read_vspecpol(my_rank, nts);

    if (vgrid_flag == 1) {
      if (nts % 2 == 0) {
        snprintf(filename, 128, "vpkt_grid_%d_%d_odd.tmp", 0, my_rank);
      } else {
        snprintf(filename, 128, "vpkt_grid_%d_%d_even.tmp", 0, my_rank);
      }

      FILE *vpktgrid_file = fopen_required(filename, "rb");

      read_vpkt_grid(vpktgrid_file);

      fclose(vpktgrid_file);
    }
  }
#endif

  while (nts < last_loop && !terminate_early) {
    globals::nts_global = nts;
#ifdef MPI_ON
    //        const time_t time_before_barrier = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    //        const time_t time_after_barrier = time(NULL);
    //        printout("timestep %d: time before barrier %d, time after barrier %d\n", nts, time_before_barrier,
    //        time_after_barrier);
#endif

#ifdef DO_TITER
    // The first time step must solve the ionisation balance in LTE
    globals::initial_iteration = (nts == 0);
#else
    /// Do 3 iterations on timestep 0-9
    /*if (nts == 0)
    {
      n_titer = 3;
      initial_iteration = true;
    }
    else if (nts < 6)
    {
      n_titer = 3;
      initial_iteration = false;
    }
    else
    {
      n_titer = 1;
      initial_iteration = false;
    }*/
    globals::n_titer = 1;
    globals::initial_iteration = (nts < globals::n_lte_timesteps);
#endif

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

  if (my_rank == 0) {
    fclose(linestat_file);
  }
  // fclose(ldist_file);
  // fclose(output_file);

  /* Spec syn. */
  // grid_init();
  // syn_gamma();

  if ((globals::ntstep != globals::ftstep) || (terminate_early)) {
    printout("RESTART_NEEDED to continue model\n");
  } else {
    printout("No need for restart\n");
  }

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  printout("simulation finished at %ld\n", time(NULL));
#if CUDA_ENABLED
  cudaDeviceReset();
#endif
  // fclose(tb_file);
  if (estimators_file != NULL) {
    fclose(estimators_file);
  }

  macroatom_close_file();
  if (NLTE_POPS_ON) nltepop_close_file();

  radfield::close_file();
  nonthermal::close_file();

#ifdef _OPENMP
#pragma omp parallel
#endif
  { fclose(output_file); }

#ifdef _OPENMP
  omp_set_dynamic(0);
#pragma omp parallel
#endif
  {
#ifndef __CUDA_ARCH__
    gsl_integration_workspace_free(gslworkspace);
#endif
  }

  free(packets);
  if (TRACK_ION_STATS) {
    stats::cleanup();
  }

  decay::cleanup();

#ifdef MPI_ON
  MPI_Finalize();
#endif

  return 0;
}

// printout should be used instead of printf throughout the whole code for output messages
extern inline void gsl_error_handler_printout(const char *reason, const char *file, int line, int gsl_errno);
extern inline FILE *fopen_required(const char *filename, const char *mode);
