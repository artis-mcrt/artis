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

#include <getopt.h>
#include <unistd.h>
// #include <filesystem>
#include "sn3d.h"
#include "emissivities.h"
#include "grey_emissivities.h"
#include "grid.h"
#include "input.h"
#include "light_curve.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "packet_init.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "stats.h"
#include "update_grid.h"
#include "update_packets.h"
#include "version.h"
#include "vpkt.h"

const bool KEEP_ALL_RESTART_FILES = false; // once a new gridsave and packets*.tmp have been written, don't delete the previous set

// threadprivate variables
int tid;
__managed__ int myGpuId = 0;
__managed__ bool use_cellhist;
__managed__ bool neutral_flag;
#ifndef __CUDA_ARCH__
gsl_rng *rng;
#else
__device__ void *rng = NULL;
#endif
gsl_integration_workspace *gslworkspace = NULL;
FILE *output_file;
static FILE *linestat_file;
static time_t real_time_start;
static time_t time_timestep_start = -1; // this will be set after the first update of the grid and before packet prop
static FILE *estimators_file;
static int nstart = 0;
static int ndo = 0;

int mpi_grid_buffer_size = 0;
char *mpi_grid_buffer = NULL;


static void initialise_linestat_file(void)
{
  linestat_file = fopen_required("linestat.out", "w");

  for (int i = 0; i < globals::nlines; i++)
    fprintf(linestat_file, "%g ", CLIGHT / globals::linelist[i].nu);
  fprintf(linestat_file,"\n");

  for (int i = 0; i < globals::nlines; i++)
    fprintf(linestat_file, "%d ", get_element(globals::linelist[i].elementindex));
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++)
    fprintf(linestat_file, "%d ", get_ionstage(globals::linelist[i].elementindex, globals::linelist[i].ionindex));
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++)
    fprintf(linestat_file, "%d ", globals::linelist[i].upperlevelindex + 1);
  fprintf(linestat_file, "\n");

  for (int i = 0; i < globals::nlines; i++)
    fprintf(linestat_file, "%d ", globals::linelist[i].lowerlevelindex + 1);
  fprintf(linestat_file, "\n");

  fflush(linestat_file);
  //setvbuf(linestat_file, NULL, _IOLBF, 1); // flush after every line makes it slow!
}


#ifdef MPI_ON
static void mpi_communicate_grid_properties(const int my_rank, const int p, const int nstart, const int ndo, const int nts, const int titer, char *mpi_grid_buffer, const int mpi_grid_buffer_size)
{
  int position = 0;
  for (int root = 0; root < p; root++)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    int root_nstart = nstart;
    MPI_Bcast(&root_nstart, 1, MPI_INT, root, MPI_COMM_WORLD);
    int root_ndo = ndo;
    MPI_Bcast(&root_ndo, 1, MPI_INT, root, MPI_COMM_WORLD);

    for (int modelgridindex = root_nstart; modelgridindex < (root_nstart + root_ndo); modelgridindex++)
    {
      radfield::do_MPI_Bcast(modelgridindex, root);

      if (get_numassociatedcells(modelgridindex) > 0)
      {
        nonthermal::nt_MPI_Bcast(modelgridindex, root);
        if (NLTE_POPS_ON)
        {
          MPI_Bcast(globals::modelgrid[modelgridindex].nlte_pops, globals::total_nlte_levels, MPI_DOUBLE, root, MPI_COMM_WORLD);
        }
      }
    }

    if (root == my_rank)
    {
      position = 0;
      MPI_Pack(&ndo, 1, MPI_INT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
      for (int mgi = nstart; mgi < (nstart + ndo); mgi++)
      //for (int nncl = 0; nncl < ndo; nncl++)
      {
        //nn = nonemptycells[my_rank+nncl*nprocs];
        MPI_Pack(&mgi, 1, MPI_INT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
        //if (globals::cell[nn].rho > MINDENSITY)
        if (get_numassociatedcells(mgi) > 0)
        {
          MPI_Pack(&globals::modelgrid[mgi].Te, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&globals::modelgrid[mgi].TR, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&globals::modelgrid[mgi].TJ, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&globals::modelgrid[mgi].W, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&globals::modelgrid[mgi].rho, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&globals::modelgrid[mgi].nne, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&globals::modelgrid[mgi].nnetot, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&globals::modelgrid[mgi].kappagrey, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&globals::modelgrid[mgi].totalcooling, 1, MPI_DOUBLE, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&globals::modelgrid[mgi].thick, 1, MPI_SHORT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);

          for (int element = 0; element < get_nelements(); element++)
          {
            MPI_Pack(globals::modelgrid[mgi].composition[element].groundlevelpop, get_nions(element), MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
            MPI_Pack(globals::modelgrid[mgi].composition[element].partfunct, get_nions(element), MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
            MPI_Pack(globals::modelgrid[mgi].cooling[element].contrib, get_nions(element), MPI_DOUBLE, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          }
        }
      }
      printout("mem_usage: MPI_BUFFER: used %d of %d bytes allocated to mpi_grid_buffer\n", position, mpi_grid_buffer_size);
      assert(position <= mpi_grid_buffer_size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(mpi_grid_buffer, mpi_grid_buffer_size, MPI_PACKED, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    position = 0;
    int nlp;
    MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &nlp, 1, MPI_INT, MPI_COMM_WORLD);
    for (int nn = 0; nn < nlp; nn++)
    {
      int mgi;
      MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &mgi, 1, MPI_INT, MPI_COMM_WORLD);
      //if (globals::cell[ncl].rho > MINDENSITY)
      if (get_numassociatedcells(mgi) > 0)
      {
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &globals::modelgrid[mgi].Te, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &globals::modelgrid[mgi].TR, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &globals::modelgrid[mgi].TJ, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &globals::modelgrid[mgi].W, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &globals::modelgrid[mgi].rho, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &globals::modelgrid[mgi].nne, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &globals::modelgrid[mgi].nnetot, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &globals::modelgrid[mgi].kappagrey, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &globals::modelgrid[mgi].totalcooling, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &globals::modelgrid[mgi].thick, 1, MPI_SHORT, MPI_COMM_WORLD);

        for (int element = 0; element < get_nelements(); element++)
        {
          MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, globals::modelgrid[mgi].composition[element].groundlevelpop, get_nions(element), MPI_FLOAT, MPI_COMM_WORLD);
          MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, globals::modelgrid[mgi].composition[element].partfunct, get_nions(element), MPI_FLOAT, MPI_COMM_WORLD);
          MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, globals::modelgrid[mgi].cooling[element].contrib, get_nions(element), MPI_DOUBLE, MPI_COMM_WORLD);
        }
      }
    }
  }

  #ifndef FORCE_LTE
    #if (!NO_LUT_PHOTOION)
      if ((!globals::simulation_continued_from_saved) || (nts - globals::itstep != 0) || (titer != 0))
      {
        MPI_Barrier(MPI_COMM_WORLD);
        /// Reduce the corrphotoionrenorm array.
        printout("nts %d, titer %d: bcast corr photoionrenorm\n", nts, titer);
        MPI_Allreduce(MPI_IN_PLACE, &corrphotoionrenorm, MMODELGRID * get_nelements() * get_max_nions(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /// Reduce the gammaestimator array. Only needed to write restart data.
        printout("nts %d, titer %d: bcast gammaestimator\n", nts, titer);
        MPI_Allreduce(MPI_IN_PLACE, &gammaestimator, MMODELGRID * get_nelements() * get_max_nions(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
    #endif
  #endif
  MPI_Barrier(MPI_COMM_WORLD);
}


static void mpi_reduce_estimators(int my_rank, int nts)
{
  radfield::reduce_estimators();
  #ifndef FORCE_LTE
    MPI_Allreduce(MPI_IN_PLACE, &globals::ffheatingestimator, MMODELGRID, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &globals::colheatingestimator, MMODELGRID, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #if (!NO_LUT_PHOTOION)
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &gammaestimator, MMODELGRID * get_nelements() * get_max_nions(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif
    #if (!NO_LUT_BFHEATING)
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &bfheatingestimator, MMODELGRID * get_nelements() * get_max_nions(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif
  #endif

  #ifdef RECORD_LINESTAT
    MPI_Allreduce(MPI_IN_PLACE, globals::ecounter, globals::nlines, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, globals::acounter, globals::nlines, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #endif

  //double deltaV = pow(wid_init * globals::time_step[nts].mid/globals::tmin, 3.0);
  //double deltat = globals::time_step[nts].width;
  if (globals::do_rlc_est != 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, &globals::rpkt_emiss, MMODELGRID, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (globals::do_comp_est)
  {
    MPI_Allreduce(MPI_IN_PLACE, &globals::compton_emiss, MMODELGRID * EMISS_MAX, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  }

  /// Communicate gamma and positron deposition and write to file
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].cmf_lum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].gamma_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &globals::time_step[nts].positron_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  globals::time_step[nts].cmf_lum /= globals::nprocs;
  globals::time_step[nts].gamma_dep /= globals::nprocs;
  globals::time_step[nts].positron_dep /= globals::nprocs;

  #if TRACK_ION_STATS
  MPI_Allreduce(MPI_IN_PLACE, ionstats, npts_model * includedions * ION_STAT_COUNT, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif

  MPI_Barrier(MPI_COMM_WORLD);
}
#endif



static void write_temp_packetsfile(const int timestep, const int my_rank, const PKT *const pkt)
{
  char filename[100];
  sprintf(filename, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  printout("Writing %s...", filename);
  FILE *packets_file = fopen_required(filename, "wb");

  fwrite(pkt, sizeof(PKT), globals::npkts, packets_file);
  fclose(packets_file);
  printout("done\n");
}


static void remove_temp_packetsfile(const int timestep, const int my_rank)
{
  char filename[100];
  sprintf(filename, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  if (!access(filename, F_OK))
  // if (std::filesystem::exists(filename))
  {
    remove(filename);
    printout("Deleted %s\n", filename);
  }
}


static void remove_grid_restart_data(const int timestep)
{
  char prevfilename[100];
  sprintf(prevfilename, "gridsave_ts%d.tmp", timestep);

  if (!access(prevfilename, F_OK))
  // if (std::filesystem::exists(prevfilename))
  {
    remove(prevfilename);
    printout("Deleted %s\n", prevfilename);
  }
}


static bool walltime_sufficient_to_continue(const int nts, const int nts_prev, const int walltimelimitseconds)
{
  #ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // time is measured from just before packet propagation from one timestep to the next
  const int estimated_time_per_timestep = time(NULL) - time_timestep_start;
  printout("TIME: time between timesteps is %d seconds (measured packet prop of ts %d and update grid of ts %d)\n", estimated_time_per_timestep, nts_prev, nts);

  bool do_this_full_loop = true;
  if (walltimelimitseconds > 0)
  {
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
      printout("TIMED_RESTARTS: Going to continue since remaining time %d s >= 1.5 * time_per_timestep\n", wallclock_remaining_seconds);
    else
      printout("TIMED_RESTARTS: Going to terminate since remaining time %d s < 1.5 * time_per_timestep\n", wallclock_remaining_seconds);
  }
  return do_this_full_loop;
}


#if CUDA_ENABLED
void* makemanaged(void* ptr, size_t curSize)
{
  if (ptr == NULL)
  {
    return NULL;
  }
  void *newptr;
  cudaMallocManaged(&newptr, curSize);
  memcpy(newptr, ptr, curSize);
  free(ptr);
  return newptr;
}
#endif


static void save_grid_and_packets(
  const int nts, const int my_rank, PKT* packets)
{
  const time_t time_write_packets_file_start = time(NULL);
  printout("time before write temporary packets file %ld\n", time_write_packets_file_start);

  // save packet state at start of current timestep (before propagation)
  write_temp_packetsfile(nts, my_rank, packets);

  #ifdef VPKT_ON
  char filename[100];
  if (nts % 2 == 0)
    sprintf(filename,"vspecpol_%d_%d_even.tmp", 0, my_rank);
  else
    sprintf(filename,"vspecpol_%d_%d_odd.tmp", 0, my_rank);

  FILE *vspecpol_file = fopen_required(filename, "wb");

  write_vspecpol(vspecpol_file);
  fclose(vspecpol_file);

  // Write temporary files for vpkt_grid
  if (vgrid_flag == 1)
  {
    if (nts % 2 == 0)
      sprintf(filename,"vpkt_grid_%d_%d_even.tmp", 0, my_rank);
    else
      sprintf(filename,"vpkt_grid_%d_%d_odd.tmp", 0, my_rank);

    FILE *vpkt_grid_file = fopen_required(filename, "wb");

    write_vpkt_grid(vpkt_grid_file);
    fclose(vpkt_grid_file);
  }
  #endif

  printout("time after write temporary packets file %ld (took %ld seconds)\n", time(NULL), time(NULL) - time_write_packets_file_start);

  if (my_rank == 0)
  {
    write_grid_restart_data(nts);
    update_parameterfile(nts);
  }

  if (!KEEP_ALL_RESTART_FILES)
  {
    // ensure new packets files have been written by all processes before we remove the old set
    #ifdef MPI_ON
      MPI_Barrier(MPI_COMM_WORLD);
    #endif

    if (my_rank == 0)
      remove_grid_restart_data(nts - 1);

    // delete temp packets files from previous timestep now that all restart data for the new timestep is available
    remove_temp_packetsfile(nts - 1, my_rank);
  }
}

static bool do_timestep(
  const int outer_iteration, const int nts, const int titer,
  const int my_rank, PKT* packets, const int walltimelimitseconds)
{
  bool do_this_full_loop = true;

  const int nts_prev = (titer != 0 || nts == 0) ? nts : nts - 1;
  if ((titer > 0) || (globals::simulation_continued_from_saved && (nts == globals::itstep)))
  {
    /// Read the packets file to reset before each additional iteration on the timestep
    read_temp_packetsfile(nts, my_rank, packets);
  }

  /// Some counters on pkt-actions need to be reset to do statistics
  stats::pkt_action_counters_reset();

  if (nts == 0)
  {
    radfield::initialise_prev_titer_photoionestimators();
  }

  #ifdef RECORD_LINESTAT
    // The same for absorption/emission of r-pkts in lines
    for (int i = 0; i < globals::nlines; i++)
    {
      globals::acounter[i] = 0;
      globals::ecounter[i] = 0;
    }
  #endif

  globals::do_comp_est = globals::do_r_lc ? false : estim_switch(nts, globals::time_step);

  // Update the matter quantities in the grid for the new timestep.

  const time_t sys_time_start_update_grid = time(NULL);
  printout("\ntimestep %d: time before update grid %ld (tstart + %ld)\n",
           nts, sys_time_start_update_grid, sys_time_start_update_grid - real_time_start);

  #ifndef FORCE_LTE
    #if (!NO_LUT_PHOTOION)
      /// Initialise globals::corrphotoionrenorm[i] to zero before update_grid is called
      /// This allows reduction after update_grid has finished
      /// unless they have been read from file and must neither be touched
      /// nor broadcasted after update_grid
      if ((!globals::simulation_continued_from_saved) || (nts - globals::itstep != 0) || (titer != 0))
      {
        printout("nts %d, titer %d: reset corr photoionrenorm\n",nts,titer);
        for (int i = 0; i < MMODELGRID * get_nelements() * get_max_nions(); i++)
        {
          globals::corrphotoionrenorm[i] = 0.;
        }
        printout("after nts %d, titer %d: reset corr photoionrenorm\n",nts,titer);
      }
    #endif
  #endif

  update_grid(estimators_file, nts, nts_prev, my_rank, nstart, ndo, titer);

  const time_t sys_time_finish_update_grid = time(NULL);
  printout("timestep %d: update_grid: process %d finished update grid at %ld (took %ld seconds)\n",
           nts, my_rank, sys_time_finish_update_grid, sys_time_finish_update_grid - sys_time_start_update_grid);

  #ifdef DO_TITER
    /// No iterations over the zeroth timestep, set titer > n_titer
    if (nts == 0)
      titer = n_titer + 1;
  #endif
  #ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
  #endif
  const time_t sys_time_finish_update_grid_all_processes = time(NULL);
  printout("timestep %d: waiting for update grid to finish on other processes took %ld seconds\n",
           nts, sys_time_finish_update_grid_all_processes - sys_time_finish_update_grid);
  printout("timestep %d: time after update grid for all processes %ld (took %ld seconds)\n",
           nts, sys_time_finish_update_grid_all_processes, sys_time_finish_update_grid_all_processes - sys_time_start_update_grid);

  const time_t sys_time_start_communicate_grid = time(NULL);

  /// Each process has now updated its own set of cells. The results now need to be communicated between processes.
  #ifdef MPI_ON
    mpi_communicate_grid_properties(my_rank, globals::nprocs, nstart, ndo, nts, titer, mpi_grid_buffer, mpi_grid_buffer_size);
  #endif

  printout("timestep %d: time after grid properties have been communicated %ld (took %ld seconds)\n",
           nts, time(NULL), time(NULL) - sys_time_start_communicate_grid);

  /// If this is not the 0th time step of the current job step,
  /// write out a snapshot of the grid properties for further restarts
  /// and update input.txt accordingly
  if (((nts - globals::itstep) != 0))
  {
    save_grid_and_packets(nts, my_rank, packets);
    do_this_full_loop = walltime_sufficient_to_continue(nts, nts_prev, walltimelimitseconds);
  }
  time_timestep_start = time(NULL);

  // set all the estimators to zero before moving packets. This is now done
  // after update_grid so that, if requires, the gamma-ray heating estimator is known there
  // and also the photoion and stimrecomb estimators
  zero_estimators();

  // MPI_Barrier(MPI_COMM_WORLD);
  if ((nts < globals::ftstep) && do_this_full_loop)
  {
    /// Now process the packets.
    const time_t time_update_packets_start = time(NULL);
    printout("timestep %d: time before update packets %ld\n", nts, time_update_packets_start);

    update_packets(my_rank, nts, packets);

    stats::pkt_action_counters_printout(packets, nts);

    #ifdef MPI_ON
      MPI_Barrier(MPI_COMM_WORLD); // hold all processes once the packets are updated
      const time_t time_communicate_estimators_start = time(NULL);
    #endif
    printout("timestep %d: time after update packets %ld (took %ld seconds)\n", nts, time(NULL), time(NULL) - time_update_packets_start);

    #ifdef MPI_ON
      // All the processes have their own versions of the estimators for this time step now.
      // Since these are going to be needed in the next time step, we will gather all the
      // estimators together now, sum them, and distribute the results

      mpi_reduce_estimators(my_rank, nts);
    #endif

    // The estimators have been summed across all proceses and distributed.
    // They will now be normalised independently on all processes
    if (globals::do_comp_est)
    {
      normalise_compton_estimators(nts, globals::time_step);
      if (my_rank == 0)
      {
        write_compton_estimators(nts);
      }
    }

    if (globals::do_rlc_est != 0)
    {
      normalise_grey(nts);
      if ((globals::do_rlc_est != 3) && (my_rank == 0))
      {
        write_grey(nts);
      }
    }

    if (my_rank == 0)
    {
      FILE *dep_file = fopen_required("deposition.out", "w");
      for (int i = 0; i <= nts; i++)
      {
        fprintf(dep_file, "%g %g %g %g\n", globals::time_step[i].mid / DAY,
                globals::time_step[i].gamma_dep / globals::time_step[i].width / LSUN,
                globals::time_step[i].positron_dep / globals::time_step[i].width / LSUN,
                (globals::time_step[i].gamma_dep + globals::time_step[i].positron_dep) / globals::time_step[i].width / LSUN);
      }
      fclose(dep_file);
    }
    write_partial_lightcurve(my_rank, nts, packets);

    #ifdef MPI_ON
      printout("timestep %d: time after estimators have been communicated %ld (took %ld seconds)\n", nts, time(NULL), time(NULL) - time_communicate_estimators_start);
    #endif

    printout("%d: During timestep %d on MPI process %d, %d pellets decayed and %d packets escaped. (t=%gd)\n",
             outer_iteration, nts, my_rank, globals::time_step[nts].pellet_decays, globals::nesc, globals::time_step[nts].mid / DAY);

    #ifdef VPKT_ON
      printout("%d: During timestep %d on MPI process %d, %d virtual packets were generated and %d escaped. \n",
               outer_iteration, nts, my_rank, nvpkt, nvpkt_esc1 + nvpkt_esc2 + nvpkt_esc3);
      printout("%d virtual packets came from an electron scattering event, %d from a kpkt deactivation and %d from a macroatom deactivation. \n",
               nvpkt_esc1, nvpkt_esc2, nvpkt_esc3);

      nvpkt = 0;
      nvpkt_esc1 = 0;
      nvpkt_esc2 = 0;
      nvpkt_esc3 = 0;
    #endif


    #ifdef RECORD_LINESTAT
      if (my_rank == 0)
      {
        /// Print net absorption/emission in lines to the linestat_file
        /// Currently linestat information is only properly implemented for MPI only runs
        /// For hybrid runs only data from thread 0 is recorded
        for (int i = 0; i < globals::nlines; i++)
          fprintf(linestat_file, "%d ", globals::ecounter[i]);
        fprintf(linestat_file, "\n");
        for (int i = 0; i < globals::nlines; i++)
          fprintf(linestat_file, "%d ", globals::acounter[i]);
        fprintf(linestat_file, "\n");
        fflush(linestat_file);
      }
    #endif

    if (nts == globals::ftstep - 1)
    {
      char filename[100];
      sprintf(filename, "packets%.2d_%.4d.out", 0, my_rank);
      //sprintf(filename,"packets%.2d_%.4d.out", middle_iteration, my_rank);
      write_packets(filename, packets);

      // write specpol of the virtual packets
      #ifdef VPKT_ON
        sprintf(filename,"vspecpol_%d-%d.out",my_rank,tid);
        FILE *vspecpol_file = fopen_required(filename, "w");

        write_vspecpol(vspecpol_file);
        fclose(vspecpol_file);

        if (vgrid_flag == 1)
        {
          sprintf(filename,"vpkt_grid_%d-%d.out",my_rank,tid);
          FILE *vpkt_grid_file = fopen_required(filename, "w");
          write_vpkt_grid(vpkt_grid_file);
          fclose(vpkt_grid_file);
        }
      #endif

      printout("time after write final packets file %ld\n", time(NULL));

      // final packets*.out have been written, so remove the temporary packets files
      // commented out because you might still want to resume the simulation
      // sprintf(filename, "packets%d_%d_odd.tmp", 0, my_rank);
      // remove(filename);
      // sprintf(filename, "packets%d_%d_even.tmp", 0, my_rank);
      // remove(filename);
    }
  }
  return !do_this_full_loop;
}


static void get_nstart_ndo(int my_rank, int nprocesses, int *nstart, int *ndo, int *maxndo)
{
  #ifndef MPI_ON
    // no MPI, single process updates all cells
    *nstart = 0;
    *ndo = globals::npts_model;
    return;
  #endif

  int n_leftover = 0;

  int nblock = globals::npts_model / nprocesses; // integer division, minimum cells for any process
  const int numtot = nblock * nprocesses; // cells counted if all processes do the minimum number of cells
  if (numtot > globals::npts_model) // LJS: should never be the case?
  {
    nblock = nblock - 1;
    *maxndo = nblock + 1;
    n_leftover = globals::npts_model - (nblock * nprocesses);
  }
  else if (numtot < globals::npts_model)
  {
    *maxndo = nblock + 1;
    n_leftover = globals::npts_model - (nblock * nprocesses);
  }
  else
  {
    *maxndo = nblock;
    n_leftover = 0;
  }

  if (my_rank < n_leftover)
  {
    *ndo = nblock + 1;
    *nstart = my_rank * (nblock + 1);
  }
  else
  {
    *ndo = nblock;
    *nstart = n_leftover * (nblock + 1) + (my_rank - n_leftover) * (nblock);
  }
}


int main(int argc, char** argv)
// Main - top level routine.
{
  //FILE *temperature_file;
  char filename[100];

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
  assert(deviceCount > 0);
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
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
  #else
    int my_rank = 0;
    int p = 1;
  #endif

  globals::nprocs = p;              /// Global variable which holds the number of MPI processes
  globals::rank_global = my_rank;   /// Global variable which holds the rank of the active MPI process

#ifdef _OPENMP
  /// Explicitly turn off dynamic threads because we use the threadprivate directive!!!
  omp_set_dynamic(0);

  #pragma omp parallel private(filename)
#endif
  {
    /// Get the current threads ID, copy it to a threadprivate variable
    tid = get_thread_num();
    /// and initialise the threads outputfile
    sprintf(filename,"output_%d-%d.txt", my_rank, tid);
    output_file = fopen_required(filename, "w");
    /// Makes sure that the output_file is written line-by-line
    setvbuf(output_file, NULL, _IOLBF, 1);

#   ifdef _OPENMP
    printout("OpenMP parallelisation active with %d threads (max %d)\n", get_num_threads(), get_max_threads());
#   endif

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
    switch (opt)
    {
      case 'w':
      {
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

      default:
      {
        fprintf(stderr, "Usage: %s [-w WALLTIMELIMITHOURS]\n", argv[0]);
        exit(EXIT_FAILURE);
      }
    }
  }

  #if CUDA_ENABLED
  PKT *packets;
  cudaMallocManaged(&packets, MPKTS * sizeof(PKT));
    #if USECUDA_UPDATEPACKETS
    cudaMemAdvise(packets, MPKTS * sizeof(PKT), cudaMemAdviseSetPreferredLocation, myGpuId);
    #endif
  #else
  PKT *const packets = (PKT *) malloc(MPKTS * sizeof(PKT));
  #endif

  assert(packets != NULL);

  #ifndef GIT_BRANCH
    #define GIT_BRANCH "UNKNOWN"
  #endif
  printout("ARTIS git branch %s\n", GIT_BRANCH);

  #ifndef GIT_VERSION
    #define GIT_VERSION "UNKNOWN"
  #endif
  printout("Current version: %s\n", GIT_VERSION);

  //printout("Hash of most recent commit: %s\n",GIT_HASH);
  printout("Compiled at %s on %s\n", __TIME__, __DATE__);

  #ifdef MPI_ON
    printout("MPI enabled with %d processes\n", globals::nprocs);
  #else
    printout("MPI disabled\n");
  #endif

  #ifdef __CUDACC__
  printout("[CUDA] NVIDIA CUDA is available\n");
  #else
  printout("[CUDA] NVIDIA CUDA is not available\n");
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
  #else
  printout("[CUDA] NVIDIA CUDA accelerated routines are disabled\n");
  #endif

  if ((globals::kappa_rpkt_cont = (rpkt_cont_opacity_struct *) calloc(get_max_threads(), sizeof(rpkt_cont_opacity_struct))) == NULL)
  {
    printout("[fatal] input: error initializing continuum opacity communication variables ... abort\n");
    abort();
  }

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


  //sprintf(filename,"tb%.4d.txt",my_rank);
  //(tb_file = fopen_required(filename, "w");
  //setvbuf(tb_file, NULL, _IOLBF, 1);

  macroatom_open_file(my_rank);

  sprintf(filename, "estimators_%.4d.out", my_rank);
  estimators_file = fopen_required(filename, "w");
  //setvbuf(estimators_file, NULL, _IOLBF, 1);

  if (NLTE_POPS_ON)
    nltepop_open_file(my_rank);

  //printout("CELLHISTORYSIZE %d\n",CELLHISTORYSIZE);

  /// Get input stuff
  input(my_rank);

  #ifdef RECORD_LINESTAT
  if (my_rank == 0)
  {
    initialise_linestat_file();
  }
  #endif

  printout("time after input %ld\n", time(NULL));
  printout("simulation propagates %d packets per process through a %d x %d x %d grid\n",
           globals::npkts, globals::ncoordgrid[0], globals::ncoordgrid[1], globals::ncoordgrid[2]);
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
//  #ifdef MPI_ON
//    const time_t time_before_barrier = time(NULL);
//    MPI_Barrier(MPI_COMM_WORLD);
//    const time_t time_after_barrier = time(NULL);
//    printout("barrier after tabulation of rate coefficients: time before barrier %d, time after barrier %d\n", time_before_barrier, time_after_barrier);
//  #endif

  stats::init();

  /// As a precaution, explicitly zero all the estimators here
  zero_estimators();

  printout("time after zero estimators %ld\n", time(NULL));

  /// Record the chosen syn_dir
  FILE *syn_file = fopen_required("syn_dir.txt", "w");
  fprintf(syn_file, "%g %g %g", globals::syn_dir[0], globals::syn_dir[1], globals::syn_dir[2]);
  fclose(syn_file);
  printout("time read syn file %ld\n", time(NULL));

  bool terminate_early = false;
  globals::file_set = false; // LJS: deprecate this switch?
  globals::debuglevel = 4;  /// Selects detail level of debug output, needs still some work.
  for (int outer_iteration = 0; outer_iteration < globals::n_out_it; outer_iteration++)
  {
    /// Initialise the grid. Call routine that sets up the initial positions
    /// and sizes of the grid cells.
    time_init();
    if (my_rank == 0)
    {
      write_timestep_file();
    }
    printout("time time init %ld\n", time(NULL));
    grid_init(my_rank);

    printout("mem_usage: packets occupy %.1f MB\n", MPKTS * sizeof(PKT) / 1024. / 1024.);

    if (!globals::simulation_continued_from_saved)
    {
      /// Next we want to initialise the packets.
      /// To overcome memory limitations for large numbers of packets, which need to be
      /// propagated on the _same_ grid, this middle_iteration loop was introduced.
      for (int middle_iteration = 0; middle_iteration < globals::n_middle_it; middle_iteration++)
      {
        /// Create a bunch of npkts packets
        /// and write them to a binary file for later readin.
        packet_init(middle_iteration, my_rank, packets);
      }
    }

    /// For the parallelisation of update_grid, the process needs to be told which cells belong to it.
    /// The next loop is over all grid cells. For parallelisation, we want to split this loop between
    /// processes. This is done by assigning each MPI process nblock cells. The residual n_leftover
    /// cells are sent to processes 0 ... process n_leftover -1.
    int maxndo = 0;
    get_nstart_ndo(my_rank, p, &nstart, &ndo, &maxndo);
    printout("process rank %d (of %d) doing %d cells", my_rank, globals::nprocs, ndo);
    if (ndo > 0)
    {
      printout(": numbers %d to %d\n", nstart, nstart + ndo - 1);
    }
    else
    {
      printout("\n");
    }

    #ifdef MPI_ON
      /// Initialise the exchange buffer
      /// The factor 4 comes from the fact that our buffer should contain elements of 4 byte
      /// instead of 1 byte chars. But the MPI routines don't care about the buffers datatype
      mpi_grid_buffer_size = 4 * ((12 + 4 * globals::includedions) * (maxndo) + 1);
      printout("reserve mpi_grid_buffer_size %d space for MPI communication buffer\n", mpi_grid_buffer_size);
      //char buffer[mpi_grid_buffer_size];
      mpi_grid_buffer  = (char *) malloc(mpi_grid_buffer_size * sizeof(char));
      if (mpi_grid_buffer == NULL)
      {
        printout("[fatal] input: not enough memory to initialize MPI grid buffer ... abort.\n");
        abort();
      }
    #endif


    /** That's the end of the initialisation. */
    /** *******************************************/

    /// Now comes the loop over timesteps. Before doing this we need to
    /// define the time steps.
    //time_init();

    /// Standard approach: for loop over given number of time steps
    //for (nts = globals::itstep; nts < ftstep; nts++)
    //{

    /// Now use while loop to allow for timed restarts
    const int last_loop = globals::ftstep;
    int nts = globals::itstep;

    // Initialise virtual packets file and vspecpol
    #ifdef VPKT_ON
    init_vspecpol();
    if (vgrid_flag == 1)
      init_vpkt_grid();
    if (globals::simulation_continued_from_saved)
    {
      // Continue simulation: read into temporary files

      if (nts % 2 == 0)
        sprintf(filename,"vspecpol_%d_%d_odd.tmp", 0, my_rank);
      else
        sprintf(filename,"vspecpol_%d_%d_even.tmp", 0 ,my_rank);

      FILE *vspecpol_file = fopen_required(filename, "rb");

      read_vspecpol(vspecpol_file);

      fclose(vspecpol_file);

      if (vgrid_flag == 1)
      {
        if (nts % 2 == 0)
        {
          sprintf(filename,"vpkt_grid_%d_%d_odd.tmp", 0, my_rank);
        }
        else
        {
          sprintf(filename,"vpkt_grid_%d_%d_even.tmp", 0, my_rank);
        }

        FILE *vpktgrid_file = fopen_required(filename, "rb");

        read_vpkt_grid(vpktgrid_file);

        fclose(vpktgrid_file);
      }
    }
    #endif

    while (nts < last_loop && !terminate_early)
    {
      globals::nts_global = nts;
      #ifdef MPI_ON
//        const time_t time_before_barrier = time(NULL);
        MPI_Barrier(MPI_COMM_WORLD);
//        const time_t time_after_barrier = time(NULL);
//        printout("timestep %d: time before barrier %d, time after barrier %d\n", nts, time_before_barrier, time_after_barrier);
      #endif

      #ifdef DO_TITER
        // The first time step must solve the ionisation balance in LTE
        initial_iteration = (nts == 0);
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

      for (int titer = 0; titer < globals::n_titer; titer++)
      {
        terminate_early = do_timestep(outer_iteration, nts, titer, my_rank, packets, walltimelimitseconds);
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

    if (terminate_early)
    {
      break;
    }
  }

  if (my_rank == 0)
  {
    fclose(linestat_file);
  }
  //fclose(ldist_file);
  //fclose(output_file);

  /* Spec syn. */
  //grid_init();
  //syn_gamma();

  if ((globals::ntstep != globals::ftstep) || (terminate_early))
  {
    printout("RESTART_NEEDED to continue model\n");
  }
  else
  {
    printout("No need for restart\n");
  }

  #ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
  #endif

  printout("simulation finished at %ld\n", time(NULL));
  #if CUDA_ENABLED
  cudaDeviceReset();
  #endif
  //fclose(tb_file);
  fclose(estimators_file);
  macroatom_close_file();
  if (NLTE_POPS_ON)
    nltepop_close_file();
  radfield::close_file();
  nonthermal::close_file();

  #ifdef _OPENMP
    #pragma omp parallel
  #endif
    {
      fclose(output_file);
    }

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
  if (TRACK_ION_STATS)
  {
    stats::cleanup();
  }

  #ifdef MPI_ON
    MPI_Finalize();
  #endif

  return 0;
}


// printout should be used instead of printf throughout the whole code for output messages
extern inline void gsl_error_handler_printout(const char *reason, const char *file, int line, int gsl_errno);
extern inline FILE *fopen_required(const char *filename, const char *mode);


/*void print_opticaldepth(int cellnumber, int timestep, int samplecell, int element)
{
  int nions,nlevels,nuptrans;
  int ion,lower,upper,lineindex;
  int i;
  double epsilon_lower,nu_trans;
  double n_u,n_l,tau_line,A_ul,B_ul,B_lu;
  FILE *tau_file;

  double t_current = globals::time_step[timestep].mid;
  float T_e = globals::cell[cellnumber].T_e;
  double T_R = globals::cell[cellnumber].T_R;
  double W = globals::cell[cellnumber].W;



  sprintf(filename,"tau%.2d_sc%.2d.out",timestep,samplecell);
  tau_file = fopen_required(filename, "w");
  setvbuf(tau_file, NULL, _IOLBF, 1);
  fprintf(tau_file,"%d %g %g %g\n",samplecell,T_e,T_R,W);

  nions = get_nions(element);
  for (ion = 0; ion < nions; ion++)
  {
    nlevels = get_nlevels(element,ion);
    for (lower = 0; lower < nlevels; lower++)
    {
      epsilon_lower = epsilon(element,ion,lower);
      nuptrans = get_nuptrans(element, ion, lower);
      int i = 0; i < nuptrans; i++)
      {
        lineindex = globals::elements[element].ions[ion].levels[lower].uptrans_lineindicies[i];
        upper = linelist[lineindex].upperlevelindex;
        nu_trans = (elements[element].ions[ion].levels[lower].uptrans[i].epsilon - epsilon_lower)/H;

        A_ul = einstein_spontaneous_emission(lineindex);
        B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans,3) * A_ul;
        B_lu = stat_weight(element,ion,upper)/stat_weight(element,ion,lower) * B_ul;

        n_u = calculate_exclevelpop(cellnumber,element,ion,upper);
        n_l = calculate_exclevelpop(cellnumber,element,ion,lower);
        tau_line = (B_lu*n_l - B_ul*n_u) * HCLIGHTOVERFOURPI * t_current;

        fprintf(tau_file,"%g %g %d\n",nu_trans,tau_line,ion);
      }
    }
  }
  fclose(tau_file);
}*/



/*void *my_malloc(size_t size)
{
  //char *adr;
  void *adr;
  #ifdef POWER6
    adr = &heap[heapcounter];
    heapcounter += size;
    if (heapcounter >= heapsize) adr = NULL;
  #else
//    adr = &heap[heapcounter];
//    heapcounter += size;
//    if (heapcounter >= heapsize) adr = NULL;
    adr = malloc(size);
  #endif
  return adr;
}*/
