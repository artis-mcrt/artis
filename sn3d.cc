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

#include <unistd.h>
#include <stdbool.h>
#include <getopt.h>
#include "sn3d.h"
#include "emissivities.h"
#include "grey_emissivities.h"
#include "grid_init.h"
#include "input.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "packet_init.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "update_grid.h"
#include "update_packets.h"
#include "version.h"
#include "vpkt.h"

#if CUDA_ENABLED
#include <cuda_runtime.h>
#endif

const bool KEEP_ALL_RESTART_FILES = false; // once a new gridsave and packets*.tmp have been written, don't delete the previous set

// threadprivate variables
int tid;
bool use_cellhist;
bool neutral_flag;
gsl_rng *rng;
gsl_integration_workspace *gslworkspace;
FILE *output_file;
static FILE *linestat_file;
static time_t real_time_start;
static time_t time_timestep_start = -1; // this will be set after the first update of the grid and before packet prop
static FILE *estimators_file;
static int nstart = 0;
static int ndo = 0;

int mpi_grid_buffer_size = 0;
char *mpi_grid_buffer = NULL;


#ifndef _OPENMP
typedef int omp_int_t;
static inline omp_int_t omp_get_thread_num(void) { return 0; }
static inline omp_int_t omp_get_num_threads(void) { return 1; }
#endif

#if TRACK_ION_STATS
static double *ionstats = NULL;
#endif

static void initialise_linestat_file(void)
{
  linestat_file = fopen_required("linestat.out", "w");

  for (int i = 0; i < nlines; i++)
    fprintf(linestat_file, "%g ", CLIGHT/linelist[i].nu);
  fprintf(linestat_file,"\n");

  for (int i = 0; i < nlines; i++)
    fprintf(linestat_file, "%d ", get_element(linelist[i].elementindex));
  fprintf(linestat_file, "\n");

  for (int i = 0; i < nlines; i++)
    fprintf(linestat_file, "%d ", get_ionstage(linelist[i].elementindex, linelist[i].ionindex));
  fprintf(linestat_file, "\n");

  for (int i = 0; i < nlines; i++)
    fprintf(linestat_file, "%d ", linelist[i].upperlevelindex + 1);
  fprintf(linestat_file, "\n");

  for (int i = 0; i < nlines; i++)
    fprintf(linestat_file, "%d ", linelist[i].lowerlevelindex + 1);
  fprintf(linestat_file, "\n");

  fflush(linestat_file);
  //setvbuf(linestat_file, NULL, _IOLBF, 1); // flush after every line makes it slow!
}


static void pkt_action_counters_reset(void)
{
  ma_stat_activation_collexc = 0;
  ma_stat_activation_collion = 0;
  ma_stat_activation_ntcollexc = 0;
  ma_stat_activation_ntcollion = 0;
  ma_stat_activation_bb = 0;
  ma_stat_activation_bf = 0;
  ma_stat_activation_fb = 0;
  ma_stat_deactivation_colldeexc = 0;
  ma_stat_deactivation_collrecomb = 0;
  ma_stat_deactivation_bb = 0;
  ma_stat_deactivation_fb = 0;
  ma_stat_internaluphigher = 0;
  ma_stat_internaluphighernt = 0;
  ma_stat_internaldownlower = 0;
  k_stat_to_ma_collexc = 0;
  k_stat_to_ma_collion = 0;
  k_stat_to_r_ff = 0;
  k_stat_to_r_fb = 0;
  k_stat_to_r_bb = 0;
  k_stat_from_ff = 0;
  k_stat_from_bf = 0;
  k_stat_from_earlierdecay = 0;

  nt_reset_stats();

  escounter = 0;
  cellcrossings = 0;
  updatecellcounter = 0;
  coolingratecalccounter = 0;
  resonancescatterings = 0;
  upscatter = 0;
  downscatter = 0;
  nesc = 0;
}


static void pkt_action_counters_printout(const PKT *const pkt, const int nts)
{
  int allpktinteractions = 0;
  for (int i = 0; i < npkts; i++)
  {
    allpktinteractions += pkt[i].interactions;
  }
  const double meaninteractions = (double) allpktinteractions / npkts;
  printout("mean number of interactions per packet = %g\n", meaninteractions);

  const double deltat = time_step[nts].width;
  double modelvolume = 0.;
  for (int mgi = 0; mgi < npts_model; mgi++)
  {
    modelvolume += vol_init_modelcell(mgi) * pow(time_step[nts].mid / tmin, 3);
  }

  /// Printout packet statistics
  printout("ma_stat_activation_collexc = %d\n", ma_stat_activation_collexc);
  printout("ma_stat_activation_collion = %d\n", ma_stat_activation_collion);
  printout("ma_stat_activation_ntcollexc = %d\n", ma_stat_activation_ntcollexc);
  printout("ma_stat_activation_ntcollion = %d\n", ma_stat_activation_ntcollion);
  printout("ma_stat_activation_bb = %d\n", ma_stat_activation_bb);
  printout("ma_stat_activation_bf = %d\n", ma_stat_activation_bf);
  printout("ma_stat_activation_fb = %d\n", ma_stat_activation_fb);
  printout("ma_stat_deactivation_colldeexc = %d\n", ma_stat_deactivation_colldeexc);
  printout("ma_stat_deactivation_collrecomb = %d\n", ma_stat_deactivation_collrecomb);
  printout("ma_stat_deactivation_bb = %d\n", ma_stat_deactivation_bb);
  printout("ma_stat_deactivation_fb = %d\n", ma_stat_deactivation_fb);
  printout("ma_stat_internaluphigher = %d\n", ma_stat_internaluphigher);
  printout("ma_stat_internaluphighernt = %d\n", ma_stat_internaluphighernt);
  printout("ma_stat_internaldownlower = %d\n", ma_stat_internaldownlower);

  printout("k_stat_to_ma_collexc = %d\n", k_stat_to_ma_collexc);
  printout("k_stat_to_ma_collion = %d\n", k_stat_to_ma_collion);
  printout("k_stat_to_r_ff = %d\n", k_stat_to_r_ff);
  printout("k_stat_to_r_fb = %d\n", k_stat_to_r_fb);
  printout("k_stat_to_r_bb = %d\n", k_stat_to_r_bb);
  printout("k_stat_from_ff = %d\n", k_stat_from_ff);
  printout("k_stat_from_bf = %d\n", k_stat_from_bf);
  printout("k_stat_from_earlierdecay = %d\n", k_stat_from_earlierdecay);

  nt_print_stats(nts, modelvolume, deltat);

  printout("escounter = %d\n", escounter);
  printout("cellcrossing  = %d\n", cellcrossings);
  printout("updatecellcounter  = %d\n", updatecellcounter);
  printout("coolingratecalccounter = %d\n", coolingratecalccounter);
  printout("resonancescatterings  = %d\n", resonancescatterings);

  printout("upscatterings  = %d\n", upscatter);
  printout("downscatterings  = %d\n", downscatter);
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
      radfield_MPI_Bcast(modelgridindex, root);

      if (get_numassociatedcells(modelgridindex) > 0)
      {
        nt_MPI_Bcast(modelgridindex, root);
        if (NLTE_POPS_ON)
        {
          MPI_Bcast(modelgrid[modelgridindex].nlte_pops, total_nlte_levels, MPI_DOUBLE, root, MPI_COMM_WORLD);
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
        //if (cell[nn].rho > MINDENSITY)
        if (get_numassociatedcells(mgi) > 0)
        {
          MPI_Pack(&modelgrid[mgi].Te, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&modelgrid[mgi].TR, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&modelgrid[mgi].TJ, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&modelgrid[mgi].W, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&modelgrid[mgi].rho, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&modelgrid[mgi].nne, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&modelgrid[mgi].nnetot, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&modelgrid[mgi].kappagrey, 1, MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&modelgrid[mgi].totalcooling, 1, MPI_DOUBLE, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          MPI_Pack(&modelgrid[mgi].thick, 1, MPI_SHORT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);

          for (int element = 0; element < nelements; element++)
          {
            MPI_Pack(modelgrid[mgi].composition[element].groundlevelpop, get_nions(element), MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
            MPI_Pack(modelgrid[mgi].composition[element].partfunct, get_nions(element), MPI_FLOAT, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
            MPI_Pack(modelgrid[mgi].cooling[element].contrib, get_nions(element), MPI_DOUBLE, mpi_grid_buffer, mpi_grid_buffer_size, &position, MPI_COMM_WORLD);
          }
        }
      }
      printout("mem_usage: MPI_BUFFER: used %d of %d bytes allocated to mpi_grid_buffer\n", position, mpi_grid_buffer_size);
      assert(position <= mpi_grid_buffer_size);
    }
    MPI_Bcast(mpi_grid_buffer, mpi_grid_buffer_size, MPI_PACKED, root, MPI_COMM_WORLD);

    position = 0;
    int nlp;
    MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &nlp, 1, MPI_INT, MPI_COMM_WORLD);
    for (int nn = 0; nn < nlp; nn++)
    {
      int mgi;
      MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &mgi, 1, MPI_INT, MPI_COMM_WORLD);
      //if (cell[ncl].rho > MINDENSITY)
      if (get_numassociatedcells(mgi) > 0)
      {
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &modelgrid[mgi].Te, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &modelgrid[mgi].TR, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &modelgrid[mgi].TJ, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &modelgrid[mgi].W, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &modelgrid[mgi].rho, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &modelgrid[mgi].nne, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &modelgrid[mgi].nnetot, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &modelgrid[mgi].kappagrey, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &modelgrid[mgi].totalcooling, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, &modelgrid[mgi].thick, 1, MPI_SHORT, MPI_COMM_WORLD);

        for (int element = 0; element < nelements; element++)
        {
          MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, modelgrid[mgi].composition[element].groundlevelpop, get_nions(element), MPI_FLOAT, MPI_COMM_WORLD);
          MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, modelgrid[mgi].composition[element].partfunct, get_nions(element), MPI_FLOAT, MPI_COMM_WORLD);
          MPI_Unpack(mpi_grid_buffer, mpi_grid_buffer_size, &position, modelgrid[mgi].cooling[element].contrib, get_nions(element), MPI_DOUBLE, MPI_COMM_WORLD);
        }
      }
    }
  }

  #ifndef FORCE_LTE
    #if (!NO_LUT_PHOTOION)
      if ((!simulation_continued_from_saved) || (nts - itstep != 0) || (titer != 0))
      {
        /// Reduce the corrphotoionrenorm array.
        printout("nts %d, titer %d: bcast corr photoionrenorm\n", nts, titer);
        MPI_Allreduce(MPI_IN_PLACE, &corrphotoionrenorm, MMODELGRID * nelements * maxion, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /// Reduce the gammaestimator array. Only needed to write restart data.
        printout("nts %d, titer %d: bcast gammaestimator\n", nts, titer);
        MPI_Allreduce(MPI_IN_PLACE, &gammaestimator, MMODELGRID * nelements * maxion, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
    #endif
  #endif
  MPI_Barrier(MPI_COMM_WORLD);
}


static void mpi_reduce_estimators(int my_rank, int nts)
{
  radfield_reduce_estimators();
  #ifndef FORCE_LTE
    MPI_Allreduce(MPI_IN_PLACE, &ffheatingestimator, MMODELGRID, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &colheatingestimator, MMODELGRID, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #if (!NO_LUT_PHOTOION)
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &gammaestimator, MMODELGRID * nelements * maxion, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif
    #if (!NO_LUT_BFHEATING)
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &bfheatingestimator, MMODELGRID * nelements * maxion, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif

    // MPI_Reduce(MPI_IN_PLACE, &ionfluxestimator, MMODELGRID*nelements*maxion, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // MPI_Reduce(MPI_IN_PLACE, &twiddle, MMODELGRID*nelements*maxion, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // MPI_Reduce(MPI_IN_PLACE, &stimrecombestimator, MMODELGRID*nelements*maxion, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // MPI_Reduce(&mabfcount, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // if (my_rank == 0)
    // {
    //   for (int i = 0; i < MMODELGRID; i++)
    //   {
    //     mabfcount[i] = redhelper[i]/p;
    //   }
    // }
    // MPI_Reduce(&mabfcount_thermal, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // if (my_rank == 0)
    // {
    //   for (int i = 0; i < MMODELGRID; i++)
    //   {
    //     mabfcount_thermal[i] = redhelper[i]/p;
    //   }
    // }
    // MPI_Reduce(&kbfcount, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // if (my_rank == 0)
    // {
    //   for (int i = 0; i < MMODELGRID; i++)
    //   {
    //     kbfcount[i] = redhelper[i]/p;
    //   }
    // }
    // MPI_Reduce(&kbfcount_ion, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // if (my_rank == 0)
    // {
    //   for (int i = 0; i < MMODELGRID; i++)
    //   {
    //     kbfcount_ion[i] = redhelper[i]/p;
    //   }
    // }
    // MPI_Reduce(&kffcount, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // if (my_rank == 0)
    // {
    //   for (int i = 0; i < MMODELGRID; i++)
    //   {
    //     kffcount[i] = redhelper[i]/p;
    //   }
    // }
    // MPI_Reduce(&kffabs, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // if (my_rank == 0)
    // {
    //   for (int i = 0; i < MMODELGRID; i++)
    //   {
    //     kffabs[i] = redhelper[i]/p;
    //   }
    // }
    // MPI_Reduce(&kbfabs, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // if (my_rank == 0)
    // {
    //   for (int i = 0; i < MMODELGRID; i++)
    //   {
    //     kbfabs[i] = redhelper[i]/p;
    //   }
    // }
    // MPI_Reduce(&kgammadep, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // if (my_rank == 0)
    // {
    //   for (int i = 0; i < MMODELGRID; i++)
    //   {
    //     kgammadep[i] = redhelper[i]/p;
    //   }
    // }
  #endif

  #ifdef RECORD_LINESTAT
    MPI_Allreduce(MPI_IN_PLACE, ecounter, nlines, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, acounter, nlines, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #endif

  //double deltaV = pow(wid_init * time_step[nts].mid/tmin, 3.0);
  //double deltat = time_step[nts].width;
  if (do_rlc_est != 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, &rpkt_emiss, MMODELGRID, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (do_comp_est)
  {
    MPI_Allreduce(MPI_IN_PLACE, &compton_emiss, MMODELGRID * EMISS_MAX, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  }

  /// Communicate gamma and positron deposition and write to file
  MPI_Allreduce(MPI_IN_PLACE, &time_step[nts].gamma_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &time_step[nts].positron_dep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  time_step[nts].gamma_dep /= nprocs;
  time_step[nts].positron_dep /= nprocs;

  #if TRACK_ION_STATS
  MPI_Allreduce(MPI_IN_PLACE, ionstats, npts_model * includedions * ION_COUNTER_COUNT, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif

  MPI_Barrier(MPI_COMM_WORLD);
}
#endif


static void read_temp_packetsfile(const int timestep, const int my_rank, PKT *const pkt)
{
  char filename[100];
  sprintf(filename, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  printout("Reading %s...", filename);
  FILE *packets_file = fopen_required(filename, "rb");
  fread(pkt, sizeof(PKT), npkts, packets_file);
  //read_packets(packets_file);
  fclose(packets_file);
  printout("done\n");
}


static void write_temp_packetsfile(const int timestep, const int my_rank, const PKT *const pkt)
{
  char filename[100];
  sprintf(filename, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  printout("Writing %s...", filename);
  FILE *packets_file = fopen_required(filename, "wb");

  fwrite(pkt, sizeof(PKT), npkts, packets_file);
  fclose(packets_file);
  printout("done\n");
}


static void remove_temp_packetsfile(const int timestep, const int my_rank)
{
  char filename[100];
  sprintf(filename, "packets_%.4d_ts%d.tmp", my_rank, timestep);

  if (!access(filename, F_OK))
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
  {
    remove(prevfilename);
    printout("Deleted %s\n", prevfilename);
  }
}


#if (TRACK_ION_STATS)
void increment_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstatscounters ion_counter_type, const double increment)
{
  if (!TRACK_ION_MASTATS && (ion_counter_type >= 18))
  {
    return;
  }

  assert(ion < get_nions(element));
  assert(ion_counter_type < ION_COUNTER_COUNT);

  const int uniqueionindex = get_uniqueionindex(element, ion);
  #ifdef _OPENMP
    #pragma omp atomic
  #endif
  ionstats[modelgridindex * includedions * ION_COUNTER_COUNT + uniqueionindex * ION_COUNTER_COUNT + ion_counter_type] += increment;
}


double get_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstatscounters ion_counter_type)
{
  assert(ion < get_nions(element));
  assert(ion_counter_type < ION_COUNTER_COUNT);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  return ionstats[modelgridindex * includedions * ION_COUNTER_COUNT + uniqueionindex * ION_COUNTER_COUNT + ion_counter_type];
}


void set_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstatscounters ion_counter_type, const double newvalue)
{
  assert(ion < get_nions(element));
  assert(ion_counter_type < ION_COUNTER_COUNT);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  #ifdef _OPENMP
    #pragma omp atomic
  #endif
  ionstats[modelgridindex * includedions * ION_COUNTER_COUNT + uniqueionindex * ION_COUNTER_COUNT + ion_counter_type] = newvalue;
}
#endif


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


void* reallocmanaged(void* ptr, size_t newSize, size_t curSize)
{
  void *newptr;
  cudaMallocManaged(&newptr, newSize);
  if (ptr == 0)
  {
    return newptr;
  }
  else
  {
    size_t transferSize = newSize > curSize ? curSize : newSize;
    cudaMemcpy(newptr, ptr, transferSize, cudaMemcpyDefault);
    cudaDeviceSynchronize();
    cudaFree(ptr);
    return newptr;
  }
}


void* makemanaged(void* ptr, size_t curSize)
{
  void *newptr;
  cudaMallocManaged(&newptr, curSize);
  cudaMemcpy(newptr, ptr, curSize, cudaMemcpyDefault);
  cudaDeviceSynchronize();
  free(ptr);
  return newptr;
}

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

  packets_file = fopen_required(filename, "wb");

  write_vspecpol(packets_file);
  fclose(packets_file);

  // Write temporary files for vpkt_grid
  if (vgrid_flag == 1)
  {
    if (nts % 2 == 0)
      sprintf(filename,"vpkt_grid_%d_%d_even.tmp", 0, my_rank);
    else
      sprintf(filename,"vpkt_grid_%d_%d_odd.tmp", 0, my_rank);

    packets_file = fopen_required(filename, "wb");

    write_vpkt_grid(packets_file);
    fclose(packets_file);
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

static bool do_timestep(const int outer_iteration, const int nts, const int titer, const int my_rank, PKT* packets, const int walltimelimitseconds)
{
  bool do_this_full_loop = true;

  const int nts_prev = (titer != 0 || nts == 0) ? nts : nts - 1;
  if ((titer > 0) || (simulation_continued_from_saved && (nts == itstep)))
  {
    /// Read the packets file to reset before each additional iteration on the timestep
    read_temp_packetsfile(nts, my_rank, packets);
  }

  /// Some counters on pkt-actions need to be reset to do statistics
  pkt_action_counters_reset();

  if (nts == 0)
  {
    initialise_prev_titer_photoionestimators();
  }

  #ifdef RECORD_LINESTAT
    // The same for absorption/emission of r-pkts in lines
    for (int i = 0; i < nlines; i++)
    {
      acounter[i] = 0;
      ecounter[i] = 0;
    }
  #endif

  do_comp_est = do_r_lc ? false : estim_switch(nts);

  // Update the matter quantities in the grid for the new timestep.

  const time_t sys_time_start_update_grid = time(NULL);
  printout("\ntimestep %d: time before update grid %ld (tstart + %ld)\n",
           nts, sys_time_start_update_grid, sys_time_start_update_grid - real_time_start);

  #ifndef FORCE_LTE
    #if (!NO_LUT_PHOTOION)
      /// Initialise corrphotoionrenorm[i] to zero before update_grid is called
      /// This allows reduction after update_grid has finished
      /// unless they have been read from file and must neither be touched
      /// nor broadcasted after update_grid
      if ((!simulation_continued_from_saved) || (nts - itstep != 0) || (titer != 0))
      {
        printout("nts %d, titer %d: reset corr photoionrenorm\n",nts,titer);
        for (int i = 0; i < MMODELGRID * nelements * maxion; i++)
        {
          corrphotoionrenorm[i] = 0.;
        }
        printout("after nts %d, titer %d: reset corr photoionrenorm\n",nts,titer);
      }
    #endif
  #endif

  update_grid(estimators_file, nts, nts_prev, my_rank, nstart, ndo, titer);

  const time_t sys_time_finish_update_grid = time(NULL);
  printout("timestep %d: update_grid: process %d finished update_grid at %ld (took %ld seconds)\n",
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
    mpi_communicate_grid_properties(my_rank, nprocs, nstart, ndo, nts, titer, mpi_grid_buffer, mpi_grid_buffer_size);
  #endif

  printout("timestep %d: time after grid properties have been communicated %ld (took %ld seconds)\n",
           nts, time(NULL), time(NULL) - sys_time_start_communicate_grid);

  /// If this is not the 0th time step of the current job step,
  /// write out a snapshot of the grid properties for further restarts
  /// and update input.txt accordingly
  if (((nts - itstep) != 0))
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
  if ((nts < ftstep) && do_this_full_loop)
  {
    /// Now process the packets.
    const time_t time_update_packets_start = time(NULL);
    printout("timestep %d: time before update packets %ld\n", nts, time_update_packets_start);

    update_packets(nts, packets);

    pkt_action_counters_printout(packets, nts);

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
    if (do_comp_est)
    {
      normalise_compton_estimators(nts);
      if (my_rank == 0)
      {
        write_compton_estimators(nts);
      }
    }

    if (do_rlc_est != 0)
    {
      normalise_grey(nts);
      if ((do_rlc_est != 3) && (my_rank == 0))
      {
        write_grey(nts);
      }
    }

    if (my_rank == 0)
    {
      FILE *dep_file = fopen_required("deposition.out", "w");
      for (int i = 0; i <= nts; i++)
      {
        fprintf(dep_file, "%g %g %g %g\n", time_step[i].mid / DAY,
                time_step[i].gamma_dep / time_step[i].width / LSUN,
                time_step[i].positron_dep / time_step[i].width / LSUN,
                (time_step[i].gamma_dep + time_step[i].positron_dep) / time_step[i].width / LSUN);
      }
      fclose(dep_file);
    }

    #ifdef MPI_ON
      printout("timestep %d: time after estimators have been communicated %d (took %ld seconds)\n", nts, time(NULL), time(NULL) - time_communicate_estimators_start);
    #endif

    printout("%d: During timestep %d on MPI process %d, %d pellets decayed and %d packets escaped. (t=%gd)\n",
             outer_iteration, nts, my_rank, time_step[nts].pellet_decays, nesc, time_step[nts].mid / DAY);

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
        for (int i = 0; i < nlines; i++)
          fprintf(linestat_file, "%d ", ecounter[i]);
        fprintf(linestat_file, "\n");
        for (int i = 0; i < nlines; i++)
          fprintf(linestat_file, "%d ", acounter[i]);
        fprintf(linestat_file, "\n");
        fflush(linestat_file);
      }
    #endif

    if (nts == ftstep - 1)
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
    *ndo = npts_model;
    return;
  #endif

  int n_leftover = 0;

  int nblock = npts_model / nprocesses; // integer division, minimum cells for any process
  const int numtot = nblock * nprocesses; // cells counted if all processes do the minimum number of cells
  if (numtot > npts_model) // LJS: should never be the case?
  {
    nblock = nblock - 1;
    *maxndo = nblock + 1;
    n_leftover = npts_model - (nblock * nprocesses);
  }
  else if (numtot < npts_model)
  {
    *maxndo = nblock + 1;
    n_leftover = npts_model - (nblock * nprocesses);
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

  nprocs = p;              /// Global variable which holds the number of MPI processes
  rank_global = my_rank;   /// Global variable which holds the rank of the active MPI process

#ifdef _OPENMP
  /// Explicitly turn off dynamic threads because we use the threadprivate directive!!!
  omp_set_dynamic(0);
  #pragma omp parallel private(filename)
#endif
  {
    /// Get the current threads ID, copy it to a threadprivate variable
    tid = omp_get_thread_num();
    /// and initialise the threads outputfile
    sprintf(filename,"output_%d-%d.txt", my_rank, tid);
    output_file = fopen_required(filename, "w");
    /// Makes sure that the output_file is written line-by-line
    setvbuf(output_file, NULL, _IOLBF, 1);

    /// Get the total number of active threads
    nthreads = omp_get_num_threads();
    if (nthreads > MTHREADS)
    {
      printout("[Fatal] too many threads. Set MTHREADS (%d) > nthreads (%d). Abort.\n", MTHREADS, nthreads);
      abort();
    }
#   ifdef _OPENMP
    printout("OpenMP parallelisation active with %d threads\n", nthreads);
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

  PKT *packets;
  #if CUDA_ENABLED
  cudaMallocManaged(&packets, MPKTS * sizeof(PKT));
  #else
  packets = (PKT *) calloc(MPKTS, sizeof(PKT));
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
    printout("MPI enabled\n");
  #else
    printout("MPI disabled\n");
  #endif

  #if CUDA_ENABLED
  printout("NVIDIA CUDA is enabled\n");
  int deviceCount = 0;
  assert(cudaGetDeviceCount(&deviceCount) == cudaSuccess);
  printout("Detected %d CUDA capable device(s)\n", deviceCount);
  assert(deviceCount > 0);
  #endif

  if ((mastate = (mastate_t *) calloc(nthreads, sizeof(mastate_t))) == NULL)
  {
    printout("[fatal] input: error initializing macro atom state variables ... abort\n");
    abort();
  }
  if ((kappa_rpkt_cont = (rpkt_cont_opacity_struct *) calloc(nthreads, sizeof(rpkt_cont_opacity_struct))) == NULL)
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

  /// Initialise linestat file
  #ifdef RECORD_LINESTAT
  if (my_rank == 0)
  {
    initialise_linestat_file();
  }
  #endif

  printout("time after input %ld\n", time(NULL));
  printout("simulation propagates %d packets through a %d x %d x %d grid\n",
           npkts, ncoordgrid[0], ncoordgrid[1], ncoordgrid[2]);
  printout("timesteps %d\n", ntstep);

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

  #if TRACK_ION_STATS
  ionstats = calloc(npts_model * includedions * ION_COUNTER_COUNT, sizeof(double));
  #endif

  /// As a precaution, explicitly zero all the estimators here
  zero_estimators();

  printout("time after zero estimators %ld\n", time(NULL));

  /// Record the chosen syn_dir
  FILE *syn_file = fopen_required("syn_dir.txt", "w");
  fprintf(syn_file, "%g %g %g", syn_dir[0], syn_dir[1], syn_dir[2]);
  fclose(syn_file);
  printout("time read syn file %ld\n", time(NULL));

  bool terminate_early = false;
  file_set = false;
  debuglevel = 4;  /// Selects detail level of debug output, needs still some work.
  for (int outer_iteration = 0; outer_iteration < n_out_it; outer_iteration++)
  {
    /// Initialise the grid. Call routine that sets up the initial positions
    /// and sizes of the grid cells.
    time_init();
    printout("time time init %ld\n", time(NULL));
    grid_init(my_rank);

    printout("mem_usage: packets occupy %.1f MB\n", MPKTS * sizeof(PKT) / 1024. / 1024.);

    if (!simulation_continued_from_saved)
    {
      /// Next we want to initialise the packets.
      /// To overcome memory limitations for large numbers of packets, which need to be
      /// propagated on the _same_ grid, this middle_iteration loop was introduced.
      for (int middle_iteration = 0; middle_iteration < n_middle_it; middle_iteration++)
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
    printout("process rank %d (of %d) doing %d cells", my_rank, nprocs, ndo);
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
      mpi_grid_buffer_size = 4 * ((12 + 4 * includedions) * (maxndo) + 1);
      printout("reserve mpi_grid_buffer_size %d space for MPI communication buffer\n", mpi_grid_buffer_size);
      //char buffer[mpi_grid_buffer_size];
      mpi_grid_buffer  = malloc(mpi_grid_buffer_size * sizeof(char));
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
    //for (nts = itstep; nts < ftstep; nts++)
    //{

    /// Now use while loop to allow for timed restarts
    const int last_loop = ftstep;
    int nts = itstep;

    // Initialise virtual packets file and vspecpol
    #ifdef VPKT_ON
      if (!simulation_continued_from_saved)
      {
        // New simulation

        init_vspecpol();

        if (vgrid_flag == 1)
          init_vpkt_grid();
      }
      else
      {
        // Continue simulation: read into temporary files

        if (nts % 2 == 0)
          sprintf(filename,"vspecpol_%d_%d_odd.tmp", 0, my_rank);
        else
          sprintf(filename,"vspecpol_%d_%d_even.tmp", 0 ,my_rank);

        FILE *packets_file = fopen_required(filename, "rb");

        read_vspecpol(packets_file);

        if (vgrid_flag == 1)
        {
          if (nts % 2 == 0)
            sprintf(filename,"vpkt_grid_%d_%d_odd.tmp", 0, my_rank);
          else
            sprintf(filename,"vpkt_grid_%d_%d_even.tmp", 0, my_rank);

          packets_file = fopen_required(filename, "rb");

          read_vpkt_grid(packets_file);
        }
      }
    #endif

    while (nts < last_loop && !terminate_early)
    {
      nts_global = nts;
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
        n_titer = 1;
        initial_iteration = (nts < n_lte_timesteps);
      #endif

      for (int titer = 0; titer < n_titer; titer++)
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

  if ((ntstep != ftstep) || (terminate_early))
  {
    printout("RESTART_NEEDED to continue model\n");
  }
  else
  {
    printout("No need for restart\n");
  }


  printout("simulation finished at %ld\n", time(NULL));
  //fclose(tb_file);
  fclose(estimators_file);
  macroatom_close_file();
  if (NLTE_POPS_ON)
    nltepop_close_file();
  radfield_close_file();
  nt_close_file();

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
    gsl_integration_workspace_free(gslworkspace);
  }

  #ifdef MPI_ON
    MPI_Finalize();
  #endif

  free(packets);
  #if TRACK_ION_STATS
  free(ionstats);
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

  double t_current = time_step[timestep].mid;
  float T_e = cell[cellnumber].T_e;
  double T_R = cell[cellnumber].T_R;
  double W = cell[cellnumber].W;



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
        lineindex = elements[element].ions[ion].levels[lower].uptrans_lineindicies[i];
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
