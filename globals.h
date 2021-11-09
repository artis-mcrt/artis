#ifndef GLOBALS_H
#define GLOBALS_H

#if CUDA_ENABLED
#include <curand_kernel.h>
#endif

#include "grid.h"
#include "types.h"
#include "cuda.h"
#ifdef MPI_ON
  #include "mpi.h"
#endif

namespace globals
{
  #if CUDA_ENABLED
  extern __managed__ curandState curandstates[MCUDATHREADS];
  #endif

  extern __managed__ double syn_dir[3]; // vector pointing from origin to observer

  extern __managed__ struct time *time_step;

  #define MSYN_TIME 100
  extern __managed__ int nsyn_time;
  extern __managed__ double time_syn[MSYN_TIME];

  #define EMISS_MAX 2
  extern __managed__ int emiss_offset;
  extern __managed__ int emiss_max;

  extern __managed__ float *compton_emiss;
  extern __managed__ double *rpkt_emiss;

  #if (!NO_LUT_PHOTOION)
    extern __managed__ double *corrphotoionrenorm;
    extern __managed__ double *gammaestimator;
  #endif
  #if (!NO_LUT_BFHEATING)
    extern __managed__ double *bfheatingestimator;
  #endif
  extern __managed__ double *ffheatingestimator;
  extern __managed__ double *colheatingestimator;
  #ifndef FORCE_LTE
    #ifdef DO_TITER
      extern __managed__ double *gammaestimator_save;
      extern __managed__ double *bfheatingestimator_save;
      extern __managed__ double *ffheatingestimator_save;
      extern __managed__ double *colheatingestimator_save;
    #endif
  #endif

  #ifdef RECORD_LINESTAT
    extern __managed__ int *ecounter;
    extern __managed__ int *acounter;
    extern __managed__ int *linestat_reduced;
  #endif

  extern __managed__ int nprocs_exspec;
  extern __managed__ bool do_emission_res;

  extern __managed__ bool file_set;

  extern __managed__ bool do_comp_est;
  extern __managed__ bool do_r_lc;
  extern __managed__ int do_rlc_est;

  extern __managed__ double CLIGHT_PROP;

  extern __managed__ double gamma_grey;

  #define GREY_OP 0.1

  extern __managed__ double max_path_step;

  extern __managed__ int opacity_case;

  extern __managed__ int nlines;
  extern __managed__ int includedions;
  extern __managed__ elementlist_entry *elements;
  extern __managed__ linelist_entry *linelist;
  extern __managed__ bflist_t *bflist;

  extern __managed__ rpkt_cont_opacity_struct *kappa_rpkt_cont;

  extern __managed__ int ncoolingterms;

  extern __managed__ double *allcont_nu_edge;
  extern __managed__ fullphixslist_t *allcont;
  extern __managed__ phixslist_t *phixslist;
  extern __managed__ int nbfcontinua;
  extern __managed__ int nbfcontinua_ground;
  extern __managed__ int NPHIXSPOINTS;
  extern __managed__ double NPHIXSNUINCREMENT;

  extern __managed__ cellhistory_struct *cellhistory;

  extern __managed__ int debuglevel;

  #ifdef MPI_ON
  extern MPI_Comm mpi_comm_node;
  extern MPI_Comm mpi_comm_internode;
  #endif

  extern __managed__ int nprocs;
  extern __managed__ int rank_global;

  extern __managed__ int node_nprocs;
  extern __managed__ int rank_in_node;

  extern __managed__ int node_count;
  extern __managed__ int node_id;

  extern __managed__ int npkts;
  extern __managed__ int nesc;

  extern __managed__ double coordmax[3];
  extern __managed__ double vmax;
  extern __managed__ double rmax;
  extern __managed__ double tmax;
  extern __managed__ double tmin;

  extern __managed__ int ntstep;
  extern __managed__ int itstep;
  extern __managed__ int ftstep;
  extern __managed__ int nts_global;

  extern __managed__ int nnubins;
  extern __managed__ double nu_min_r;
  extern __managed__ double nu_max_r;

  extern __managed__ double nusyn_min, nusyn_max;
  extern __managed__ int nfake_gam;

  extern __managed__ double opcase3_normal;
  extern __managed__ double rho_crit_para;
  extern __managed__ double rho_crit;

  extern __managed__ int debug_packet;

  extern __managed__ int total_nlte_levels;

  extern __managed__ bool homogeneous_abundances;

  extern __managed__ bool simulation_continued_from_saved;
  extern __managed__ double nu_rfcut;
  extern __managed__ int n_lte_timesteps;
  extern __managed__ double cell_is_optically_thick;
  extern __managed__ int n_grey_timesteps;
  extern __managed__ int n_titer;
  extern __managed__ bool initial_iteration;
  extern __managed__ int max_bf_continua;
  extern __managed__ int n_kpktdiffusion_timesteps;
  extern __managed__ float kpktdiffusion_timescale;

} // namespace globals

#endif // GLOBALS_H