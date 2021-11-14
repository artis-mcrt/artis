#include "types.h"
#include "sn3d.h"
#include "globals.h"
#ifdef MPI_ON
  #include "mpi.h"
#endif

namespace globals
{

  #if CUDA_ENABLED
  __managed__ curandState curandstates[MCUDATHREADS];
  #endif

  __managed__ double syn_dir[3]; // vector pointing from origin to observer

  //#define NRAYS_SYN 1 // number of rays traced in a syn calculation

  //RAY rays[NRAYS_SYN];
  __managed__ struct time *time_step = NULL;

  __managed__ int nsyn_time;
  __managed__ double time_syn[MSYN_TIME];

  __managed__ int emiss_offset;   // the index in the line list of the 1st line for which
                                  // an emissivity estimator is recorded
  __managed__ int emiss_max;      // actual number of frequency points in emissivity grid


  /// THESE ARE THE GRID BASED ESTIMATORS
  __managed__ float *compton_emiss = NULL;  /// Volume estimator for the compton emissivity
  __managed__ double *rpkt_emiss = NULL;                /// Volume estimator for the rpkt emissivity


  #if (!NO_LUT_PHOTOION)
    __managed__ double *corrphotoionrenorm = NULL;
    __managed__ double *gammaestimator = NULL;
  #endif
  #if (!NO_LUT_BFHEATING)
    __managed__ double *bfheatingestimator = NULL;
  #endif
  __managed__ double *ffheatingestimator = NULL;
  __managed__ double *colheatingestimator = NULL;
  #ifndef FORCE_LTE
    #ifdef DO_TITER
      __managed__ double *gammaestimator_save = NULL;
      __managed__ double *bfheatingestimator_save = NULL;
      __managed__ double *ffheatingestimator_save = NULL;
      __managed__ double *colheatingestimator_save = NULL;
    #endif
  #endif

  #ifdef RECORD_LINESTAT
    __managed__ int *ecounter = NULL;
    __managed__ int *acounter = NULL;
    __managed__ int *linestat_reduced = NULL;
  #endif


  __managed__ int nprocs_exspec = 1;
  __managed__ bool do_emission_res = 1;

  __managed__ bool file_set; // 1 if the output files already exist. 0 otherwise.

  __managed__ bool do_comp_est; // 1 = compute compton emissivity estimators. 0 = don't
  __managed__ bool do_r_lc;     // If not set to 1 then the opacity for r-packets is 0.
  __managed__ int do_rlc_est;  // 1 = compute estimators for the r-pkt light curve.
                               // 2 = compute estimators with opacity weights
                               // 3 = compute estimators, but use only for gamma-heating rate


  __managed__ double CLIGHT_PROP; // Speed of light for ray travel. Physically = CLIGHT but
                                  // can be changed for testing.

  __managed__ double gamma_grey; // set to -ve for proper treatment. If possitive, then
                                 // gamma_rays are treated as grey with this opacity.

  __managed__ double max_path_step;

  __managed__ int opacity_case; // 0 normally, 1 for Fe-grp dependence.
                    ///MK: 2 for Fe-grp dependence and proportional to 1/rho
                    ///MK: 3 combination of 1 & 2 depending on a rho_crit
                    ///MK: 4 non-grey treatment



  /// ATOMIC DATA

  __managed__ int nlines;
  __managed__ int includedions;
  __managed__ elementlist_entry *elements = NULL;
  __managed__ linelist_entry *linelist = NULL;
  __managed__ bflist_t *bflist = NULL;

  __managed__ rpkt_cont_opacity_struct *kappa_rpkt_cont = NULL;

  /// Coolinglist
  __managed__ int ncoolingterms;

  /// PHIXSLIST

  __managed__ double *allcont_nu_edge = NULL;
  __managed__ fullphixslist_t *allcont = NULL;
  __managed__ phixslist_t *phixslist = NULL;
  __managed__ int nbfcontinua;
  __managed__ int nbfcontinua_ground; ///number of bf-continua
  __managed__ int NPHIXSPOINTS;
  __managed__ double NPHIXSNUINCREMENT;

  __managed__ cellhistory_struct *cellhistory = NULL;

  __managed__ int debuglevel;

  #ifdef MPI_ON
  MPI_Comm mpi_comm_node = NULL;
  MPI_Comm mpi_comm_internode = NULL;
  #endif

  __managed__ int nprocs = -1;         // number of MPI processes
  __managed__ int rank_global = -1;    // rank of the active MPI process

  __managed__ int node_nprocs = -1;    // number of MPI processes on this node
  __managed__ int rank_in_node = -1;   // local rank within this node

  __managed__ int node_count = -1;      // number of MPI nodes
  __managed__ int node_id = -1;         // unique number for each node

  __managed__ int npkts = -1;
  __managed__ int nesc = 0; //number of packets that escape during current timestep

  __managed__ double coordmax[3];
  __managed__ double vmax;
  __managed__ double rmax;  /// Total mass and outer velocity/radius
  __managed__ double tmax = -1.;              /// End time of current simulation
  __managed__ double tmin = -1.;              /// Start time of current simulation

  __managed__ int ntstep = -1;       /// Number of timesteps
  __managed__ int itstep = -1;       /// Initial timestep's number
  __managed__ int ftstep = -1;       /// Final timestep's number
  __managed__ int nts_global = -1;   /// Current time step

  __managed__ int nnubins = -1; //number of bins for spectrum
  __managed__ double nu_min_r = -1.;
  __managed__ double nu_max_r = -1.; //limits on frequency range for r-pkt spectrum

  __managed__ double nusyn_min;
  __managed__ double nusyn_max; //limits on range for syn
  __managed__ int nfake_gam; //# of fake gamma ray lines for syn

  /// New variables for other opacity cases, still grey.
  __managed__ double opcase3_normal;           ///MK: normalisation factor for opacity_case 3
  __managed__ double rho_crit_para;            ///MK: free parameter for the selection of the critical opacity in opacity_case 3
  __managed__ double rho_crit;                 ///MK: critical opacity in opacity_case 3 (could now be declared locally)

  /// New variables for the non-grey case
  __managed__ int debug_packet;                /// activate debug output for this packet if non negative

  __managed__ int total_nlte_levels;            ///total number of nlte levels

  __managed__ bool homogeneous_abundances;

  __managed__ bool simulation_continued_from_saved;
  __managed__ double nu_rfcut;
  __managed__ int n_lte_timesteps;
  __managed__ double cell_is_optically_thick;
  __managed__ int n_grey_timesteps;
  __managed__ int n_titer;
  __managed__ bool initial_iteration;
  __managed__ int max_bf_continua;
  __managed__ int n_kpktdiffusion_timesteps;
  __managed__ float kpktdiffusion_timescale;

}