#include "globals.h"

#include "sn3d.h"
#ifdef MPI_ON
#include <mpi.h>
#endif

namespace globals {

#if CUDA_ENABLED
__managed__ curandState curandstates[MCUDATHREADS];
#endif

__managed__ double syn_dir[3];  // vector pointing from origin to observer

// #define NRAYS_SYN 1 // number of rays traced in a syn calculation

// RAY rays[NRAYS_SYN];
__managed__ struct time *time_step = nullptr;

__managed__ int nsyn_time;
__managed__ double time_syn[MSYN_TIME];

__managed__ int emiss_offset;  // the index in the line list of the 1st line for which
                               // an emissivity estimator is recorded
__managed__ int emiss_max;     // actual number of frequency points in emissivity grid

/// THESE ARE THE GRID BASED ESTIMATORS
__managed__ float *compton_emiss = nullptr;  /// Volume estimator for the compton emissivity
__managed__ double *rpkt_emiss = nullptr;    /// Volume estimator for the rpkt emissivity

#if (!NO_LUT_PHOTOION)
__managed__ double *corrphotoionrenorm = nullptr;
__managed__ double *gammaestimator = nullptr;
#endif
#if (!NO_LUT_BFHEATING)
__managed__ double *bfheatingestimator = nullptr;
#endif
__managed__ double *ffheatingestimator = nullptr;
__managed__ double *colheatingestimator = nullptr;
#ifndef FORCE_LTE
#ifdef DO_TITER
__managed__ double *gammaestimator_save = nullptr;
__managed__ double *bfheatingestimator_save = nullptr;
__managed__ double *ffheatingestimator_save = nullptr;
__managed__ double *colheatingestimator_save = nullptr;
#endif
#endif

__managed__ int *ecounter = nullptr;
__managed__ int *acounter = nullptr;

__managed__ int nprocs_exspec = 1;
__managed__ bool do_emission_res = 1;

__managed__ bool file_set;  // 1 if the output files already exist. 0 otherwise.
__managed__ std::unique_ptr<bool[]> startofline;
__managed__ bool do_comp_est;  // 1 = compute compton emissivity estimators. 0 = don't
__managed__ bool do_r_lc;      // If not set to 1 then the opacity for r-packets is 0.
__managed__ int do_rlc_est;    // 1 = compute estimators for the r-pkt light curve.
                               // 2 = compute estimators with opacity weights
                               // 3 = compute estimators, but use only for gamma-heating rate

__managed__ double gamma_grey;  // set to -ve for proper treatment. If possitive, then
                                // gamma_rays are treated as grey with this opacity.

__managed__ double max_path_step;

__managed__ int opacity_case;  // 0 normally, 1 for Fe-grp dependence.
                               /// MK: 2 for Fe-grp dependence and proportional to 1/rho
                               /// MK: 3 combination of 1 & 2 depending on a rho_crit
                               /// MK: 4 non-grey treatment

/// ATOMIC DATA

__managed__ int nlines = -1;
__managed__ struct elementlist_entry *elements = nullptr;
__managed__ const struct linelist_entry *linelist = nullptr;
__managed__ struct bflist_t *bflist = nullptr;
__managed__ double *spontrecombcoeff = nullptr;
#if (!NO_LUT_PHOTOION)
__managed__ double *corrphotoioncoeff = nullptr;
#endif
#if (!NO_LUT_BFHEATING)
__managed__ double *bfheating_coeff = nullptr;
#endif
__managed__ double *bfcooling_coeff = nullptr;
;

__managed__ struct rpkt_cont_opacity *kappa_rpkt_cont = nullptr;

/// Coolinglist
__managed__ int ncoolingterms;

/// PHIXSLIST

__managed__ double *allcont_nu_edge = nullptr;
__managed__ const struct fullphixslist *allcont = nullptr;
#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
__managed__ struct groundphixslist *groundcont = nullptr;
#endif
__managed__ struct phixslist *phixslist = nullptr;
__managed__ int nbfcontinua = -1;
__managed__ int nbfcontinua_ground = -1;  /// number of bf-continua
__managed__ int NPHIXSPOINTS = -1;
__managed__ double NPHIXSNUINCREMENT = -1;

__managed__ struct cellhistory *cellhistory = nullptr;

#ifdef MPI_ON
MPI_Comm mpi_comm_node = MPI_COMM_NULL;
MPI_Comm mpi_comm_internode = MPI_COMM_NULL;
#endif

__managed__ int nprocs = -1;       // number of MPI processes
__managed__ int rank_global = -1;  // rank of the active MPI process

__managed__ int node_nprocs = -1;   // number of MPI processes on this node
__managed__ int rank_in_node = -1;  // local rank within this node

__managed__ int node_count = -1;  // number of MPI nodes
__managed__ int node_id = -1;     // unique number for each node

__managed__ int npkts = -1;
__managed__ int nesc = 0;  // number of packets that escape during current timestep

__managed__ double coordmax[3];
__managed__ double vmax;
__managed__ double rmax;        /// Total mass and outer velocity/radius
__managed__ double tmax = -1.;  /// End time of current simulation
__managed__ double tmin = -1.;  /// Start time of current simulation

__managed__ int ntstep = -1;      /// Number of timesteps
__managed__ int itstep = -1;      /// Initial timestep's number
__managed__ int ftstep = -1;      /// Final timestep's number
__managed__ int nts_global = -1;  /// Current time step

__managed__ double nusyn_min;
__managed__ double nusyn_max;  // limits on range for syn
__managed__ int nfake_gam;     // # of fake gamma ray lines for syn

/// New variables for other opacity cases, still grey.
__managed__ double opcase3_normal;  /// MK: normalisation factor for opacity_case 3
__managed__ double rho_crit_para;   /// MK: free parameter for the selection of the critical opacity in opacity_case 3
__managed__ double rho_crit;        /// MK: critical opacity in opacity_case 3 (could now be declared locally)

/// New variables for the non-grey case
__managed__ int debug_packet;  /// activate debug output for this packet if non negative

__managed__ int total_nlte_levels;  /// total number of nlte levels

__managed__ bool homogeneous_abundances;

__managed__ bool simulation_continued_from_saved;
__managed__ double nu_rfcut;
__managed__ int num_lte_timesteps;
__managed__ double cell_is_optically_thick;
__managed__ int num_grey_timesteps;
__managed__ int n_titer;
__managed__ bool initial_iteration;
__managed__ int n_kpktdiffusion_timesteps;
__managed__ float kpktdiffusion_timescale;

}  // namespace globals