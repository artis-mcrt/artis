#include "globals.h"

#include "sn3d.h"
#ifdef MPI_ON
#include <mpi.h>
#endif

namespace globals {

double syn_dir[3];  // vector pointing from origin to observer

// #define NRAYS_SYN 1 // number of rays traced in a syn calculation

// RAY rays[NRAYS_SYN];
struct time *time_step = nullptr;

int nsyn_time;
double time_syn[MSYN_TIME];

int emiss_offset;  // the index in the line list of the 1st line for which
                   // an emissivity estimator is recorded
int emiss_max;     // actual number of frequency points in emissivity grid

/// THESE ARE THE GRID BASED ESTIMATORS
float *compton_emiss = nullptr;  /// Volume estimator for the compton emissivity
double *rpkt_emiss = nullptr;    /// Volume estimator for the rpkt emissivity

// for NO_LUT_PHOTOION = false
double *corrphotoionrenorm = nullptr;
double *gammaestimator = nullptr;

// for NO_LUT_BFHEATING = false
double *bfheatingestimator = nullptr;

double *ffheatingestimator = nullptr;
double *colheatingestimator = nullptr;
#ifndef FORCE_LTE
#ifdef DO_TITER
double *gammaestimator_save = nullptr;
double *bfheatingestimator_save = nullptr;
double *ffheatingestimator_save = nullptr;
double *colheatingestimator_save = nullptr;
#endif
#endif

int *ecounter = nullptr;
int *acounter = nullptr;

int nprocs_exspec = 1;
bool do_emission_res = 1;

bool file_set;  // 1 if the output files already exist. 0 otherwise.
std::unique_ptr<bool[]> startofline;
bool do_comp_est;  // 1 = compute compton emissivity estimators. 0 = don't
bool do_r_lc;      // If not set to 1 then the opacity for r-packets is 0.
int do_rlc_est;    // 1 = compute estimators for the r-pkt light curve.
                   // 2 = compute estimators with opacity weights
                   // 3 = compute estimators, but use only for gamma-heating rate

double gamma_grey;  // set to -ve for proper treatment. If possitive, then
                    // gamma_rays are treated as grey with this opacity.

double max_path_step;

int opacity_case;  // 0 normally, 1 for Fe-grp dependence.
                   /// MK: 2 for Fe-grp dependence and proportional to 1/rho
                   /// MK: 3 combination of 1 & 2 depending on a rho_crit
                   /// MK: 4 non-grey treatment

/// ATOMIC DATA

int nlines = -1;
struct elementlist_entry *elements = nullptr;
const struct linelist_entry *linelist = nullptr;
struct bflist_t *bflist = nullptr;
double *spontrecombcoeff = nullptr;

// for NO_LUT_PHOTOION = false
double *corrphotoioncoeff = nullptr;

// for NO_LUT_BFHEATING = false
double *bfheating_coeff = nullptr;

double *bfcooling_coeff = nullptr;

struct rpkt_cont_opacity *kappa_rpkt_cont = nullptr;

/// Coolinglist
int ncoolingterms;

/// PHIXSLIST

double *allcont_nu_edge = nullptr;
const struct fullphixslist *allcont = nullptr;

// for either NO_LUT_PHOTOION = false or NO_LUT_BFHEATING = false
struct groundphixslist *groundcont = nullptr;

struct phixslist *phixslist = nullptr;
int nbfcontinua = -1;
int nbfcontinua_ground = -1;  /// number of bf-continua
int NPHIXSPOINTS = -1;
double NPHIXSNUINCREMENT = -1;

struct cellhistory *cellhistory = nullptr;

#ifdef MPI_ON
MPI_Comm mpi_comm_node = MPI_COMM_NULL;
MPI_Comm mpi_comm_internode = MPI_COMM_NULL;
#endif

int nprocs = -1;       // number of MPI processes
int rank_global = -1;  // rank of the active MPI process

int node_nprocs = -1;   // number of MPI processes on this node
int rank_in_node = -1;  // local rank within this node

int node_count = -1;  // number of MPI nodes
int node_id = -1;     // unique number for each node

constexpr int npkts = MPKTS;
int nesc = 0;  // number of packets that escape during current timestep

double coordmax[3];
double vmax;
double rmax;        /// Total mass and outer velocity/radius
double tmax = -1.;  /// End time of current simulation
double tmin = -1.;  /// Start time of current simulation

int ntstep = -1;      /// Number of timesteps
int itstep = -1;      /// Initial timestep's number
int ftstep = -1;      /// Final timestep's number
int nts_global = -1;  /// Current time step

double nusyn_min;
double nusyn_max;  // limits on range for syn
int nfake_gam;     // # of fake gamma ray lines for syn

/// New variables for other opacity cases, still grey.
double opcase3_normal;  /// MK: normalisation factor for opacity_case 3
double rho_crit_para;   /// MK: free parameter for the selection of the critical opacity in opacity_case 3
double rho_crit;        /// MK: critical opacity in opacity_case 3 (could now be declared locally)

int total_nlte_levels;  /// total number of nlte levels

bool homogeneous_abundances;

bool simulation_continued_from_saved;
double nu_rfcut;
int num_lte_timesteps;
double cell_is_optically_thick;
int num_grey_timesteps;
int n_titer;
bool initial_iteration;
int n_kpktdiffusion_timesteps;
float kpktdiffusion_timescale;

}  // namespace globals