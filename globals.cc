#include "globals.h"

#include <atomic>
#include <ctime>
#include <memory>
#include <vector>

#include "artisoptions.h"

#ifdef MPI_ON
#include <mpi.h>
#endif

namespace globals {

std::array<double, 3> syn_dir{};  // vector pointing from origin to observer

std::unique_ptr<struct time[]> timesteps = nullptr;

double *rpkt_emiss = nullptr;  /// Volume estimator for the rpkt emissivity

// for USE_LUT_PHOTOION = true
double *corrphotoionrenorm = nullptr;
double *gammaestimator = nullptr;

// for USE_LUT_BFHEATING = true
double *bfheatingestimator = nullptr;

double *ffheatingestimator = nullptr;
double *colheatingestimator = nullptr;
#ifdef DO_TITER
double *gammaestimator_save = nullptr;
double *bfheatingestimator_save = nullptr;
double *ffheatingestimator_save = nullptr;
double *colheatingestimator_save = nullptr;
#endif

int *ecounter = nullptr;
int *acounter = nullptr;

int nprocs_exspec = 1;
bool do_emission_res = true;

std::unique_ptr<bool[]> startofline;

double gamma_kappagrey;  // set to -ve for proper treatment. If possitive, then
                         // gamma_rays are treated as grey with this opacity.

double max_path_step;

int opacity_case;  // 0 normally, 1 for Fe-grp dependence.
                   /// MK: 2 for Fe-grp dependence and proportional to 1/rho
                   /// MK: 3 combination of 1 & 2 depending on a rho_crit
                   /// MK: 4 non-grey treatment

/// ATOMIC DATA

int nlines = -1;
std::vector<struct elementlist_entry> elements;
const struct linelist_entry *linelist = nullptr;
struct bflist_t *bflist = nullptr;

// for USE_LUT_BFHEATING = true
double *bfheating_coeff = nullptr;

struct rpkt_continuum_absorptioncoeffs *chi_rpkt_cont = nullptr;

/// Coolinglist
int ncoolingterms;

/// PHIXSLIST

double *allcont_nu_edge = nullptr;
const struct fullphixslist *allcont = nullptr;

// for either USE_LUT_PHOTOION = true or !USE_LUT_BFHEATING = false
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
std::atomic<int> nesc = 0;  // number of packets that escape during current timestep

double vmax;
double rmax;        /// Total mass and outer velocity/radius
double tmax = -1.;  /// End time of current simulation
double tmin = -1.;  /// Start time of current simulation

int ntimesteps = -1;        /// Number of timesteps
int timestep_initial = -1;  /// Initial timestep's number
int timestep_finish = -1;   /// Final timestep's number
int timestep = -1;          /// Current time step

/// New variables for other opacity cases, still grey.
double opcase3_normal;  /// MK: normalisation factor for opacity_case 3
double rho_crit_para;   /// MK: free parameter for the selection of the critical opacity in opacity_case 3
double rho_crit;        /// MK: critical opacity in opacity_case 3 (could now be declared locally)

int total_nlte_levels;  /// total number of nlte levels

bool simulation_continued_from_saved;
double nu_rfcut;
int num_lte_timesteps;
double cell_is_optically_thick;
int num_grey_timesteps;
int n_titer;
bool lte_iteration;
int n_kpktdiffusion_timesteps;
float kpktdiffusion_timescale;

void setup_mpi_vars() {
#ifdef MPI_ON
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
}

}  // namespace globals