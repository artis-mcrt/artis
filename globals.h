#include <limits>
#ifndef GLOBALS_H
#define GLOBALS_H

#include <mpi.h>

#include <array>
#include <cstddef>
#include <deque>
#include <mutex>
#include <vector>

#include "artisoptions.h"

struct TimeStep {
  double start{0.};                   // time at start of this timestep. [s]
  double width{0.};                   // Width of timestep. [s]
  double mid{0.};                     // Mid time in step - computed logarithmically. [s]
  double gamma_dep{0.};               // cmf gamma ray energy deposition from packet trajectories [erg]
  double gamma_dep_discrete{0.};      // cmf gamma ray energy deposition from absorption events [erg]
  double positron_dep{0.};            // cmf positron energy deposition from packet trajectories [erg]
  double positron_dep_discrete{0.};   // cmf positron energy deposition from absorption events [erg]
  double positron_emission{0.};       // cmf positron KE energy generation [erg]
  double eps_positron_ana_power{0.};  // cmf positron KE energy generation rate analytical [erg/s]
  double electron_dep{0.};            // cmf electron energy deposition from packet trajectories [erg]
  double electron_dep_discrete{0.};   // cmf electron energy deposition from absorption events [erg]
  double electron_emission{0.};       // cmf electron KE energy generation [erg]
  double eps_electron_ana_power{0.};  // cmf electron KE energy generation rate analytical [erg/s]
  double alpha_dep{0.};               // cmf alpha energy deposition from packet trajectories [erg]
  double alpha_dep_discrete{0.};      // cmf alpha energy deposition from absorption events [erg]
  double alpha_emission{0.};          // cmf alpha KE energy generation [erg]
  double eps_alpha_ana_power{0.};     // cmf alpha KE energy generation rate analytical [erg/s]
  double gamma_emission{0.};          // gamma decay energy generation in this timestep [erg]
  double qdot_betaminus{0.};          // energy generation from beta-minus decays (including neutrinos) [erg/s/g]
  double qdot_alpha{0.};              // energy generation from alpha decays (including neutrinos) [erg/s/g]
  double qdot_total{0.};              // energy generation from all decays (including neutrinos) [erg/s/g]
  double cmf_lum{0.};                 // cmf luminosity light curve [erg]
  int pellet_decays{0};               // Number of pellets that decay in this time step.
};

struct BFListEntry {
  int elementindex;
  int ionindex;
  int levelindex;
  int phixstargetindex;
};

struct FullPhotoionTransition {
  double nu_edge;
  int element;
  int ion;
  int level;
  int phixstargetindex;
  int upperlevel;
  const float *photoion_xs;
  double probability;
  int index_in_groundphixslist;
  int bfestimindex;
};

struct GroundPhotoion {
  double nu_edge;
  int element;
  int ion;
};

struct LevelTransition {
  int lineindex;
  int targetlevelindex;
  float einstein_A;
  float coll_str;
  float osc_strength;
  bool forbidden;
};

struct PhotoionTarget {
  double probability;  // fraction of phixs cross section leading to this final level
  int levelindex;      // index of upper ion level after photoionisation
};

struct EnergyLevel {
  double epsilon{-1};        // Excitation energy of this level relative to the neutral ground level.
  int alltrans_startdown{};  // index into globals::alltrans for first down transition from this level
  int ndowntrans{0};         // Number of down transitions from this level
  int nuptrans{0};           // Number of up transitions to this level
  int phixsstart{-1};        // index to start of photoionisation cross-sections table in global::allphixs
  int nphixstargets{0};      // number of target levels for photoionisation
  float stat_weight{0.};     // statistical weight of this level
  int phixstargetstart{};    // index into globals::allphixstargets
  int cont_index{-1};        // index of the bound-free continuum (for first target) sorted by
                             // element/ion/level/phixstargetindex
                             // (not an index into the nu_edge-sorted allcont list!)
  int closestgroundlevelcont{-1};
};

struct Ion {
  EnergyLevel *levels;      // Carries information for each level: 0,1,...,nlevels-1
  int ionstage;             // Which ionisation stage: XI=0, XII=1, XIII=2, ...
  int nlevels;              // Number of levels for this ionisation stage
  int nlevels_nlte;         // number of nlte levels for this ion
  int first_nlte;           // index into nlte_pops array of a grid cell
  int ionisinglevels;       // Number of levels which have a bf-continuum
  int maxrecombininglevel;  // level index of the highest level with a non-zero recombination rate
  int nlevels_groundterm;
  int coolingoffset;
  int ncoolingterms;
  int uniquelevelindexstart;
  int groundcontindex;
  double ionpot;  // Ionisation threshold to the next ionstage
};

struct Element {
  Ion *ions{};                         // Carries information for each ion: 0,1,...,nions-1
  int nions{0};                        // Number of ions for the current element
  int anumber{-1};                     // Atomic number
  int uniqueionindexstart{-1};         /// uniqueionindex index of the lowest ionisation stage of this element
  float initstablemeannucmass = {0.};  // Atomic mass number in multiple of MH
  bool has_nlte_levels{false};
};

struct TransitionLine {
  double nu;  // Frequency of the line transition
  float einstein_A;
  int elementindex;     // It's a transition of element (not its atomic number,
                        // but the (x-1)th element included in the simulation.
  int ionindex;         // The same for the elements ion
  int upperlevelindex;  // And the participating upper
  int lowerlevelindex;  // and lower levels
};

struct GSLIntegrationParas {
  double nu_edge;
  float T;
  const float *photoion_xs;
};

template <bool separatestimrecomb>
struct chphixstargets {
  double corrphotoioncoeff;
};

template <>
struct chphixstargets<true> {
  double corrphotoioncoeff;
  double separatestimrecomb;
};

using CellCachePhixsTargets = chphixstargets<SEPARATE_STIMRECOMB>;

enum ma_action {
  // Radiative deexcitation rate from this level.
  MA_ACTION_RADDEEXC = 0,
  // Collisional deexcitation rate from this level.
  MA_ACTION_COLDEEXC = 1,
  // Radiative recombination from this level.
  MA_ACTION_RADRECOMB = 2,
  // Collisional recombination rate from this level.
  MA_ACTION_COLRECOMB = 3,
  // Rate for internal downward transitions to same ionisation stage.
  MA_ACTION_INTERNALDOWNSAME = 4,
  // Rate for internal upward transitions to same ionisation stage.
  MA_ACTION_INTERNALDOWNLOWER = 5,
  // Rate for internal downward transitions to lower ionisation stage.
  MA_ACTION_INTERNALUPSAME = 6,
  // Rate for internal upward transitions to higher ionisation stage.
  MA_ACTION_INTERNALUPHIGHER = 7,
  // Rate for internal upward transitions to higher ionisation stage due to non-thermal collisions.
  MA_ACTION_INTERNALUPHIGHERNT = 8,
  MA_ACTION_COUNT = 9,
};

struct CellCacheLevels {
  std::array<double, MA_ACTION_COUNT> processrates;
  CellCachePhixsTargets *chphixstargets;
  double population;
  double *sum_epstrans_rad_deexc;
  double *sum_internal_down_same;
  double *sum_internal_up_same;
};

struct CellCacheIons {
  CellCacheLevels *chlevels;  // Pointer to the ions levellist.
};

struct CellCacheElements {
  CellCacheIons *chions;  // Pointer to the elements ionlist.
};

struct CellCache {
  double *cooling_contrib{};  // Cooling contributions by the different processes.
  CellCacheElements *chelements{};
  std::vector<CellCacheLevels> ch_all_levels;
  std::vector<double> ch_allcont_departureratios;
  std::vector<double> ch_allcont_nnlevel;
  std::vector<bool> ch_keep_this_cont;
  double chi_ff_nnionpart{-1};
  int nonemptymgi{-1};  // Identifies the cell the data is valid for.
};

namespace globals {

inline std::array<double, 3> syn_dir{};  // vector pointing from origin to observer

inline std::vector<TimeStep> timesteps;

// deposition estimators index by non-empty modelgridindex
// after normalisation factor has been applied, units will be erg/s/cm3
inline std::vector<double> dep_estimator_gamma;
inline std::vector<double> dep_estimator_positron;
inline std::vector<double> dep_estimator_electron;
inline std::vector<double> dep_estimator_alpha;

inline int bfestimcount{0};

// for USE_LUT_PHOTOION = true
inline double *corrphotoionrenorm{};
inline MPI_Win win_corrphotoionrenorm{MPI_WIN_NULL};

inline std::vector<double> gammaestimator;

// for USE_LUT_BFHEATING = true
inline std::vector<double> bfheatingestimator{};

inline std::vector<double> ffheatingestimator{};
inline std::vector<double> colheatingestimator{};
#ifdef DO_TITER
inline std::vector<double> gammaestimator_save{};
inline std::vector<double> bfheatingestimator_save{};
inline std::vector<double> ffheatingestimator_save{};
inline std::vector<double> colheatingestimator_save{};
#endif

inline std::vector<int> ecounter{};
inline std::vector<int> acounter{};

inline int nprocs_exspec{1};
inline bool do_emission_res{true};

inline double gamma_kappagrey{};  // set to -ve for proper treatment. If positive, then
                                  // gamma_rays are treated as grey with this opacity.

constexpr double GREY_OP = 0.1;

inline double max_path_step;

inline int opacity_case{};  // 0 grey, 1 for Fe-grp dependence.
                            // MK: 2 for Fe-grp dependence and proportional to 1/rho
                            // MK: 3 combination of 1 & 2 depending on a rho_crit
                            // MK: 4 non-grey treatment

// ATOMIC DATA

inline std::vector<float> ion_alpha_sp;  // alpha_sp for each ion and temperature table value

inline float *allphixs{};
inline LevelTransition *alltrans;
inline std::vector<PhotoionTarget> allphixstargets;

inline std::vector<Element> elements;

inline int nlines{-1};
inline const TransitionLine *linelist{};
inline std::vector<BFListEntry> bflist;

inline double *bfheating_coeff{};  // for USE_LUT_BFHEATING = true

inline std::vector<double> bfestim_nu_edge;
inline std::vector<double> allcont_nu_edge;
inline const FullPhotoionTransition *allcont{};

// for either USE_LUT_PHOTOION = true or !USE_LUT_BFHEATING = false
inline std::vector<GroundPhotoion> groundcont{};

inline int nbfcontinua{-1};         // number of bf-continua
inline int nbfcontinua_ground{-1};  // number of bf-continua from ground levels

inline int NPHIXSPOINTS{-1};
inline double NPHIXSNUINCREMENT{-1};

inline std::vector<CellCache> cellcache{};

inline MPI_Comm mpi_comm_node{MPI_COMM_NULL};
inline MPI_Comm mpi_comm_internode{MPI_COMM_NULL};

inline int nprocs{-1};
inline int my_rank{-1};

inline int node_nprocs{-1};
inline int rank_in_node{-1};

inline int node_count{-1};
inline int node_id{-1};

inline constexpr int npkts = MPKTS;
inline int nesc{0};

inline double vmax;
inline double rmax;
inline double tmax{-1};
inline double tmin{-1};

inline int ntimesteps{-1};
inline int timestep_initial{-1};
inline int timestep_finish{-1};
inline int timestep{-1};  // Current time step during the simulation

inline double opcase3_normal;  // MK: normalisation factor for opacity_case 3
inline double rho_crit_para;   // MK: free parameter for the selection of the critical opacity in opacity_case 3
inline double rho_crit;        // MK: critical opacity in opacity_case 3 (could now be declared locally)

inline int total_nlte_levels;

inline bool simulation_continued_from_saved;
inline double nu_rfcut;
inline int num_lte_timesteps;
inline double cell_is_optically_thick;
inline int num_grey_timesteps;
inline int n_titer;
inline bool lte_iteration;

inline bool packet_setup_phase{false};

inline std::deque<std::mutex> mutex_cellcachemacroatom;

inline void setup_mpi_vars() {
  MPI_Comm_rank(MPI_COMM_WORLD, &globals::my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &globals::nprocs);

  // make an intra-node communicator (group ranks that can share memory)
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globals::my_rank, MPI_INFO_NULL, &globals::mpi_comm_node);

  // get the local rank within this node
  MPI_Comm_rank(globals::mpi_comm_node, &globals::rank_in_node);

  // get the number of ranks on the node
  MPI_Comm_size(globals::mpi_comm_node, &globals::node_nprocs);

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef MAX_NODE_SIZE
  if (MAX_NODE_SIZE > 0 && globals::node_nprocs > MAX_NODE_SIZE) {
    // limit the number of ranks that can share memory
    MPI_Comm_split(globals::mpi_comm_node, globals::rank_in_node / MAX_NODE_SIZE, globals::my_rank,
                   &globals::mpi_comm_node);

    MPI_Comm_rank(globals::mpi_comm_node, &globals::rank_in_node);
    MPI_Comm_size(globals::mpi_comm_node, &globals::node_nprocs);
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // make an inter-node communicator (using local rank as the key for group membership)
  MPI_Comm_split(MPI_COMM_WORLD, globals::rank_in_node, globals::my_rank, &globals::mpi_comm_internode);

  // take the node id from the local rank 0 (node master) and broadcast it
  if (globals::rank_in_node == 0) {
    MPI_Comm_rank(globals::mpi_comm_internode, &globals::node_id);
    MPI_Comm_size(globals::mpi_comm_internode, &globals::node_count);
  }

  MPI_Bcast(&globals::node_id, 1, MPI_INT, 0, globals::mpi_comm_node);
  MPI_Bcast(&globals::node_count, 1, MPI_INT, 0, globals::mpi_comm_node);
}

}  // namespace globals

#endif  // GLOBALS_H
