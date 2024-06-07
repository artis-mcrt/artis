#pragma once
#ifndef GLOBALS_H
#define GLOBALS_H

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <array>
#include <cstddef>
#include <deque>
#include <mutex>
#include <vector>

#include "artisoptions.h"

struct TimeStep {
  double start;                   // time at start of this timestep. [s]
  double width;                   // Width of timestep. [s]
  double mid;                     // Mid time in step - computed logarithmically. [s]
  double gamma_dep;               // cmf gamma ray energy deposition from absorption events [erg]
  double gamma_dep_pathint;       // cmf gamma ray energy deposition from packet trajectories [erg]
  double positron_dep;            // cmf positron energy deposition [erg]
  double eps_positron_ana_power;  // cmf positron KE energy generation rate analytical [erg/s]
  double electron_dep;            // cmf electron energy deposition [erg]
  double electron_emission;       // cmf electron KE energy generation [erg]
  double eps_electron_ana_power;  // cmf electron KE energy generation rate analytical [erg/s]
  double alpha_dep;               // cmf alpha energy deposition [erg]
  double alpha_emission;          // cmf alpha KE energy generation [erg]
  double eps_alpha_ana_power;     // cmf alpha KE energy generation rate analytical [erg/s]
  double gamma_emission;          // gamma decay energy generation in this timestep [erg]
  double qdot_betaminus;          // energy generation from beta-minus decays (including neutrinos) [erg/s/g]
  double qdot_alpha;              // energy generation from alpha decays (including neutrinos) [erg/s/g]
  double qdot_total;              // energy generation from all decays (including neutrinos) [erg/s/g]
  double cmf_lum;                 // cmf luminosity light curve [erg]
  int pellet_decays;              // Number of pellets that decay in this time step.
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
  float *photoion_xs;
  double probability;
  int index_in_groundphixslist;
  int bfestimindex;
};

struct GroundPhotoion {
  double nu_edge;
  int element;
  int ion;
};

struct PhotoionTarget {
  double probability;  // fraction of phixs cross section leading to this final level
  int levelindex;      // index of upper ion level after photoionisation
};

struct LevelTransition {
  int lineindex;
  int targetlevelindex;
  float einstein_A;
  float coll_str;
  float osc_strength;
  bool forbidden;
};

struct EnergyLevel {
  double epsilon{-1};            /// Excitation energy of this level relative to the neutral ground level.
  LevelTransition *uptrans{};    /// Allowed upward transitions from this level
  LevelTransition *downtrans{};  /// Allowed downward transitions from this level
  int nuptrans{0};
  int ndowntrans{0};
  PhotoionTarget *phixstargets{};  /// pointer to table of target states and probabilities
  float *photoion_xs{};  /// Pointer to a lookup-table providing photoionisation cross-sections for this level.
  int nphixstargets{0};  /// length of phixstargets array:
  float stat_weight{0};  /// Statistical weight of this level.

  int cont_index{-1};  /// Index of the continuum associated to this level. Negative number.
  int closestgroundlevelcont{-1};
  bool metastable{};  ///
};

struct Ion {
  EnergyLevel *levels;      /// Carries information for each level: 0,1,...,nlevels-1
  int ionstage;             /// Which ionisation stage: XI=0, XII=1, XIII=2, ...
  int nlevels;              /// Number of levels for this ionisation stage
  int nlevels_nlte;         /// number of nlte levels for this ion
  int first_nlte;           /// index into nlte_pops array of a grid cell
  int ionisinglevels;       /// Number of levels which have a bf-continuum
  int maxrecombininglevel;  /// level index of the highest level with a non-zero recombination rate
  int nlevels_groundterm;
  int coolingoffset;
  int ncoolingterms;
  int uniquelevelindexstart;
  int groundcontindex;
  float *Alpha_sp;
  double ionpot;  /// Ionisation threshold to the next ionstage
  // int nbfcontinua;
  // ionsphixslist *phixslist;
};

struct Element {
  Ion *ions{};      /// Carries information for each ion: 0,1,...,nions-1
  int nions{0};     /// Number of ions for the current element
  int anumber{-1};  /// Atomic number
  //  int uppermost_ion;                       /// Highest ionisation stage which has a decent population for a given
  //  cell
  /// Be aware that this must not be used outside of the update_grid routine
  /// and their daughters. Neither it will work with OpenMP threads.
  int uniqueionindexstart{-1};         /// Index of the lowest ionisation stage of this element
  float abundance{0.};                 ///
  float initstablemeannucmass = {0.};  /// Atomic mass number in multiple of MH
  bool has_nlte_levels{false};
};

struct TransitionLine {
  double nu;  /// Frequency of the line transition
  float einstein_A;
  int elementindex;     /// It's a transition of element (not its atomic number,
                        /// but the (x-1)th element included in the simulation.
  int ionindex;         /// The same for the elements ion
  int upperlevelindex;  /// And the participating upper
  int lowerlevelindex;  /// and lower levels
};

struct GSLIntegrationParas {
  double nu_edge;
  float T;
  float *photoion_xs;
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
  /// Radiative deexcitation rate from this level.
  MA_ACTION_RADDEEXC = 0,
  /// Collisional deexcitation rate from this level.
  MA_ACTION_COLDEEXC = 1,
  /// Radiative recombination from this level.
  MA_ACTION_RADRECOMB = 2,
  /// Collisional recombination rate from this level.
  MA_ACTION_COLRECOMB = 3,
  /// Rate for internal downward transitions to same ionisation stage.
  MA_ACTION_INTERNALDOWNSAME = 4,
  /// Rate for internal upward transitions to same ionisation stage.
  MA_ACTION_INTERNALDOWNLOWER = 5,
  /// Rate for internal downward transitions to lower ionisation stage.
  MA_ACTION_INTERNALUPSAME = 6,
  /// Rate for internal upward transitions to higher ionisation stage.
  MA_ACTION_INTERNALUPHIGHER = 7,
  /// Rate for internal upward transitions to higher ionisation stage due to non-thermal collisions.
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
  CellCacheLevels *chlevels;  /// Pointer to the ions levellist.
};

struct CellCacheElements {
  CellCacheIons *chions;  /// Pointer to the elements ionlist.
};

struct CellCache {
  double *cooling_contrib{};  /// Cooling contributions by the different processes.
  CellCacheElements *chelements{};
  CellCacheLevels *ch_all_levels{};
  double *ch_allcont_departureratios{};
  double chi_ff_nnionpart{-1};
  int cellnumber{-1};  /// Identifies the cell the data is valid for.
};

namespace globals {

inline std::array<double, 3> syn_dir{};  // vector pointing from origin to observer

inline std::vector<TimeStep> timesteps;

inline std::vector<double> dep_estimator_gamma;
inline std::vector<double> dep_estimator_positron;
inline std::vector<double> dep_estimator_electron;
inline std::vector<double> dep_estimator_alpha;

inline int bfestimcount{0};

// for USE_LUT_PHOTOION = true
inline double *corrphotoionrenorm{};
#ifdef MPI_ON
inline MPI_Win win_corrphotoionrenorm{MPI_WIN_NULL};
#endif

inline double *gammaestimator{};

// for USE_LUT_BFHEATING = true
inline double *bfheatingestimator{};

inline double *ffheatingestimator{};
inline double *colheatingestimator{};
#ifdef DO_TITER
inline double *gammaestimator_save{};
inline double *bfheatingestimator_save{};
inline double *ffheatingestimator_save{};
inline double *colheatingestimator_save{};
#endif

inline int *ecounter{};
inline int *acounter{};

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

/// ATOMIC DATA

inline int nlines{-1};
inline std::vector<Element> elements;

inline const TransitionLine *linelist{};
inline std::vector<BFListEntry> bflist;

inline double *bfheating_coeff{};  // for USE_LUT_BFHEATING = true

inline std::vector<double> bfestim_nu_edge;
inline std::vector<double> allcont_nu_edge;
inline const FullPhotoionTransition *allcont{};

// for either USE_LUT_PHOTOION = true or !USE_LUT_BFHEATING = false
inline GroundPhotoion *groundcont{};

inline int nbfcontinua{-1};         // number of bf-continua
inline int nbfcontinua_ground{-1};  // number of bf-continua from ground levels

inline int NPHIXSPOINTS{-1};
inline double NPHIXSNUINCREMENT{-1};

inline CellCache *cellcache{};

#ifdef MPI_ON
inline MPI_Comm mpi_comm_node{MPI_COMM_NULL};
inline MPI_Comm mpi_comm_internode{MPI_COMM_NULL};
#endif

inline int nprocs{-1};
inline int rank_global{-1};

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

inline std::deque<std::mutex> mutex_cellcachemacroatom;

inline void setup_mpi_vars() {
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

#ifdef MAX_NODE_SIZE
  if (MAX_NODE_SIZE > 0 && globals::node_nprocs > MAX_NODE_SIZE) {
    // limit the number of ranks that can share memory
    MPI_Comm_split(globals::mpi_comm_node, globals::rank_in_node / MAX_NODE_SIZE, globals::rank_global,
                   &globals::mpi_comm_node);

    MPI_Comm_rank(globals::mpi_comm_node, &globals::rank_in_node);
    MPI_Comm_size(globals::mpi_comm_node, &globals::node_nprocs);
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

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

#endif  // GLOBALS_H
