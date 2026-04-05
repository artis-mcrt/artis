#ifndef GLOBALS_H
#define GLOBALS_H

#ifdef STDPAR_ON
#include <thread>
#endif
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

#include "mpi_logging.h"

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

struct Ion {
  int nlevels{0};  // Number of levels for this ionisation stage
  int nlevels_excited_nlte{0};  // number of nlte levels for this ion
  int allnltelevelsindexstart{-1};  // index into nlte_pops array for first excited nlte level
  int nlevels_ionising{0};  // Number of levels which have a bf-continuum
  int maxrecombininglevel{-1};  // level index of the highest level with a non-zero recombination rate
  int nlevels_autoion{0};  // Number of levels that can autoionise
  int nlevels_groundterm{0};
  int coolingoffset{-1};
  int ncoolingterms{0};
  int uniquelevelindexstart{-1};  // index of the first level in the alllevels list
  int groundcontindex{-1};
  double ionpot{NAN};  // Ionisation threshold to the next ionstage
};

struct Element {
  std::span<Ion> ions;  // subspan of the allions array for this element
  int anumber{-1};  // Atomic number
  int lowest_ionstage{-1};  // ionisation stage (charge - 1) of ion 0 for this element
  int uniqueionindexstart{-1};  /// uniqueionindex index of the lowest ionisation stage of this element
  float initstablemeannucmass = {0.};  // Atomic mass number in multiple of MH
  bool has_nlte_levels{false};
};

namespace globals {

struct TimeStep {
  double start{0.};  // time at start of this timestep. [s]
  double width{0.};  // Width of timestep. [s]
  double mid{0.};  // Mid time in step - computed logarithmically. [s]
  double gamma_dep{0.};  // cmf gamma ray energy deposition from packet trajectories [erg]
  double gamma_dep_discrete{0.};  // cmf gamma ray energy deposition from absorption events [erg]
  double positron_dep{0.};  // cmf positron energy deposition from packet trajectories [erg]
  double positron_dep_discrete{0.};  // cmf positron energy deposition from absorption events [erg]
  double positron_emission{0.};  // cmf positron KE energy generation [erg]
  double eps_positron_ana_power{0.};  // cmf positron KE energy generation rate analytical [erg/s]
  double electron_dep{0.};  // cmf electron energy deposition from packet trajectories [erg]
  double electron_dep_discrete{0.};  // cmf electron energy deposition from absorption events [erg]
  double electron_emission{0.};  // cmf electron KE energy generation [erg]
  double eps_electron_ana_power{0.};  // cmf electron KE energy generation rate analytical [erg/s]
  double alpha_dep{0.};  // cmf alpha energy deposition from packet trajectories [erg]
  double alpha_dep_discrete{0.};  // cmf alpha energy deposition from absorption events [erg]
  double alpha_emission{0.};  // cmf alpha KE energy generation [erg]
  double eps_alpha_ana_power{0.};  // cmf alpha KE energy generation rate analytical [erg/s]
  double spfission_dep_discrete{0.};  // cmf spontaneous fission energy deposition from absorption events [erg]
  double eps_spfission_ana_power{0.};  // cmf spontaneous fission energy generation rate analytical [erg/s]
  double gamma_emission{0.};  // gamma decay energy generation in this timestep [erg]
  double qdot_betaminus{0.};  // energy generation from beta-minus decays (including neutrinos) [erg/s/g]
  double qdot_alpha{0.};  // energy generation from alpha decays (including neutrinos) [erg/s/g]
  double qdot_spfission{0.};  // energy generation from spontaneous fission decays (including neutrinos) [erg/s/g]
  double qdot_total{0.};  // energy generation from all decays (including neutrinos) [erg/s/g]
  int pellet_decays{0};  // Number of pellets that decay in this time step.
};
inline std::vector<TimeStep> timesteps;

// deposition estimators index by non-empty modelgridindex
// after normalisation factor has been applied, units will be erg/s/cm3
inline std::vector<double> dep_estimator_gamma;
inline std::vector<double> dep_estimator_positron;
inline std::vector<double> dep_estimator_electron;
inline std::vector<double> dep_estimator_alpha;

// for USE_LUT_PHOTOION = true
inline MPI_shared_array<double> corrphotoionrenorm{};

inline std::vector<double> gammaestimator;

// for USE_ION_BFHEATING_ESTIMATORS = true
inline std::vector<double> bfheatingestimator{};

inline std::vector<double> ffheatingestimator{};
inline std::vector<double> colheatingestimator{};
#ifdef DO_TITER
inline std::vector<double> gammaestimator_save{};
inline std::vector<double> bfheatingestimator_save{};
inline std::vector<double> ffheatingestimator_save{};
inline std::vector<double> colheatingestimator_save{};
#endif

inline int nprocs_exspec{1};

constexpr double GREY_OP = 0.1;

inline double max_path_step;

// ATOMIC DATA

inline MPI_shared_array<float> allphixs{};

struct AllTransitions {
  MPI_shared_array<const int> lineindex;
  MPI_shared_array<const int> targetlevelindex;
  MPI_shared_array<const float> einstein_A;
  MPI_shared_array<const float> coll_str;
  MPI_shared_array<const float> osc_strength;
  MPI_shared_array<const bool> forbidden;
};
inline AllTransitions alltrans;

struct LevelAutoion {
  float autoion_A;  // Autoionisation A-value
  int elementindex;  // index (not atomic number) for the element involved
  int lowerionindex;
  int lowerlevelindex;  // this will be for a level index of the lower ion
  int upperionindex;
  int upperlevelindex;  // this will be for a level index of the upper ion.
                        // Note: level of the lower ion should also be at higher energy than of the higher ion
};
inline std::span<LevelAutoion> allautoion;

inline MPI_shared_array<const int> allphixstargets_levelindex;  // index of upper ion level after photoionisation
inline MPI_shared_array<const double>
    allphixstargets_probability;  // fraction of phixs cross section leading to associated final level

struct AllLevels {
  // these arrays are indexed by uniquelevelindex, which can be derived from the element, ion, level

  MPI_shared_array<const double> epsilon;
  MPI_shared_array<const float> statweight;

  // index into globals::alltrans for first down transition from each level
  std::span<const int> alltrans_startdown;

  // Number of down transitions from each level
  std::span<const int> ndowntrans;

  // Number of up transitions from each level
  std::span<const int> nuptrans;

  // Number of autoionizing transition from this level
  std::span<int> nautoiondowntrans;

  // Number of di-el captures up from this level
  std::span<int> nautoionuptrans;

  // index into globals::allautoion for first autoion from this level
  std::span<int> allautoion_start;

  std::span<int> closestgroundlevelcont;

  // index to start of photoionisation cross-sections table in global::allphixs
  std::span<int> phixsstart;

  // number of target levels for photoionisation
  MPI_shared_array<int> nphixstargets;

  // index into globals::allphixstargets for the first target level
  std::span<int> phixstargetstart;

  // index of the bound-free continuum (for first target) sorted by element/ion/level/phixstargetindex (not an index
  // into the nu_edge-sorted allcont list!)
  std::span<int> bflist_start;

  // index into cellcache allmacroatomictransitions for each level. This is
  // different to the alltrans index because two types of down transitions are stored separately
  // per level as well as the up transitions
  MPI_shared_array<const int> matransblock_start;
};

inline AllLevels alllevels{};

inline std::vector<Element> elements;
inline MPI_shared_array<Ion> allions;

struct TransitionLines {
  MPI_shared_array<const double> nu;  // Frequency of the line transition
  std::span<const float> einstein_A;
  std::span<const int> elementindex;  // It's a transition of element (not its atomic number,
                                      // but the (x-1)th element included in the simulation.
  std::span<const int> ionindex;  // The same for the elements ion
  std::span<const int> upperlevelindex;  // And the participating upper
  std::span<const int> lowerlevelindex;  // and lower levels
};
inline TransitionLines linelist{};
inline int nlines{-1};

struct BFListEntry {
  int elementindex;
  int ionindex;
  int levelindex;
  int phixstargetindex;
};
// the bound-free list sorted by element/ion/level/phixstargetindex (not nu_edge)
inline std::vector<BFListEntry> bflist;

inline std::span<const double> bfestim_nu_edge{};

struct AllCont {
  std::span<const double> nu_edge;
  std::span<const int> element;
  std::span<const int> ion;
  std::span<const int> level;
  std::span<const int> phixstargetindex;
  std::span<const int> upperlevel;
  std::span<const int> uniquelevelindex;
  std::span<const double> probability;
  std::span<const int> index_in_groundphixslist;
  std::span<const int> bfestimindex;
};
inline AllCont allcont{};

// Used when USE_LUT_PHOTOION or USE_ION_BFHEATING_ESTIMATORS is enabled
inline std::span<const double> groundcont_nu_edge{};
inline std::span<const int> groundcont_element{};
inline std::span<const int> groundcont_ion{};

inline int nbfcontinua{-1};  // number of bf-continua
inline int nbfcontinua_ground{-1};  // number of bf-continua from ground levels

inline int NPHIXSPOINTS{-1};  // number of photoionisation cross-section points per level
inline double NPHIXSNUINCREMENT{-1};  // frequency increment between points as a fraction of nu_edge

struct CellCache {
  int nonemptymgi{-1};  // non-empty model grid index for this cache slot
  std::vector<double> cooling_contrib;  // Cooling contributions by the different processes.
  std::vector<double> alllevels_pops;
  std::vector<double> alllevels_maprocessrates;  // rates for macroatom processes
  std::vector<double> allmacroatomictransitions;  // cumulative macroatom transition rates for all levels
  std::vector<double> allcont_modified_departureratios;
  std::vector<double> allcont_nnlevel;
  std::vector<bool> allcont_keep;
  double chi_ff_nnionpart{-1};
  std::vector<double> allphixstargets_corrphotoioncoeff;
  std::vector<int> cooling_contrib_locks;
  std::vector<int> allmacroatomictransitions_locks;
};
inline std::vector<CellCache> cellcache{};

inline double vmax{NAN};
inline double rmax{NAN};
inline double tmax{-1};
inline double tmin{-1};

inline int ntimesteps{-1};
inline int timestep_initial{-1};
inline int timestep_finish{-1};
inline int timestep{-1};  // Current time step during the simulation

inline int total_nlte_levels{0};

inline bool simulation_continued_from_saved{false};
inline int num_lte_timesteps{-1};
inline double cell_is_optically_thick{NAN};
inline int num_grey_timesteps{-1};
inline int n_titer{1};
inline bool lte_iteration{false};

}  // namespace globals

constexpr int cellcacheslotid = 0;

[[nodiscard]] inline auto get_max_threads() -> int {
#ifdef STDPAR_ON
  return static_cast<int>(std::thread::hardware_concurrency());
#else
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
#endif
}

[[nodiscard]] inline auto get_thread_num() -> int {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}
#endif  // GLOBALS_H
