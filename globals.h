// Global simulation state: the core data structs (elements, ions, levels, transitions,
// timesteps, cell cache) and the globals namespace of shared arrays, estimators, and run
// parameters used throughout the code.

#ifndef GLOBALS_H
#define GLOBALS_H

#include "constants.h"
#ifdef STDPAR_ON
#include <thread>
#endif
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "mpi_logging.h"

// ma_action enumerates the types of macroatom transitions that can occur from a given level
enum ma_action {
  // Rate for radiative deexcitation
  MA_ACTION_RADDEEXC = 0,
  // Rate for thermal collisional deexcitation
  MA_ACTION_COLDEEXC = 1,
  // Rate for radiative recombination
  MA_ACTION_RADRECOMB = 2,
  // Rate for collisional recombination
  MA_ACTION_COLRECOMB = 3,
  // Rate for internal downward transitions to same ionisation stage
  MA_ACTION_INTERNALDOWNSAME = 4,
  // Rate for internal downward transitions to lower ionisation stage
  MA_ACTION_INTERNALDOWNLOWER = 5,
  // Rate for internal upward transitions to same ionisation stage
  MA_ACTION_INTERNALUPSAME = 6,
  // Rate for internal upward transitions to higher ionisation stage
  MA_ACTION_INTERNALUPHIGHER = 7,
  // Rate for internal upward transitions to any higher ionisation stage due to non-thermal collisions
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
  int lowest_ionstage{-1};  // ionisation stage (charge + 1) of ion 0 for this element
  int uniqueionindexstart{-1};  // uniqueionindex of the lowest ionisation stage of this element
  float initstablemeannucmass = {0.};  // mean nuclear mass [g] of the element's stable component (mass_amu * MH)
  bool has_nlte_levels{false};
};

namespace globals {

// The *_discrete, *_emission, and pellet_decays members are accumulated atomically by all threads during packet
// propagation, so they are cache-line aligned in CPU multithreaded modes to keep them from false sharing with each
// other and with the constantly-read start/width/mid members.
struct TimeStep {
  double start{0.};  // time at start of this timestep. [s]
  double width{0.};  // Width of timestep. [s]
  double mid{0.};  // Mid time in step - computed logarithmically. [s]
  double gamma_dep{0.};  // cmf gamma ray energy deposition from packet trajectories [erg]
  ALIGNAS_AVOID_FALSE_SHARING double gamma_dep_discrete{
      0.,
  };  // cmf gamma ray energy deposition from absorption events [erg]
  double positron_dep{0.};  // cmf positron energy deposition from packet trajectories [erg]
  ALIGNAS_AVOID_FALSE_SHARING double positron_dep_discrete{
      0.,
  };  // cmf positron energy deposition from absorption events [erg]
  ALIGNAS_AVOID_FALSE_SHARING double positron_emission{0.};  // cmf positron KE energy generation [erg]
  double eps_positron_ana_power{0.};  // cmf positron KE energy generation rate analytical [erg/s]
  double electron_dep{0.};  // cmf electron energy deposition from packet trajectories [erg]
  ALIGNAS_AVOID_FALSE_SHARING double electron_dep_discrete{
      0.,
  };  // cmf electron energy deposition from absorption events [erg]
  ALIGNAS_AVOID_FALSE_SHARING double electron_emission{0.};  // cmf electron KE energy generation [erg]
  double eps_electron_ana_power{0.};  // cmf electron KE energy generation rate analytical [erg/s]
  double alpha_dep{0.};  // cmf alpha energy deposition from packet trajectories [erg]
  ALIGNAS_AVOID_FALSE_SHARING double alpha_dep_discrete{
      0.,
  };  // cmf alpha energy deposition from absorption events [erg]
  ALIGNAS_AVOID_FALSE_SHARING double alpha_emission{0.};  // cmf alpha KE energy generation [erg]
  double eps_alpha_ana_power{0.};  // cmf alpha KE energy generation rate analytical [erg/s]
  ALIGNAS_AVOID_FALSE_SHARING double spfission_dep_discrete{
      0.,
  };  // cmf spontaneous fission energy deposition from absorption events [erg].
      // Unlike the other *_dep_discrete counters, this is accumulated at pellet decay in update_pellet()
      // rather than at a separate deposition event, because a spontaneous fission pellet converts straight
      // to TYPE_NTALPHA_FISPROD_DEPOSITED. Emission and deposition therefore coincide, which is why
      // write_deposition_file() uses this same member in both the total deposition and the total emission sums.
  double eps_spfission_ana_power{0.};  // cmf spontaneous fission energy generation rate analytical [erg/s]
  ALIGNAS_AVOID_FALSE_SHARING double gamma_emission{0.};  // gamma decay energy generation in this timestep [erg]
  double qdot_betaminus{0.};  // energy generation from beta-minus decays (including neutrinos) [erg/s/g]
  double qdot_alpha{0.};  // energy generation from alpha decays (including neutrinos) [erg/s/g]
  double qdot_spfission{0.};  // energy generation from spontaneous fission decays (including neutrinos) [erg/s/g]
  double qdot_total{0.};  // energy generation from all decays (including neutrinos) [erg/s/g]
  ALIGNAS_AVOID_FALSE_SHARING int pellet_decays{0};  // Number of pellets that decay in this time step.
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

inline double max_path_step{NAN};  // set at the end of update_grid(); NaN disables the cap until then

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
inline MPI_shared_array<LevelAutoion> allautoion;

inline MPI_shared_array<const int> allphixstargets_levelindex;  // index of upper ion level after photoionisation
inline MPI_shared_array<const double>
    allphixstargets_probability;  // fraction of phixs cross section leading to associated final level

struct AllLevels {
  // these arrays are indexed by uniquelevelindex, which can be derived from the element, ion, level

  MPI_shared_array<const double> epsilon;
  MPI_shared_array<const float> statweight;

  // index into globals::alltrans for first down transition from each level
  MPI_shared_array<const int> alltrans_startdown;

  // Number of down transitions from each level
  MPI_shared_array<const int> ndowntrans;

  // Number of up transitions from each level
  MPI_shared_array<const int> nuptrans;

  // Number of autoionizing transition from this level
  MPI_shared_array<int> nautoiondowntrans;

  // Number of di-el captures up from this level
  MPI_shared_array<int> nautoionuptrans;

  // index into globals::allautoion for first autoion from this level
  MPI_shared_array<int> allautoion_start;

  MPI_shared_array<int> closestgroundlevelcont;

  // index to start of photoionisation cross-sections table in global::allphixs
  MPI_shared_array<int> phixsstart;

  // number of target levels for photoionisation
  MPI_shared_array<int> nphixstargets;

  // index into globals::allphixstargets for the first target level
  MPI_shared_array<int> phixstargetstart;

  // index of the bound-free continuum (for first target) sorted by element/ion/level/phixstargetindex (not an index
  // into the nu_edge-sorted allcont list!)
  MPI_shared_array<int> bflist_start;

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
  MPI_shared_array<const int> elementindex;  // index into the list of elements included in the simulation
                                             // (not the atomic number)
  MPI_shared_array<const int> ionindex;  // The same for the elements ion
  MPI_shared_array<const int> uniquelevelindex_lower;  // globally unique index of the lower level of the transition
  MPI_shared_array<const int> uniquelevelindex_upper;  // globally unique index of the upper level of the transition
  MPI_shared_array<const float> B_ul;  // Einstein B coefficient for stimulated emission
  MPI_shared_array<const float> B_lu;  // Einstein B coefficient for absorption (B_lu = (g_u/g_l) * B_ul)
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

inline MPI_shared_array<const double> bfestim_nu_edge{};

struct AllCont {
  MPI_shared_array<const double> nu_edge;
  MPI_shared_array<const int> element;
  MPI_shared_array<const int> ion;
  MPI_shared_array<const int> level;
  MPI_shared_array<const int> phixstargetindex;
  MPI_shared_array<const int> upperlevel;
  MPI_shared_array<const int> uniquelevelindex;
  MPI_shared_array<const double> probability;
  // index into the ground-level continuum estimator arrays, or -1 for continua that do not feed them
  // (only a ground level's first photoionisation target does). This is the nearest-nu_edge slot found by
  // search_groundphixslist(), the same value closestgroundlevelcont holds, which is not necessarily the
  // ion's own get_groundcontindex() slot when two ions share a ground-state threshold.
  MPI_shared_array<const int> groundcontestimindex;
  MPI_shared_array<const int> bfestimindex;
};
inline AllCont allcont{};

// Used when USE_LUT_PHOTOION or USE_ION_BFHEATING_ESTIMATORS is enabled
inline MPI_shared_array<const double> groundcont_nu_edge{};
inline MPI_shared_array<const int> groundcont_element{};
inline MPI_shared_array<const int> groundcont_ion{};

inline int nbfcontinua{-1};  // number of bf-continua
inline int nbfcontinua_ground{-1};  // number of bf-continua from ground levels

inline int NPHIXSPOINTS{-1};  // number of photoionisation cross-section points per level
inline double NPHIXSNUINCREMENT{-1};  // frequency increment between points as a fraction of nu_edge

// A cell cache slot holds pre-calculated quantities for a single model grid cell. The large per-cell
// arrays are non-owning views (spans) into shared-memory backing storage held in
// globals::cellcache_backing, with each slot viewing its own sub-range.
struct CellCache {
  int nonemptymgi{-1};  // non-empty model grid index for this cache slot
  std::span<double> cooling_contrib;  // Cooling contributions by the different processes.
  std::span<double> alllevels_pops;
  std::span<double> alllevels_maprocessrates;  // rates for macroatom processes
  std::span<double> allmacroatomictransitions;  // cumulative macroatom transition rates for all levels
  std::span<double> allcont_modified_departureratios;
  // modified departure ratio times exp(H * nu_edge / (KB * T_e)), so that the stimulated correction
  // factor at any nu needs only a multiply by the per-call exp(-H * nu / (KB * T_e)). Negative means
  // not cached for this cell (either not evaluated yet, or the product would overflow), in which case
  // the factor is computed directly from nu - nu_edge (see calculate_chi_bf_gammacontr())
  std::span<double> allcont_stimfactor_edgepart;
  // level population of each bound-free continuum's lower level
  std::span<double> allcont_nnlevel;
  // whether each bound-free continuum contributes to this cell's opacity at all: a populated lower level
  // of a species the cell contains (see keep_this_cont). Bit (i % keepbits_wordbits) of word
  // (i / keepbits_wordbits), in ascending continuum order. This duplicates "allcont_nnlevel[i] > 0", but is
  // held separately
  // and packed because it is what calculate_chi_bf_gammacontr() scans across the whole window on every
  // evaluation: with only ~11% of continua contributing in a typical cell (classic; the
  // DETAILED_BF_ESTIMATORS_ON presets keep nearly all of them), packing lets that scan clear a whole word
  // of absent species per test.
  std::span<std::uint64_t> allcont_keepbits;
  std::span<double> chi_ff_nnionpart;  // single element per cell (stored as a span to allow shared backing)
  std::span<double> allphixstargets_corrphotoioncoeff;
  // The lazy mutex-guarded macroatom/cooling rate calculation is only used when cellcache_singleslot is
  // true. When it is false, these rates are pre-calculated for every cell (see cellcacheslot_populate) so
  // no locking is required during packet propagation and these vectors stay empty.
  std::vector<PaddedMutex> cooling_contrib_locks;
  std::vector<PaddedMutex> allmacroatomictransitions_locks;

  [[nodiscard]] auto get_mem_usage() const {
    auto mem_usage = (cooling_contrib.size() * sizeof(double));
    mem_usage += (alllevels_pops.size() * sizeof(double));
    mem_usage += (alllevels_maprocessrates.size() * sizeof(double));
    mem_usage += (allmacroatomictransitions.size() * sizeof(double));
    mem_usage += (allcont_modified_departureratios.size() * sizeof(double));
    mem_usage += (allcont_stimfactor_edgepart.size() * sizeof(double));
    mem_usage += (allcont_nnlevel.size() * sizeof(double));
    mem_usage += (allcont_keepbits.size() * sizeof(allcont_keepbits[0]));
    mem_usage += (chi_ff_nnionpart.size() * sizeof(double));
    mem_usage += (allphixstargets_corrphotoioncoeff.size() * sizeof(double));
    mem_usage += (cooling_contrib_locks.size() * sizeof(PaddedMutex));
    mem_usage += (allmacroatomictransitions_locks.size() * sizeof(PaddedMutex));
    return mem_usage;
  }
};
inline std::vector<CellCache> cellcache{};

// Backing storage for the cell cache arrays, allocated in node-shared memory. When cellcache_singleslot
// is false, each array spans every non-empty cell and cellcache[nonemptymgi] views the relevant sub-range,
// shared by all MPI ranks on the node. When it is true, each array holds one reusable slot per node rank
// and each rank uses only its own cellcache[rank_in_node] view.
struct CellCacheBacking {
  MPI_shared_array<double> cooling_contrib;
  MPI_shared_array<double> alllevels_pops;
  MPI_shared_array<double> alllevels_maprocessrates;
  MPI_shared_array<double> allmacroatomictransitions;
  MPI_shared_array<double> allcont_modified_departureratios;
  MPI_shared_array<double> allcont_stimfactor_edgepart;
  MPI_shared_array<double> allcont_nnlevel;
  MPI_shared_array<std::uint64_t> allcont_keepbits;
  MPI_shared_array<double> chi_ff_nnionpart;
  MPI_shared_array<double> allphixstargets_corrphotoioncoeff;
};
inline CellCacheBacking cellcache_backing{};

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
inline double optical_depth_is_thick{NAN};
inline int num_grey_timesteps{-1};
inline int n_titer{1};
inline bool lte_iteration{false};

}  // namespace globals

// DO_TITER mode: average an estimator with its saved value from the previous timestep iteration
// (if one exists) and store the result as the new saved value.
inline void titer_average(double& value, double& saved) {
  if (saved >= 0) {
    value = (value + saved) / 2;
  }
  saved = value;
}

[[nodiscard]] inline auto get_max_threads() -> int {
#ifdef _OPENMP
  return omp_get_max_threads();
#elif defined(STDPAR_ON) && !defined(GPU_ON)
  // for GPU_ON mode, this call would return the number of CPU threads, which are not used.
  return static_cast<int>(std::thread::hardware_concurrency());
#else
  return 1;
#endif
}

[[nodiscard]] inline auto get_thread_num() -> int {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

// Word size of the CellCache::allcont_keepbits packing. Changing it means auditing every user of the
// layout: the two helpers below, and calculate_chi_bf_gammacontr()'s walk in rpkt.cc, which derives the
// word index and both window-edge masks itself so that it can keep the rest of a word in a register.
inline constexpr auto keepbits_wordbits = std::numeric_limits<std::uint64_t>::digits;

[[nodiscard]] inline auto get_allcont_keepwordcount() -> ptrdiff_t {
  return get_chunk_count(globals::nbfcontinua, keepbits_wordbits);
}

// mark a bound-free continuum as contributing to this cell
inline void set_allcont_keepbit(const std::span<std::uint64_t> keepbits, const int contindex) {
  keepbits[contindex / keepbits_wordbits] |= UINT64_C(1) << static_cast<unsigned>(contindex % keepbits_wordbits);
}

inline auto get_cellcache(const int nonemptymgi) -> globals::CellCache& {
  assert_testmodeonly(nonemptymgi >= 0);
  const int cacheslotid = cellcache_singleslot ? globals::rank_in_node : nonemptymgi;
  return globals::cellcache[cacheslotid];
}

#endif  // GLOBALS_H
