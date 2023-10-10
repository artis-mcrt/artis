#ifndef GLOBALS_H
#define GLOBALS_H

#include <memory>
#include <vector>

#ifdef MPI_ON
#include <mpi.h>
#endif

#include "artisoptions.h"

struct time {
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

struct bflist_t {
  int elementindex;
  int ionindex;
  int levelindex;
  int phixstargetindex;
};

struct fullphixslist {
  double nu_edge;
  int element;
  int ion;
  int level;
  int phixstargetindex;
  int upperlevel;
  float *photoion_xs;
  double probability;
  int index_in_groundphixslist;
};

struct groundphixslist {
  double nu_edge;
  int element;
  int ion;
  int level;
  int phixstargetindex;
};

struct phixslist {
  double *groundcont_gamma_contr = nullptr;  // for either USE_LUT_PHOTOION = true or !USE_LUT_BFHEATING = false
  double *chi_bf_sum = nullptr;
  double *gamma_contr = nullptr;  // needed for DETAILED_BF_ESTIMATORS_ON
};

struct phixstarget_entry {
  double probability;  // fraction of phixs cross section leading to this final level
  int levelindex;      // index of upper ion level after photoionisation
};

struct level_transition {
  int lineindex;
  int targetlevelindex;
  float einstein_A;
  float coll_str;
  float osc_strength;
  bool forbidden;
};

struct levellist_entry {
  double epsilon;                      /// Excitation energy of this level relative to the neutral ground level.
  struct level_transition *uptrans;    /// Allowed upward transitions from this level
  struct level_transition *downtrans;  /// Allowed downward transitions from this level
  int nuptrans;
  int ndowntrans;
  double phixs_threshold;                            /// Energy of first point in the photion_xs table
  struct phixstarget_entry *phixstargets = nullptr;  /// pointer to table of target states and probabilities
  float *photoion_xs = nullptr;  /// Pointer to a lookup-table providing photoionisation cross-sections for this level.
  int nphixstargets;             /// length of phixstargets array:
  float stat_weight;             /// Statistical weight of this level.

  int cont_index;  /// Index of the continuum associated to this level. Negative number.
  int closestgroundlevelcont;

  int uniquelevelindex;
  bool metastable;  ///
};

struct ionlist_entry {
  struct levellist_entry *levels;  /// Carries information for each level: 0,1,...,nlevels-1
  int ionstage;                    /// Which ionisation stage: XI=0, XII=1, XIII=2, ...
  int nlevels;                     /// Number of levels for this ionisation stage
  int nlevels_nlte;                /// number of nlte levels for this ion
  int first_nlte;                  /// index into nlte_pops array of a grid cell
  int ionisinglevels;              /// Number of levels which have a bf-continuum
  int maxrecombininglevel;         /// level index of the highest level with a non-zero recombination rate
  int nlevels_groundterm;
  int coolingoffset;
  int ncoolingterms;
  int uniqueionindex;
  float *Alpha_sp;
  double ionpot;  /// Ionisation threshold to the next ionstage
  // int nbfcontinua;
  // ionsphixslist *phixslist;
};

struct elementlist_entry {
  ionlist_entry *ions;  /// Carries information for each ion: 0,1,...,nions-1
  int nions;            /// Number of ions for the current element
  int anumber;          /// Atomic number
  //  int uppermost_ion;                       /// Highest ionisation stage which has a decent population for a given
  //  cell
  /// Be aware that this must not be used outside of the update_grid routine
  /// and their daughters. Neither it will work with OpenMP threads.
  float abundance;              ///
  float initstablemeannucmass;  /// Atomic mass number in multiple of MH
};

struct linelist_entry {
  double nu;  /// Frequency of the line transition
  float einstein_A;
  int elementindex;     /// It's a transition of element (not its atomic number,
                        /// but the (x-1)th element included in the simulation.
  int ionindex;         /// The same for the elements ion
  int upperlevelindex;  /// And the participating upper
  int lowerlevelindex;  /// and lower levels
};

struct gslintegration_paras {
  const double nu_edge;
  const float T;
  const float *const photoion_xs;
};

struct rpkt_continuum_absorptioncoeffs {
  double nu = NAN;  // frequency at which opacity was calculated
  double total = 0.;
  double es = 0.;
  double ff = 0.;
  double bf = 0.;
  double ffheating = 0.;
  // double bfheating;
  int modelgridindex = -1;
  bool recalculate_required = true;  // e.g. when cell or timestep has changed
};

template <bool separatestimrecomb>
struct _chphixstargets {
  double corrphotoioncoeff;
};

template <>
struct _chphixstargets<true> {
  double corrphotoioncoeff;
  double separatestimrecomb;
};

using chphixstargets_t = struct _chphixstargets<SEPARATE_STIMRECOMB>;

#include "macroatom.h"

struct chlevels {
  double processrates[MA_ACTION_COUNT];
  chphixstargets_t *chphixstargets;
  double bfheatingcoeff;
  double population;
  double *sum_epstrans_rad_deexc;
  double *sum_internal_down_same;
  double *sum_internal_up_same;
};

struct chions {
  struct chlevels *chlevels;  /// Pointer to the ions levellist.
};

struct chelements {
  struct chions *chions;  /// Pointer to the elements ionlist.
};

struct cellhistory {
  double *cooling_contrib;  /// Cooling contributions by the different processes.
  struct chelements *chelements;
  struct chlevels *ch_all_levels;
  double *ch_allcont_departureratios;
  int cellnumber;  /// Identifies the cell the data is valid for.
  int bfheating_mgi;
};

namespace globals {

extern double syn_dir[3];  // vector pointing from origin to observer

extern struct time *timesteps;

extern double *rpkt_emiss;

// for USE_LUT_PHOTOION = true
extern double *corrphotoionrenorm;
extern double *gammaestimator;

// for USE_LUT_BFHEATING = true
extern double *bfheatingestimator;

extern double *ffheatingestimator;
extern double *colheatingestimator;
#ifdef DO_TITER
extern double *gammaestimator_save;
extern double *bfheatingestimator_save;
extern double *ffheatingestimator_save;
extern double *colheatingestimator_save;
#endif

extern int *ecounter;
extern int *acounter;

extern int nprocs_exspec;
extern bool do_emission_res;

extern std::unique_ptr<bool[]> startofline;

extern double gamma_kappagrey;

constexpr double GREY_OP = 0.1;

extern double max_path_step;

extern int opacity_case;

extern int nlines;
extern std::vector<struct elementlist_entry> elements;

extern const struct linelist_entry *linelist;
extern struct bflist_t *bflist;

// for USE_LUT_BFHEATING = true
extern double *bfheating_coeff;

extern struct rpkt_continuum_absorptioncoeffs *chi_rpkt_cont;

extern int ncoolingterms;

extern double *allcont_nu_edge;
extern const struct fullphixslist *allcont;

// for either USE_LUT_PHOTOION = true or !USE_LUT_BFHEATING = false
extern struct groundphixslist *groundcont;

extern struct phixslist *phixslist;
extern int nbfcontinua;
extern int nbfcontinua_ground;
extern int NPHIXSPOINTS;
extern double NPHIXSNUINCREMENT;

extern struct cellhistory *cellhistory;

#ifdef MPI_ON
extern MPI_Comm mpi_comm_node;
extern MPI_Comm mpi_comm_internode;
#endif

extern int nprocs;
extern int rank_global;

extern int node_nprocs;
extern int rank_in_node;

extern int node_count;
extern int node_id;

extern const int npkts;
extern int nesc;

extern double vmax;
extern double rmax;
extern double tmax;
extern double tmin;

extern int ntimesteps;
extern int timestep_initial;
extern int timestep_finish;
extern int timestep;

extern double opcase3_normal;
extern double rho_crit_para;
extern double rho_crit;

extern int total_nlte_levels;

extern bool simulation_continued_from_saved;
extern double nu_rfcut;
extern int num_lte_timesteps;
extern double cell_is_optically_thick;
extern int num_grey_timesteps;
extern int n_titer;
extern bool lte_iteration;
extern int n_kpktdiffusion_timesteps;
extern float kpktdiffusion_timescale;

}  // namespace globals

#endif  // GLOBALS_H