#ifndef GLOBALS_H
#define GLOBALS_H

#include "cuda.h"

#if CUDA_ENABLED
#include <curand_kernel.h>
#endif

#ifdef MPI_ON
#include "mpi.h"
#endif

#include "artisoptions.h"
#include "boundary.h"
#include "macroatom.h"

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

struct fullphixslist {
  double nu_edge;
  int element;
  int ion;
  int level;
  int phixstargetindex;
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
#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
  double *groundcont_gamma_contr;
#endif
  double *kappa_bf_sum;
#if (DETAILED_BF_ESTIMATORS_ON)
  double *gamma_contr;
#endif
};

struct phixstarget_entry {
  double *spontrecombcoeff;
#if (!NO_LUT_PHOTOION)
  double *corrphotoioncoeff;
#endif
#if (!NO_LUT_BFHEATING)
  double *bfheating_coeff;
#endif
  double *bfcooling_coeff;

  double probability;  // fraction of phixs cross section leading to this final level
  int levelindex;      // index of upper ion level after photoionisation
};

struct levellist_entry {
  double epsilon;               /// Excitation energy of this level relative to the neutral ground level.
  int *uptrans_lineindicies;    /// Allowed upward transitions from this level
  int *downtrans_lineindicies;  /// Allowed downward transitions from this level
  int nuptrans;
  int ndowntrans;
  double phixs_threshold;                  /// Energy of first point in the photion_xs table
  struct phixstarget_entry *phixstargets;  /// pointer to table of target states and probabilities
  float *photoion_xs;  /// Pointer to a lookup-table providing photoionisation cross-sections for this level.
  int nphixstargets;   /// length of phixstargets array:
  float stat_weight;   /// Statistical weight of this level.

  int cont_index;  /// Index of the continuum associated to this level. Negative number.
#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
  int closestgroundlevelcont;
#endif

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
  float osc_strength;
  float coll_str;
  int elementindex;     /// It's a transition of element (not its atomic number,
                        /// but the (x-1)th element included in the simulation.
  int ionindex;         /// The same for the elements ion
  int upperlevelindex;  /// And the participating upper
  int lowerlevelindex;  /// and lower levels
  bool forbidden;
};

struct bflist_t {
  int elementindex;
  int ionindex;
  int levelindex;
  int phixstargetindex;
};

struct nne_solution_paras {
  int cellnumber;
};

struct gslintegration_paras {
  double nu_edge;
  double T;
  float *photoion_xs;
};

struct rpkt_cont_opacity {
  double nu;  // frequency at which opacity was calculated
  double total;
  double es;
  double ff;
  double bf;
  double ffheating;
  // double bfheating;
  int modelgridindex;
  bool recalculate_required;  // e.g. when cell or timestep has changed
};

struct chphixstargets {
  double corrphotoioncoeff;
#if (SEPARATE_STIMRECOMB)
  double stimrecombcoeff;
#endif
};

struct chlevels {
  double processrates[MA_ACTION_COUNT];
  struct chphixstargets *chphixstargets;
  double bfheatingcoeff;
  double population;
  double *individ_rad_deexc;
  double *individ_internal_down_same;
  double *individ_internal_up_same;
};

struct chions {
  struct chlevels *chlevels;  /// Pointer to the ions levellist.
};

struct chelements {
  struct chions *chions;  /// Pointer to the elements ionlist.
};

struct cellhistory {
  double *cooling_contrib;        /// Cooling contributions by the different processes.
  struct chelements *chelements;  /// Pointer to a nested list which helds compositional
                                  /// information for all the elements=0,1,...,nelements-1
  int cellnumber;                 /// Identifies the cell the data is valid for.
  int bfheating_mgi;
};

namespace globals {
#if CUDA_ENABLED
extern __managed__ curandState curandstates[MCUDATHREADS];
#endif

extern __managed__ double syn_dir[3];  // vector pointing from origin to observer

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
extern __managed__ struct elementlist_entry *elements;
extern __managed__ struct linelist_entry *linelist;
extern __managed__ struct bflist_t *bflist;

extern __managed__ struct rpkt_cont_opacity *kappa_rpkt_cont;

extern __managed__ int ncoolingterms;

extern __managed__ double *allcont_nu_edge;
extern __managed__ struct fullphixslist *allcont;
#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
extern __managed__ struct groundphixslist *groundcont;
#endif
extern __managed__ struct phixslist *phixslist;
extern __managed__ int nbfcontinua;
extern __managed__ int nbfcontinua_ground;
extern __managed__ int NPHIXSPOINTS;
extern __managed__ double NPHIXSNUINCREMENT;

extern __managed__ struct cellhistory *cellhistory;

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

}  // namespace globals

#endif  // GLOBALS_H