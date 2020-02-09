#include "types.h"
#include "sn3d.h"
#include "globals.h"

double syn_dir[3]; // vector pointing from origin to observer

//#define NRAYS_SYN 1 // number of rays traced in a syn calculation

//RAY rays[NRAYS_SYN];
struct time time_step[MTSTEP];

int nsyn_time;
double time_syn[MSYN_TIME];

int emiss_offset;   // the index in the line list of the 1st line for which
                    // an emissivity estimator is recorded
int emiss_max;      // actual number of frequency points in emissivity grid


__managed__ modelgrid_t modelgrid[MMODELGRID + 1];

/// THESE ARE THE GRID BASED ESTIMATORS
float compton_emiss[MMODELGRID+1][EMISS_MAX];  /// Volume estimator for the compton emissivity
double rpkt_emiss[MMODELGRID+1];                /// Volume estimator for the rpkt emissivity


#if (!NO_LUT_PHOTOION)
  double corrphotoionrenorm[MMODELGRID * MELEMENTS * MIONS];
  double gammaestimator[MMODELGRID * MELEMENTS * MIONS];
#endif
#if (!NO_LUT_BFHEATING)
  double bfheatingestimator[MMODELGRID * MELEMENTS * MIONS];
#endif
#ifdef FORCE_LTE
  // don't use the variables below in LTE mode, just declare them here so the code compiles
  double *ffheatingestimator;
#else
  double ffheatingestimator[MMODELGRID + 1];
  double colheatingestimator[MMODELGRID + 1];

  #ifdef DO_TITER
    double ffheatingestimator_save[MMODELGRID];
    double colheatingestimator_save[MMODELGRID];
    double gammaestimator_save[MMODELGRID * MELEMENTS * MIONS];
    double bfheatingestimator_save[MMODELGRID * MELEMENTS * MIONS];
  #endif
#endif

#ifdef RECORD_LINESTAT
  int *ecounter;
  int *acounter;
  int *linestat_reduced;
#endif


bool file_set; // 1 if the output files already exist. 0 otherwise.

bool do_comp_est; // 1 = compute compton emissivity estimators. 0 = don't
bool do_r_lc;     // If not set to 1 then the opacity for r-packets is 0.
int do_rlc_est;  // 1 = compute estimators for the r-pkt light curve.
                 // 2 = compute estimators with opacity weights
                 // 3 = compute estimators, but use only for gamma-heating rate


int n_out_it; // # of sets of 1,000,000 photons to run.

int npts_model; // number of points in 1-D input model
double vout_model[MMODELGRID];
double t_model; // time at which densities in input model are correct.
int ncoord1_model, ncoord2_model; // For 2D model, the input grid dimensions
double dcoord1, dcoord2; // spacings of a 2D model grid - must be uniform grid

double CLIGHT_PROP; // Speed of light for ray travel. Physically = CLIGHT but
                    // can be changed for testing.

double gamma_grey; // set to -ve for proper treatment. If possitive, then
                   // gamma_rays are treated as grey with this opacity.

double min_den;

double max_path_step;

int opacity_case; // 0 normally, 1 for Fe-grp dependence.
                  ///MK: 2 for Fe-grp dependence and proportional to 1/rho
                  ///MK: 3 combination of 1 & 2 depending on a rho_crit
                  ///MK: 4 non-grey treatment


double dlogt;


/// ATOMIC DATA
///============================================================================
int maxion;

short elements_uppermost_ion[MTHREADS][MELEMENTS]; /// Highest ionisation stage which has a decent population for a particular element
                                                   /// in a given cell. Be aware that this must not be used outside of the update_grid
                                                   /// routine and their daughters.
__managed__ int nelements;
__managed__ int nlines;
__managed__ int includedions;
__managed__ elementlist_entry *elements;
__managed__ linelist_entry *linelist;
__managed__ bflist_t *bflist;

__managed__ rpkt_cont_opacity_struct *kappa_rpkt_cont;

/// Coolinglist
int ncoolingterms;
int importantcoolingterms;

/// PHIXSLIST
///============================================================================

__managed__ phixslist_t *phixslist;
__managed__ int nbfcontinua;
__managed__ int nbfcontinua_ground; ///number of bf-continua
__managed__ int NPHIXSPOINTS;
__managed__ double NPHIXSNUINCREMENT;

cellhistory_struct *cellhistory;

/// Debug/Analysis Counters
int ma_stat_activation_collexc;
int ma_stat_activation_collion;
int ma_stat_activation_ntcollexc;
int ma_stat_activation_ntcollion;
int ma_stat_activation_bb;
int ma_stat_activation_bf;
int ma_stat_activation_fb;
int ma_stat_deactivation_colldeexc;
int ma_stat_deactivation_collrecomb;
int ma_stat_deactivation_bb;
int ma_stat_deactivation_fb;
int ma_stat_internaluphigher;
int ma_stat_internaluphighernt;
int ma_stat_internaldownlower;
int k_stat_to_ma_collexc;
int k_stat_to_ma_collion;
int k_stat_to_r_ff;
int k_stat_to_r_fb;
int k_stat_to_r_bb;
int k_stat_from_ff;
int k_stat_from_bf;
int nt_stat_from_gamma;
int k_stat_from_earlierdecay;
int escounter;
int resonancescatterings;
int cellcrossings;
int upscatter;
int downscatter;
int updatecellcounter;
int coolingratecalccounter;

int debuglevel;

int ncoordgrid[3]; /// propagration grid dimensions
int ngrid;
int grid_type;
char coordlabel[3];

enum model_types model_type;

int nprocs;      /// Global variable which holds the number of MPI processes
int rank_global; /// Global variable which holds the rank of the active MPI process
int npkts;
int nesc; //number of packets that escape during current timestep

double coordmax[3];
double mtot;
double vmax;
double rmax;  /// Total mass and outer velocity/radius
double totmassradionuclide[RADIONUCLIDE_COUNT]; /// total mass of each radionuclide in the ejecta
double mfeg;              /// Total mass of Fe group elements in ejecta
double tmax;              /// End time of current simulation
double tmin;              /// Start time of current simulation

int ntstep;       /// Number of timesteps
int itstep;       /// Initial timestep's number
int ftstep;       /// Final timestep's number
__managed__ int nts_global;   /// Current time step

int ntbins, nnubins; //number of bins for spectrum
double nu_min_r, nu_max_r; //limits on frequency range for r-pkt spectrum

int ntlcbins; //number of bins for light curve

double nusyn_min, nusyn_max; //limits on range for syn
int nfake_gam; //# of fake gamma ray lines for syn

/// New variables for other opacity cases, still grey.
double opcase3_normal;           ///MK: normalisation factor for opacity_case 3
double rho_crit_para;            ///MK: free parameter for the selection of the critical opacity in opacity_case 3
double rho_crit;                 ///MK: critical opacity in opacity_case 3 (could now be declared locally)


/// New variables for the non-grey case
int debug_packet;                /// activate debug output for this packet if non negative
int n_middle_it;

int total_nlte_levels;            ///total number of nlte levels
int n_super_levels;

mastate_t *mastate;

CELL cell[MGRID+1];


bool homogeneous_abundances;

bool simulation_continued_from_saved;
int nthreads;
double nu_rfcut;
int n_lte_timesteps;
double cell_is_optically_thick;
int n_grey_timesteps;
int n_titer;
__managed__ bool initial_iteration;
int max_bf_continua;
int n_kpktdiffusion_timesteps;
float kpktdiffusion_timescale;
