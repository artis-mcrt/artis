#ifndef SN3D_H
#define SN3D_H

#include <unistd.h>
#include <stdarg.h>  /// MK: needed for printout()
#include <stdbool.h>
#include <gsl/gsl_integration.h>

#define DEBUG_ON
// #define DO_TITER
// #define FORCE_LTE

/// non-thermal ionisation
static const bool NT_ON = false;

/// use the detailed Spencer-Fano solver instead of the work function approximation
static const bool NT_SOLVE_SPENCERFANO = false;

// non-LTE population solver
static const bool NLTE_POPS_ON = false;

// solve the NLTE population matrix equation simultaneously for levels of all ions of an element
static const bool NLTE_POPS_ALL_IONS_SIMULTANEOUS = false;

// maximum number of NLTE/Te/Spencer-Fano iterations
static const int NLTEITER = 30;

// if using this, avoid look up tables and switch on the direct integration options below
// (since LUTs created with Planck function J_nu)
static const bool MULTIBIN_RADFIELD_MODEL_ON = false;

// store Jb_lu estimators for particular lines chosen in radfield.c:radfield_init()
#define DETAILED_LINE_ESTIMATORS_ON false

// store detailed bound-free rate estimators
#define DETAILED_BF_ESTIMATORS_ON false

// extremely slow and memory consuming - for debugging only
// not safe for MPI or OpenMP - single process and single thread only!
// this will output a list of contributions to each bound-free rate estimator
// with each packet emission type ranked by their contribution to the rate
#define DETAILED_BF_ESTIMATORS_BYTYPE false

// dynamically calculate photoionization rates for the current radiation field
// instead of interpolating values from a lookup table for a blackbody radiation field
#define NO_LUT_PHOTOION true

// as above for bound-free heating
#define NO_LUT_BFHEATING true

#define SEPARATE_STIMRECOMB false // if false, then stimulated recombination is treated as negative photoionisation

#define DIRECT_COL_HEAT
#define NO_INITIAL_PACKETS
#define RECORD_LINESTAT

/// Rate coefficients
///============================================================================
#define TABLESIZE 100 //200 //100
#define MINTEMP 1000.
#define MAXTEMP 30000. //1000000.

// Polarisation for real packets
// #define DIPOLE
// #define POL_ON

// Polarisation for virtual packets
// #define VPKT_ON

static const size_t GSLWSIZE = 16384;

#include "types.h"

#if (DETAILED_BF_ESTIMATORS_ON && !NO_LUT_PHOTOION)
  #error Must use NO_LUT_PHOTOION with DETAILED_BF_ESTIMATORS_ON
#endif

#if !defined DO_EXSPEC && !defined MPI_ON
  // #define MPI_ON //only needed for debugging MPI, the makefile will switch this on
#endif

#ifdef MPI_ON
  #include "mpi.h"
#endif

//#define _OPENMP
#ifdef _OPENMP
  #include "omp.h"
#endif


/// fundamental constants
#define CLIGHT        2.99792458e+10    /// Speed of light [cm/s]
#define H             6.6260755e-27     /// Planck constant [erg s]
#define MSUN          1.98855e+33       /// Solar mass [g]
#define LSUN          3.826e+33         /// Solar luminosity [erg/s]
#define MH            1.67352e-24       /// Mass of hydrogen atom [g]
#define ME            9.1093897e-28     /// Mass of free electron [g]
#define QE            4.80325E-10       /// elementary charge in cgs units [statcoulomb]
#define PI            3.1415926535987
#define EV            1.6021772e-12     /// eV to ergs [eV/erg]
#define MEV           1.6021772e-6      /// MeV to ergs [MeV/erg]
#define DAY           86400.0           /// day to seconds [s/day]
#define SIGMA_T       6.6524e-25        /// Thomson cross-section
#define THOMSON_LIMIT 1e-2              /// Limit below which e-scattering is Thomson
#define PARSEC        3.0857e+18        /// pc to cm [pc/cm]
#define KB            1.38064852e-16    /// Boltzmann constant [erg/K]
#define STEBO         5.670400e-5       /// Stefan-Boltzmann constant [erg cm^−2 s^−1 K^−4.]
                                        /// (data from NIST http://physics.nist.gov/cgi-bin/cuu/Value?eqsigma)
#define SAHACONST     2.0706659e-16     /// Saha constant

/// numerical constants
#define CLIGHTSQUARED         8.9875518e+20   /// Speed of light squared [cm^2/s^2]
#define TWOOVERCLIGHTSQUARED  2.2253001e-21
#define TWOHOVERCLIGHTSQUARED 1.4745007e-47
#define CLIGHTSQUAREDOVERTWOH 6.7819570e+46

#define ONEOVERH              1.509188961e+26
#define HOVERKB               4.799243681748932e-11
#define FOURPI                1.256637061600000e+01
#define ONEOVER4PI            7.957747153555701e-02
#define HCLIGHTOVERFOURPI     1.580764662876770e-17
#define OSCSTRENGTHCONVERSION 1.3473837e+21

#define H_ionpot (13.5979996 * EV)


// warning! when these are lower case, the variables represent total masses in the ejecta
// in upper case, these constants are the masses of each nucleus
#define MNI57 (57 * MH)                            /// Mass of Ni57
#define MCO57 (57 * MH)                            /// Mass of Co57
#define MNI56 (56 * MH)                            /// Mass of Ni56
#define MCO56 (56 * MH)                            /// Mass of Co56
#define MFE52 (52 * MH)                            /// Mass of Fe52
#define MCR48 (48 * MH)                            /// Mass of Cr48

#define MTSTEP 200       // Max number of time steps.
#define MLINES 500000    // Increase linelist by this blocksize

#define GRID_UNIFORM 1 // Simple cuboidal cells.
#define GRID_SPHERICAL1D 2 // radial shells

double E56NI;
double E56CO;
double E56CO_GAMMA;
double E57NI;
double E57NI_GAMMA;
double E57CO;
double E48CR;
double E48V;
#define E52FE (0.86*MEV)
#define E52MN (3.415*MEV)

/* mean lifetimes */
#define T56NI   (8.80*DAY)
#define T56CO   (113.7*DAY)
#define T57NI   (51.36*60)
#define T57CO   (392.03*DAY)
#define T48CR   (1.29602*DAY)
#define T48V    (23.0442*DAY)
#define T52FE   (0.497429*DAY)
#define T52MN   (0.0211395*DAY)

#define COOLING_UNDEFINED       -99

#define COOLINGCUT              0.99 //1.01
//#define TAKE_N_BFCONTINUA       900 //40 //900 //20 //900 //40 //20

#define RPKT_EVENTTYPE_BB 550
#define RPKT_EVENTTYPE_CONT 551

#define PACKET_SAME -929 //MUST be negative

#define MAX_RSCAT 50000
#define MIN_XS 1e-40

int ncoordgrid[3]; /// actual grid dimensions to use
int ngrid;
int grid_type;
char coordlabel[3];

enum model_types {
  RHO_UNIFORM = 1,  // Constant density.
  RHO_1D_READ = 2,  // Read model.
  RHO_2D_READ = 4,  // Read model.
  RHO_3D_READ = 3,  // Read model.
};

enum model_types model_type;

int nprocs;      /// Global variable which holds the number of MPI processes
int rank_global; /// Global variable which holds the rank of the active MPI process
int npkts;
int nesc; //number of packets that escape during current timestep  ///ATOMIC

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
int nts_global;   /// Current time step

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
//FILE *ldist_file;
int debug_packet;                /// activate debug output for this packet if non negative
//double global_threshold;         /// global variable to transfer a threshold frequency to another function
int n_middle_it;
#define MINDENSITY 1e-40         /// Minimum cell density. Below cells are treated as empty.
#define MINPOP 1e-40

int total_nlte_levels;            ///total number of nlte levels
int n_super_levels;

mastate_t *restrict mastate;

CELL cell[MGRID+1];
//int *nonemptycells;  /// Array which contains all the non-empty cells cellnumbers
//int nnonemptycells;  /// Total number of non-empty cells


struct time
{
  double start; // time at start of this timestep.
  double width; // Width of timestep.
  double mid; // Mid time in step - computed logarithmically.
  double gamma_dep; // cmf gamma ray energy deposition rate             ///ATOMIC
  double positron_dep; // cmf positron energy deposition rate           ///ATOMIC
  double dep; // instantaneous energy deposition rate in all decays     ///ATOMIC
  double cmf_lum; // cmf luminosity light curve                         ///ATOMIC
  int pellet_decays; // Number of pellets that decay in this time step. ///ATOMIC
} time_step[MTSTEP];


double syn_dir[3]; // vector pointing from origin to observer

//#define NRAYS_SYN 1 // number of rays traced in a syn calculation

//RAY rays[NRAYS_SYN];

#define MSYN_TIME 100
int nsyn_time;
double time_syn[MSYN_TIME];

#define EMISS_MAX 2 // Maxmimum number of frequency points in
                    // grid used to store emissivity estimators.
int emiss_offset;   // the index in the line list of the 1st line for which
                    // an emissivity estimator is recorded
int emiss_max;      // actual number of frequency points in emissivity grid


modelgrid_t modelgrid[MMODELGRID + 1];

/// THESE ARE THE GRID BASED ESTIMATORS
float compton_emiss[MMODELGRID+1][EMISS_MAX];  /// Volume estimator for the compton emissivity                     ///ATOMIC
double rpkt_emiss[MMODELGRID+1];                /// Volume estimator for the rpkt emissivity                        ///ATOMIC


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

  //double mabfcount[MGRID],mabfcount_thermal[MGRID], kbfcount[MGRID],kbfcount_ion[MGRID],kffcount[MGRID], kffabs[MGRID],kbfabs[MGRID],kgammadep[MGRID];
  //double matotem[MGRID],maabs[MGRID];
#endif

#ifdef RECORD_LINESTAT
  int *restrict ecounter;
  int *restrict acounter;
  int *restrict linestat_reduced;
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

//#define MPTS_MODEL_3D 8000000

double CLIGHT_PROP; // Speed of light for ray travel. Physically = CLIGHT but
                    // can be changed for testing.

double gamma_grey; // set to -ve for proper treatment. If possitive, then
                   // gamma_rays are treated as grey with this opacity.

double min_den;

#define GREY_OP 0.1

double max_path_step;

int opacity_case; // 0 normally, 1 for Fe-grp dependence.
                  ///MK: 2 for Fe-grp dependence and proportional to 1/rho
                  ///MK: 3 combination of 1 & 2 depending on a rho_crit
                  ///MK: 4 non-grey treatment


double dlogt;


/// ATOMIC DATA
///============================================================================
int nelements,nlines,includedions;
elementlist_entry *restrict elements;
linelist_entry *restrict linelist;
bflist_t *restrict bflist;



rpkt_cont_opacity_struct *restrict kappa_rpkt_cont;



/// Coolinglist
///============================================================================
//coolinglist_entry *globalcoolinglist;
//coolinglist_entry *globalheatinglist;
//double totalcooling;
int ncoolingterms;
int importantcoolingterms;

//double heating_col;


/// PHIXSLIST
///============================================================================

phixslist_t *restrict phixslist;
int nbfcontinua;
int nbfcontinua_ground; ///number of bf-continua
int NPHIXSPOINTS;
double NPHIXSNUINCREMENT;

/// Cell history
///============================================================================

cellhistory_struct *restrict cellhistory;          /// Global pointer to the beginning of the cellhistory stack.

/// Debug/Analysis Counter
int ma_stat_activation_collexc;
int ma_stat_activation_collion;
int ma_stat_activation_ntcollexc;
int ma_stat_activation_ntcollion;
int ma_stat_activation_bb;
int ma_stat_activation_bf;
double stat_bf_photon_absorptions;
double stat_bf_photon_emissions;
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
int nt_stat_to_ionization;
int nt_stat_to_excitation;
int nt_stat_to_kpkt;
int escounter;
int resonancescatterings;
int cellcrossings;
int upscatter;
int downscatter;
int updatecellcounter;
int coolingratecalccounter;

int debuglevel;


bool homogeneous_abundances;

bool simulation_continued_from_saved;
int nthreads;
double nu_rfcut;
int n_lte_timesteps;
double cell_is_optically_thick;
int n_grey_timesteps;
int n_titer;
bool initial_iteration;
int max_bf_continua;
int n_kpktdiffusion_timesteps;
float kpktdiffusion_timescale;

int maxion;
// FILE *restrict tau_file;
// FILE *restrict tb_file;
// FILE *restrict heating_file;

short elements_uppermost_ion[MTHREADS][MELEMENTS]; /// Highest ionisation stage which has a decent population for a particular element
                                                   /// in a given cell. Be aware that this must not be used outside of the update_grid
                                                   /// routine and their daughters.

extern int tid;
extern bool use_cellhist;
extern bool neutral_flag;
extern gsl_rng *rng;  // pointer for random number generator
extern gsl_integration_workspace *gslworkspace;
extern FILE *restrict output_file;

#ifdef _OPENMP
  #pragma omp threadprivate(tid, use_cellhist, neutral_flag, rng, gslworkspace, output_file)
#endif

inline int printout(const char *restrict format, ...)
{
   va_list args;
   va_start(args, format);
   const int ret_status = vfprintf(output_file, format, args);
   // fprintf(output_file, "vfprintf return code %d\n", va_arg(args, char*) == NULL);
   va_end(args);

   return ret_status;
}

#ifdef DEBUG_ON
  #define assert(e) if (!(e)) { printout("%s:%u: failed assertion `%s' in function %s\n", __FILE__, __LINE__, #e, __PRETTY_FUNCTION__); abort(); }
#else
  #define	assert(e)	((void)0)
#endif


inline void gsl_error_handler_printout(const char *reason, const char *file, int line, int gsl_errno)
{
  printout("WARNING: gsl (%s:%d): %s (Error code %d)\n", file, line, reason, gsl_errno);
  // abort();
}


inline FILE *fopen_required(const char *filename, const char *mode)
{
  FILE *file = fopen(filename, mode);
  if (file == NULL)
  {
    printout("ERROR: Could not open file '%s' for mode '%s'.\n", filename, mode);
    abort();
  }
  else
    return file;
}


#endif // SN3D_H
