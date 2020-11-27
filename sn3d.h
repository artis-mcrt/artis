#ifndef TYPES
  #include "types.h"
#endif
#include <unistd.h>

//#define MPI_ON
#ifdef MPI_ON
  #include "mpi.h"
#endif

//#define TB_CUTS
#define DEBUG_ON
#define SILENT 1
//#define DO_TITER
//#define FORCE_LTE
#define RECORD_LINESTAT

// Polarisation for r-packets
#define DIPOLE
#define POL_ON

// Virtual packets
//#define VPKT_ON

// Generalised energy input
#define USE_ENERGYINPUTFILE

//#define HUGEE 2000018
//#define HUGEE 50000
//double buffer[HUGEE];
//float *buffer;



///fundamental constants-1
#define CLIGHT 2.99792458e10 /*Speed of light. */
#define CLIGHT2 8.9875518e20 /*Speed of light squared. */
#define H 6.6260755e-27 /* Planck constant */
#define MSUN 1.989e33 /* Solar mass */
#define LSUN 3.826e33 /* Solar luminosity */
#define MH 1.67352e-24 /* Mass of hydrogen atom. */
#define ME 9.1093897e-28 /* Mass of free electron. */
#define QE 4.80325E-10 /* //MK: Elementary charge in cgs units*/
#define MNI56 (56*MH) /*Mass of ni56*/
#define MCO56 (56*MH) /*Mass of co56*/
#define MFE52 (52*MH) /*Mass of fe52*/
#define MCR48 (48*MH) /*Mass of cr48*/
#define PI   3.141592654 /* PI - obviously. */
#define MEV  1.6021772e-6 /* MeV to ergs */
#define DAY  86400.0 /* day to seconds */
#define SIGMA_T 6.6524e-25 /* Thomson cross-section */
#define THOMSON_LIMIT 1e-2 /* Limit below which e-scattering is Thomson.*/
#define PARSEC 3.086e18 /* 1 pc */
#define EV  1.6021772e-12 /// eV to ergs
#define KB 1.38065e-16 /// Boltzmann constant
#define STEBO 5.670400e-5 /// Stefan-Boltzmann constant (data taken from NIST http://physics.nist.gov/cgi-bin/cuu/Value?eqsigma)
#define SAHACONST 2.0706659e-16 /// Saha constant

///numerical constants
#define CLIGHTSQUARED 8.9875518e20
#define TWOOVERCLIGHTSQUARED 2.2253001e-21
#define TWOHOVERCLIGHTSQUARED 1.4745007e-47
#define CLIGHTSQUAREDOVERTWOH 6.7819570e46
#define HOVERKB 4.799243681748932e-11
#define FOURPI 1.256637061600000e+01
#define ONEOVER4PI 7.957747153555701e-02
#define HCLIGHTOVERFOURPI 1.580764662876770e-17
#define OSCSTRENGTHCONVERSION 1.3473837e+21

#define isothermal_homogeneous_grid 0
//#define MPTS_MODEL 10000


#define MXGRID 100     /* Max number of grid cells in x-direction.*/
#define MYGRID 100     /* Max number of grid cells in y-direction.*/
#define MZGRID 100     /* Max number of grid cells in z-direction.*/
//#define MGRID 1000000  /* Max number of grid cells.*/
#define MTSTEP 200     /* Max number of time steps. */
#define MLINES 500000  /// Increase linelist by this blocksize

#define GRID_UNIFORM 1 /* Simple cuboidal cells. */
#define RHO_UNIFORM 1 /* Constant density. */
#define RHO_1D_READ 2 /* Read model. */
#define RHO_2D_READ 4 /* Read model. */
#define RHO_3D_READ 3 /* Read model. */


double ENICKEL;
double ECOBALT;
double ECOBALT_GAMMA;
double E48CR;
double E48V;
#define E52FE (0.86*MEV)
#define E52MN (3.415*MEV)

#define TNICKEL (8.80*DAY)
#define TCOBALT (113.7*DAY)
#define T48CR   (1.29602*DAY)
#define T48V    (23.0442*DAY)
#define T52FE   (0.497429*DAY)
#define T52MN   (0.0211395*DAY)


#define TYPE_ESCAPE        32
#define TYPE_NICKEL_PELLET 100
#define TYPE_COBALT_PELLET 101
#define TYPE_48CR_PELLET   102
#define TYPE_48V_PELLET    103
#define TYPE_52FE_PELLET   104
#define TYPE_52MN_PELLET   105
#define TYPE_COBALT_POSITRON_PELLET   106
#define TYPE_GENERIC_ENERGY_PELLET 200
#define TYPE_GAMMA         10
#define TYPE_RPKT          11
#define TYPE_KPKT          12
#define TYPE_MA            13
#define TYPE_EMINUS        20
#define TYPE_PRE_KPKT      120
#define TYPE_GAMMA_KPKT      121

#define COOLING_UNDEFINED       -99
#define COOLINGTYPE_FF          880
#define COOLINGTYPE_FB          881
#define COOLINGTYPE_COLLEXC     882
#define COOLINGTYPE_COLLION     883
#define HEATINGTYPE_FF          884
#define HEATINGTYPE_BF          885
#define HEATINGTYPE_COLLDEEXC   886
#define HEATINGTYPE_COLLRECOMB  887
#define COOLINGCUT              0.99 //1.01
#define TAKE_N_COOLINGTERMS     100
//#define TAKE_N_BFCONTINUA       900 //40 //900 //20 //900 //40 //20

#define RPKT_EVENTTYPE_BB 550
#define RPKT_EVENTTYPE_CONT 551


#define NEG_X  101
#define POS_X  102
#define NEG_Y  103
#define POS_Y  104
#define NEG_Z  105
#define POS_Z  106
#define NONE   107

#define PACKET_SAME -929 //MUST be negative

#define NI_GAM_LINE_ID   1
#define CO_GAM_LINE_ID   2
#define FAKE_GAM_LINE_ID 3
#define CR48_GAM_LINE_ID 4
#define V48_GAM_LINE_ID  5

#define MAX_RSCAT 50000
#define MIN_XS 1e-40

extern gsl_rng *rng; /* pointer for random number generator */

int nxgrid, nygrid, nzgrid; /* actual grid dimensions to use */
int ngrid;
int grid_type, model_type;

int nprocs;      /// Global variable which holds the number of MPI processes
int rank_global; /// Global variable which holds the rank of the active MPI process
//int nprocs_exspec;
int npkts;
int nesc; //number of packets that escape during current timestep  ///ATMOIC

double xmax, ymax, zmax;
double mtot, vmax, rmax; /* Total mass and outer velocity/radius. */
double fe_sum, rho_sum; //MK: could now be declared locally
double mni56; /*Total mass of Ni56 in the ejecta. */
double mco56; /*Total mass of Co56 in the ejecta. */
double mfe52; /*Total mass of Fe52 in the ejecta. */
double mcr48; /*Total mass of Cr48 in the ejecta. */
double mfeg; /*Total mass of Fe group elements in ejecta. */
double tmax;      /// End time of current simulation.
double tmin;      /// Start time of current simulation.

int ntstep;       /// Number of timesteps
int itstep;       /// Inital timestep's number
int ftstep;       /// Final timestep's number

int ntbins, nnubins; //number of bins for spectrum
double nu_min_r, nu_max_r; //limits on frequency range for r-pkt spectrum

int ntlcbins; //number of bins for light curve

double dlognusyn; //logarithmic interval for syn
double nusyn_min, nusyn_max; //limits on range for syn
int nfake_gam; //# of fake gamma ray lines for syn

/// New variables for other opacity cases, still grey.
double opcase3_normal;           ///MK: normalisation factor for opacity_case 3
double rho_crit_para, rho_crit;  ///MK: free parameter for the selection of the critical opacity in opacity_case 3
double rho_crit;                 ///MK: critical opacity in opacity_case 3 (could now be declared locally)

/// Some often used helper functions are declared globally
void printout(char *fmt, ...);   ///MK: printout should be used instead of printf throughout the whole code
                                 ///for output messages. Therefore it must be available globally
int get_nions(int element);
int get_nlevels(int element, int ion);
int get_ionisinglevels(int element, int ion);
int get_bfcontinua(int element, int ion);



double get_volinit_modelcell(int modelgridindex);
float get_rhoinit(int modelgridindex);
float get_rho(int modelgridindex);
float get_nne(int modelgridindex);
float get_nnetot(int modelgridindex);
float get_fni(int modelgridindex);
float get_fco(int modelgridindex);
float get_f52fe(int modelgridindex);
float get_f48cr(int modelgridindex);
float get_ffe(int modelgridindex);
float get_fnistable(int modelgridindex);
float get_fcostable(int modelgridindex);
float get_ffestable(int modelgridindex);
float get_fmnstable(int modelgridindex);
float get_fcrstable(int modelgridindex);
float get_fvstable(int modelgridindex);
float get_ftirstable(int modelgridindex);
float get_kappagrey(int modelgridindex);
float get_Te(int modelgridindex);
float get_TR(int modelgridindex);
float get_TJ(int modelgridindex);
float get_W(int modelgridindex);
#ifdef USE_ENERGYINPUTFILE
  float get_modelcell_energydensity_init(int modelgridindex);
#endif
void set_rhoinit(int modelgridindex, float x);
void set_rho(int modelgridindex, float x);
void set_nne(int modelgridindex, float x);
void set_nnetot(int modelgridindex, float x);
void set_fni(int modelgridindex, float x);
void set_fco(int modelgridindex, float x);
void set_f52fe(int modelgridindex, float x);
void set_f48cr(int modelgridindex, float x);
void set_ffe(int modelgridindex, float x);
void set_fnistable(int modelgridindex, float x);
void set_fcostable(int modelgridindex, float x);
void set_ffestable(int modelgridindex, float x);
void set_fmnstable(int modelgridindex, float x);
void set_fcrstable(int modelgridindex, float x);
void set_fvstable(int modelgridindex, float x);
void set_ftistable(int modelgridindex, float x);
void set_kappagrey(int modelgridindex, float x);
void set_Te(int modelgridindex, float x);
void set_TR(int modelgridindex, float x);
void set_TJ(int modelgridindex, float x);
void set_W(int modelgridindex, float x);
#ifdef USE_ENERGYINPUTFILE
  void set_modelcell_energydensity_init(int modelgridindex, float x);
#endif






/// New variables for the non-grey case
//FILE *ldist_file;
int debug_packet;                /// activate debug output for this packet if non negative
//double global_threshold;         /// global variable to transfer a threshold frequency to another function
//char filename[100];               /// this must be long enough to hold "packetsxx.tmp" where xx is the number of "middle" iterations
//int pktnumberoffset;
int n_middle_it;
#define MINDENSITY 1e-40         /// Minimum cell density. Below cells are treated as empty.
#define MINPOP 1e-30

PKT pkt[MPKTS];
//PKT *exchangepkt_ptr;            ///global variable to transfer a local pkt to another function



mastate_t *mastate;

#define MA_ACTION_NONE 0
#define MA_ACTION_RADDEEXC 1
#define MA_ACTION_COLDEEXC 2
#define MA_ACTION_RADRECOMB 3
#define MA_ACTION_COLRECOMB 4
#define MA_ACTION_INTERNALDOWNSAME 5
#define MA_ACTION_INTERNALDOWNLOWER 6
#define MA_ACTION_INTERNALUPSAME 7
#define MA_ACTION_INTERNALUPHIGHER 8



double wid_init;





CELL cell[MGRID+1];
//int *nonemptycells;  /// Array which contains all the non-empty cells cellnumbers
//int nnonemptycells;  /// Total number of non-empty cells

//samplegrid_t samplegrid[MSAMPLEGRID];
//extern float rhosum[MSAMPLEGRID];
//extern float T_Rsum[MSAMPLEGRID];
//extern float T_esum[MSAMPLEGRID];
//extern float Wsum[MSAMPLEGRID];
//extern float T_Dsum[MSAMPLEGRID];
//extern int associatedcells[MSAMPLEGRID];
//float rhosum2[MSAMPLEGRID];
//float T_Rsum2[MSAMPLEGRID];
//float T_esum2[MSAMPLEGRID];
//float Wsum2[MSAMPLEGRID];
//float T_Dsum2[MSAMPLEGRID];
//int associatedcells2[MSAMPLEGRID];



struct time
{
  double start; /*time at start of this timestep.*/
  double width; /*Width of timestep.*/
  double mid; /*Mid time in step - computed logarithmically.*/
  int pellet_decays; /*Number of pellets that decay in this time step. */ ///ATOMIC
  double gamma_dep; /* cmf gamma ray energy deposition rate */            ///ATOMIC
  double positron_dep; /* cmf positron energy deposition rate */            ///ATOMIC
  double dep; /* instantaenous energy deposition rate in all decays*/     ///ATOMIC
  double cmf_lum; /* cmf luminosity light curve */                        ///ATOMIC
} time_step[MTSTEP];


struct gamma_spec
{
  double energy[MGAM_LINES];
  double probability[MGAM_LINES];
  int nlines;
} cobalt_spec, nickel_spec, fakeg_spec, cr48_spec, v48_spec;

LIST gam_line_list;

#define RED_OF_LIST -956  //must be negative

double syn_dir[3]; /* vector pointing from origin to observer */

#define NRAYS_SYN 1 /* number of rays traced in a syn claculation */


RAY rays[NRAYS_SYN];



#define ACTIVE 1
#define WAITING 2
#define FINISHED 3

#define MSYN_TIME 100
int nsyn_time;
double time_syn[MSYN_TIME];

double DeltaA; /* area used for syn */ ///possible to do this as local variable?

#define EMISS_MAX 2 /* Maxmimum number of frequency points in
			     grid used to store emissivity estimators. */
int emiss_offset; /*the index in the line list of the 1st line for which
		    an emissivity estimator is recorded*/
int emiss_max; /* actual number of frequency points in emissivity grid */


modelgrid_t modelgrid[MMODELGRID+1];
/*
float modelgrid_rhoinit[MMODELGRID];
float modelgrid_rho[MMODELGRID];
//modelgrid nn_tot
float modelgrid_nne[MMODELGRID];
float modelgrid_nnetot[MMODELGRID];
float modelgrid_fni[MMODELGRID];
float modelgrid_fco[MMODELGRID];
float modelgrid_ffe[MMODELGRID];
float modelgrid_fnistable[MMODELGRID];
float modelgrid_fcostable[MMODELGRID];
float modelgrid_kappagrey[MMODELGRID];
float modelgrid_Te[MMODELGRID];
float modelgrid_TR[MMODELGRID];
float modelgrid_TJ[MMODELGRID];
float modelgrid_W[MMODELGRID];
float modelgrid_composition[MMODELGRID];

float modelgrid_abundance[MMODELGRID*MELEMENTS];         /// Abundance of the element (by mass!).
float modelgrid_groundlevelpop[MMODELGRID*MELEMENTS*MIONS];   /// Pointer to an array of floats which contains the groundlevel populations
                          /// of all included ionisation stages for the element.
float modelgrid_partfunct[MMODELGRID*MELEMENTS*MIONS];        /// Pointer to an array of floats which contains the partition functions
*/

/// THESE ARE THE GRID BASED ESTIMATORS
float compton_emiss[MMODELGRID+1][EMISS_MAX];  /// Volume estimator for the compton emissivity                     ///ATOMIC
double rpkt_emiss[MMODELGRID+1];                /// Volume estimator for the rpkt emissivity                        ///ATOMIC
double J[MMODELGRID+1];
double energy_deposition[MMODELGRID+1];
#ifdef DO_TITER
  double J_reduced_save[MMODELGRID+1];
#endif


#ifdef FORCE_LTE
  double redhelper[MMODELGRID+1];
#endif
#ifndef FORCE_LTE
  double nuJ[MMODELGRID+1];
  double ffheatingestimator[MMODELGRID+1];
  double colheatingestimator[MMODELGRID+1];
  double gammaestimator[(MMODELGRID+1)*MELEMENTS*MIONS];
  double bfheatingestimator[(MMODELGRID+1)*MELEMENTS*MIONS];
  double corrphotoionrenorm[(MMODELGRID+1)*MELEMENTS*MIONS];
  double redhelper[(MMODELGRID+1)*MELEMENTS*MIONS];


  #ifdef DO_TITER
    double nuJ_reduced_save[MMODELGRID];
    double ffheatingestimator_save[MMODELGRID];
    double colheatingestimator_save[MMODELGRID];
    double gammaestimator_save[MMODELGRID*MELEMENTS*MIONS];
    double bfheatingestimator_save[MMODELGRID*MELEMENTS*MIONS];
  #endif

  //double mabfcount[MGRID],mabfcount_thermal[MGRID], kbfcount[MGRID],kbfcount_ion[MGRID],kffcount[MGRID], kffabs[MGRID],kbfabs[MGRID],kgammadep[MGRID];
  //double matotem[MGRID],maabs[MGRID];
#endif

#ifdef RECORD_LINESTAT
  int *ecounter,*acounter,*linestat_reduced;
#endif


int do_comp_est; /* 1 = compute compton emissivity estimators. 0 = don't */
int do_r_lc;     /* If not set to 1 then the opacity for r-packets is 0. */
int do_rlc_est;  /* 1= compute estimators for the r-pkt light curve.
                    2 = compute estimators with opacity weights
                    3 = compute estimators, but use only for gamma-heating rate */


int n_out_it; /* # of sets of 1,000,000 photons to run. */

int file_set; /* 1 if the output files already exist. 0 otherwise. */

int npts_model; /*number of points in 1-D input model */
double vout_model[MMODELGRID], rho_model[MMODELGRID], fni_model[MMODELGRID], fco_model[MMODELGRID], ffegrp_model[MMODELGRID], f48cr_model[MMODELGRID],f52fe_model[MMODELGRID];
double abund_model[MMODELGRID][30];
double t_model; /* time at which densities in input model are correct. */
int ncoord1_model, ncoord2_model; /*For 2d model, the input grid dimensions */
double dcoord1, dcoord2; /*spacings of a 2d model grid - must be uniform grid */

#ifdef USE_ENERGYINPUTFILE
  double modelcell_energydensity[MMODELGRID]; // energy in model grid cell read from energydistribution.txt / initial model cell volume
  double etot_fromenergyfile; // total model energy -- used to initialise pellets. Read from energydistribution.txt

  int ntimes_energydep; // number of times included in energyrate.txt
  float time_energydep[MTSTEP]; // times in seconds from energyrate.txt
  float energy_fraction_deposited[MTSTEP]; // fraction of energy deposited by time from energyrate.txt
#endif

//#define MPTS_MODEL_3D 8000000


double CLIGHT_PROP; /* Speed of light for ray travel. Physically = CLIGHT but
		can be changed for testing. */

double gamma_grey; /* set to -ve for proper treatment. If possitive, then
		      gamma_rays are treated as grey with this opacity. */

double min_den;

#define GREY_OP 0.1

//float cont[MGRID+1];
double max_path_step;

int opacity_case; /* 0 normally, 1 for Fe-grp dependence. */
                  ///MK: 2 for Fe-grp dependence and proportional to 1/rho
                  ///MK: 3 combination of 1 & 2 depending on a rho_crit
                  ///MK: 4 non-grey treatment


double dlogt;




/// ATOMIC DATA
///============================================================================
int nelements,nlines,includedions;


/// Global pointer to beginning of atomic data. This is used as the starting point to fill up
/// the atomic data in input.c after it was read in from the database.
elementlist_entry *elements;
/// Global pointer to beginning of linelist
linelist_entry *linelist;
/// Global pointer to beginning of the bound-free list
bflist_t *bflist;







rpkt_cont_opacity_struct *kappa_rpkt_cont;






/// Coolinglist
///============================================================================



//coolinglist_entry *globalcoolinglist;
//coolinglist_entry *globalheatinglist;
//double totalcooling;
int ncoolingterms;
int importantcoolingterms;                /// Number of important cooling terms


coolingrates_t *coolingrates;


heatingrates_t *heatingrates;
//double heating_col;

/// PHIXSLIST
///============================================================================

phixslist_t *phixslist;
//extern groundphixslist_t *groundphixslist;
int nbfcontinua,nbfcontinua_ground; ///number of bf-continua
//int importantbfcontinua;
int nphixspoints;

/// Cell history
///============================================================================

cellhistory_struct *cellhistory;          /// Global pointer to the beginning of the cellhistory stack.
//extern int histindex;                            /// Global index variable to acces the cellhistory stack.
//#define CELLHISTORYSIZE 2                /// Size of the cellhistory stack.


/// Debug/Analysis Counter
int ma_stat_activation_collexc;
int ma_stat_activation_collion;
int ma_stat_activation_bb;
int ma_stat_activation_bf;
int ma_stat_deactivation_colldeexc;
int ma_stat_deactivation_collrecomb;
int ma_stat_deactivation_bb;
int ma_stat_deactivation_fb;
int k_stat_to_ma_collexc;
int k_stat_to_ma_collion;
int k_stat_to_r_ff;
int k_stat_to_r_fb;
int k_stat_to_r_bb;
int k_stat_from_ff;
int k_stat_from_bf;
int k_stat_from_gamma;
int k_stat_from_eminus;
int k_stat_from_earlierdecay;
int escounter;
int resonancescatterings;
int cellcrossings;
int upscatter;
int downscatter;
int updatecellcounter;
int coolingratecalccounter;
//extern int propagationcounter;

//extern int debuglevel;
int debuglevel;


//int currentcell; ///was used for an fdf-solver


/// Rate coefficients
///============================================================================
#define TABLESIZE 100 //200 //100
#define MINTEMP 3500.
#define MAXTEMP 140000. //1000000.
double T_step;
double T_step_log;

int homogeneous_abundances;

#define NPHIXSPOINTS 100

///Constants for van-Regmorter approximation. Defined in input.c
double C_0;
double H_ionpot;

int continue_simulation;
extern int tid;
int nthreads;
double nu_rfcut;
int n_lte_timesteps;
double cell_is_optically_thick;
int n_grey_timesteps;
int n_titer;
short initial_iteration;
int max_bf_continua;
int n_kpktdiffusion_timesteps;
float kpktdiffusion_timescale;

extern FILE *output_file;
//extern short output_file_open;

transitions_t *transitions;
int maxion;
FILE *tau_file;
FILE *tb_file;
FILE *heating_file;
FILE *estimators_file;

//double *J_below_table,*J_above_table,*nuJ_below_table,*nuJ_above_table;
extern short neutral_flag;

short elements_uppermost_ion[MTHREADS][MELEMENTS];  /// Highest ionisation stage which has a decent population for a particular element
                                                   /// in a given cell. Be aware that this must not be used outside of the update_grid
                                                   /// routine and their daughters.

extern short use_cellhist;

#ifdef _OPENMP
  #pragma omp threadprivate(tid,use_cellhist,neutral_flag,rng,output_file)
#endif
