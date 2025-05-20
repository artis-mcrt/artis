#include "sn3d.h"

double ENICKEL;
double ECOBALT;
double ECOBALT_GAMMA;
double E48CR;
double E48V;

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
double rho_crit_para;  ///MK: free parameter for the selection of the critical opacity in opacity_case 3
double rho_crit;                 ///MK: critical opacity in opacity_case 3 (could now be declared locally)



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




struct time time_step[MTSTEP];


struct gamma_spec cobalt_spec, nickel_spec, fakeg_spec, cr48_spec, v48_spec;

LIST gam_line_list;

double syn_dir[3]; /* vector pointing from origin to observer */

RAY rays[NRAYS_SYN];


int nsyn_time;
double time_syn[MSYN_TIME];

double DeltaA; /* area used for syn */ ///possible to do this as local variable?

int emiss_offset; /*the index in the line list of the 1st line for which
		    an emissivity estimator is recorded*/
int emiss_max; /* actual number of frequency points in emissivity grid */


modelgrid_t modelgrid[MMODELGRID+1];

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

//#define MPTS_MODEL_3D 8000000


double CLIGHT_PROP; /* Speed of light for ray travel. Physically = CLIGHT but
		can be changed for testing. */

double gamma_grey; /* set to -ve for proper treatment. If possitive, then
		      gamma_rays are treated as grey with this opacity. */

double min_den;


float cont[MGRID+1];
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
double T_step;
double T_step_log;

int homogeneous_abundances;

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
