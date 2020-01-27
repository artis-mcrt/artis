#include "types.h"
#ifndef GLOBALS_H
#define GLOBALS_H

#include "types.h"
#include "sn3d.h"
#include "globals.h"

extern double syn_dir[3]; // vector pointing from origin to observer

//#define NRAYS_SYN 1 // number of rays traced in a syn calculation

//RAY rays[NRAYS_SYN];

extern struct time time_step[MTSTEP];

#define MSYN_TIME 100
extern int nsyn_time;
extern double time_syn[MSYN_TIME];

#define EMISS_MAX 2 // Maxmimum number of frequency points in
                    // grid used to store emissivity estimators.
extern int emiss_offset;   // the index in the line list of the 1st line for which
                    // an emissivity estimator is recorded
extern int emiss_max;      // actual number of frequency points in emissivity grid


extern modelgrid_t modelgrid[MMODELGRID + 1];

/// THESE ARE THE GRID BASED ESTIMATORS
extern float compton_emiss[MMODELGRID+1][EMISS_MAX];  /// Volume estimator for the compton emissivity                     ///ATOMIC
extern double rpkt_emiss[MMODELGRID+1];                /// Volume estimator for the rpkt emissivity                        ///ATOMIC


#if (!NO_LUT_PHOTOION)
  extern double corrphotoionrenorm[MMODELGRID * MELEMENTS * MIONS];
  extern double gammaestimator[MMODELGRID * MELEMENTS * MIONS];
#endif
#if (!NO_LUT_BFHEATING)
  extern double bfheatingestimator[MMODELGRID * MELEMENTS * MIONS];
#endif
#ifdef FORCE_LTE
  // don't use the variables below in LTE mode, just declare them here so the code compiles
  extern double *ffheatingestimator;
#else
  extern double ffheatingestimator[MMODELGRID + 1];
  extern double colheatingestimator[MMODELGRID + 1];

  #ifdef DO_TITER
    extern double ffheatingestimator_save[MMODELGRID];
    extern double colheatingestimator_save[MMODELGRID];
    extern double gammaestimator_save[MMODELGRID * MELEMENTS * MIONS];
    extern double bfheatingestimator_save[MMODELGRID * MELEMENTS * MIONS];
  #endif

  //double mabfcount[MGRID],mabfcount_thermal[MGRID], kbfcount[MGRID],kbfcount_ion[MGRID],kffcount[MGRID], kffabs[MGRID],kbfabs[MGRID],kgammadep[MGRID];
  //double matotem[MGRID],maabs[MGRID];
#endif

#ifdef RECORD_LINESTAT
  extern int *ecounter;
  extern int *acounter;
  extern int *linestat_reduced;
#endif


extern bool file_set; // 1 if the output files already exist. 0 otherwise.

extern bool do_comp_est; // 1 = compute compton emissivity estimators. 0 = don't
extern bool do_r_lc;     // If not set to 1 then the opacity for r-packets is 0.
extern int do_rlc_est;  // 1 = compute estimators for the r-pkt light curve.
                 // 2 = compute estimators with opacity weights
                 // 3 = compute estimators, but use only for gamma-heating rate


extern int n_out_it; // # of sets of 1,000,000 photons to run.

extern int npts_model; // number of points in 1-D input model
extern double vout_model[MMODELGRID];
extern double t_model; // time at which densities in input model are correct.
extern int ncoord1_model, ncoord2_model; // For 2D model, the input grid dimensions
extern double dcoord1, dcoord2; // spacings of a 2D model grid - must be uniform grid

//#define MPTS_MODEL_3D 8000000

extern double CLIGHT_PROP; // Speed of light for ray travel. Physically = CLIGHT but
                    // can be changed for testing.

extern double gamma_grey; // set to -ve for proper treatment. If possitive, then
                   // gamma_rays are treated as grey with this opacity.

extern double min_den;

#define GREY_OP 0.1

extern double max_path_step;

extern int opacity_case; // 0 normally, 1 for Fe-grp dependence.
                  ///MK: 2 for Fe-grp dependence and proportional to 1/rho
                  ///MK: 3 combination of 1 & 2 depending on a rho_crit
                  ///MK: 4 non-grey treatment


extern double dlogt;


/// ATOMIC DATA
///============================================================================
extern int maxion;
extern short elements_uppermost_ion[MTHREADS][MELEMENTS]; /// Highest ionisation stage which has a decent population for a particular element
                                                   /// in a given cell. Be aware that this must not be used outside of the update_grid
                                                   /// routine and their daughters.
extern int nelements;
extern int nlines;
extern int includedions;
extern elementlist_entry *elements;
extern linelist_entry *linelist;
extern bflist_t *bflist;

extern rpkt_cont_opacity_struct *kappa_rpkt_cont;

/// Coolinglist
///============================================================================
//coolinglist_entry *globalcoolinglist;
//coolinglist_entry *globalheatinglist;
//double totalcooling;
extern int ncoolingterms;
extern int importantcoolingterms;

//double heating_col;


/// PHIXSLIST
///============================================================================

extern phixslist_t *phixslist;
extern int nbfcontinua;
extern int nbfcontinua_ground; ///number of bf-continua
extern int NPHIXSPOINTS;
extern double NPHIXSNUINCREMENT;

/// Cell history
///============================================================================

extern cellhistory_struct *cellhistory;          /// Global pointer to the beginning of the cellhistory stack.

/// Debug/Analysis Counter
extern int ma_stat_activation_collexc;
extern int ma_stat_activation_collion;
extern int ma_stat_activation_ntcollexc;
extern int ma_stat_activation_ntcollion;
extern int ma_stat_activation_bb;
extern int ma_stat_activation_bf;
extern int ma_stat_activation_fb;
extern int ma_stat_deactivation_colldeexc;
extern int ma_stat_deactivation_collrecomb;
extern int ma_stat_deactivation_bb;
extern int ma_stat_deactivation_fb;
extern int ma_stat_internaluphigher;
extern int ma_stat_internaluphighernt;
extern int ma_stat_internaldownlower;
extern int k_stat_to_ma_collexc;
extern int k_stat_to_ma_collion;
extern int k_stat_to_r_ff;
extern int k_stat_to_r_fb;
extern int k_stat_to_r_bb;
extern int k_stat_from_ff;
extern int k_stat_from_bf;
extern int nt_stat_from_gamma;
extern int k_stat_from_earlierdecay;
extern int escounter;
extern int resonancescatterings;
extern int cellcrossings;
extern int upscatter;
extern int downscatter;
extern int updatecellcounter;
extern int coolingratecalccounter;

extern int debuglevel;

extern int ncoordgrid[3]; /// actual grid dimensions to use
extern int ngrid;
extern int grid_type;
extern char coordlabel[3];

extern enum model_types model_type;

extern int nprocs;      /// Global variable which holds the number of MPI processes
extern int rank_global; /// Global variable which holds the rank of the active MPI process
extern int npkts;
extern int nesc; //number of packets that escape during current timestep  ///ATOMIC

extern double coordmax[3];
extern double mtot;
extern double vmax;
extern double rmax;  /// Total mass and outer velocity/radius
extern double totmassradionuclide[RADIONUCLIDE_COUNT]; /// total mass of each radionuclide in the ejecta
extern double mfeg;              /// Total mass of Fe group elements in ejecta
extern double tmax;              /// End time of current simulation
extern double tmin;              /// Start time of current simulation

extern int ntstep;       /// Number of timesteps
extern int itstep;       /// Initial timestep's number
extern int ftstep;       /// Final timestep's number
extern int nts_global;   /// Current time step

extern int ntbins, nnubins; //number of bins for spectrum
extern double nu_min_r, nu_max_r; //limits on frequency range for r-pkt spectrum

extern int ntlcbins; //number of bins for light curve

extern double nusyn_min, nusyn_max; //limits on range for syn
extern int nfake_gam; //# of fake gamma ray lines for syn

/// New variables for other opacity cases, still grey.
extern double opcase3_normal;           ///MK: normalisation factor for opacity_case 3
extern double rho_crit_para;            ///MK: free parameter for the selection of the critical opacity in opacity_case 3
extern double rho_crit;                 ///MK: critical opacity in opacity_case 3 (could now be declared locally)


/// New variables for the non-grey case
//FILE *ldist_file;
extern int debug_packet;                /// activate debug output for this packet if non negative
//double global_threshold;         /// global variable to transfer a threshold frequency to another function
extern int n_middle_it;

extern int total_nlte_levels;            ///total number of nlte levels
extern int n_super_levels;

extern mastate_t *mastate;

extern CELL cell[MGRID+1];
//int *nonemptycells;  /// Array which contains all the non-empty cells cellnumbers
//int nnonemptycells;  /// Total number of non-empty cells



extern bool homogeneous_abundances;

extern bool simulation_continued_from_saved;
extern int nthreads;
extern double nu_rfcut;
extern int n_lte_timesteps;
extern double cell_is_optically_thick;
extern int n_grey_timesteps;
extern int n_titer;
extern bool initial_iteration;
extern int max_bf_continua;
extern int n_kpktdiffusion_timesteps;
extern float kpktdiffusion_timescale;


#endif // GLOBALS_H