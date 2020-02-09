#ifndef GLOBALS_H
#define GLOBALS_H

#include "types.h"
#include "sn3d.h"
#include "globals.h"

extern double syn_dir[3]; // vector pointing from origin to observer

extern struct time time_step[MTSTEP];

#define MSYN_TIME 100
extern int nsyn_time;
extern double time_syn[MSYN_TIME];

#define EMISS_MAX 2
extern int emiss_offset;
extern int emiss_max;


extern __managed__ modelgrid_t modelgrid[MMODELGRID + 1];

extern float compton_emiss[MMODELGRID+1][EMISS_MAX];
extern double rpkt_emiss[MMODELGRID+1];


#if (!NO_LUT_PHOTOION)
  extern double corrphotoionrenorm[MMODELGRID * MELEMENTS * MIONS];
  extern double gammaestimator[MMODELGRID * MELEMENTS * MIONS];
#endif
#if (!NO_LUT_BFHEATING)
  extern double bfheatingestimator[MMODELGRID * MELEMENTS * MIONS];
#endif
#ifdef FORCE_LTE
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
#endif

#ifdef RECORD_LINESTAT
  extern int *ecounter;
  extern int *acounter;
  extern int *linestat_reduced;
#endif


extern bool file_set;

extern bool do_comp_est;
extern bool do_r_lc;
extern int do_rlc_est;

extern int n_out_it;

extern int npts_model;
extern double vout_model[MMODELGRID];
extern double t_model;
extern int ncoord1_model;
extern int ncoord2_model;
extern double dcoord1, dcoord2;

extern double CLIGHT_PROP;

extern double gamma_grey;

extern double min_den;

#define GREY_OP 0.1

extern double max_path_step;

extern int opacity_case;

extern double dlogt;
extern int maxion;
extern short elements_uppermost_ion[MTHREADS][MELEMENTS];

extern __managed__ int nelements;
extern __managed__ int nlines;
extern __managed__ int includedions;
extern __managed__ elementlist_entry *elements;
extern __managed__ linelist_entry *linelist;
extern __managed__ bflist_t *bflist;

extern __managed__ rpkt_cont_opacity_struct *kappa_rpkt_cont;

extern int ncoolingterms;
extern int importantcoolingterms;

extern __managed__ phixslist_t *phixslist;
extern __managed__ int nbfcontinua;
extern __managed__ int nbfcontinua_ground;
extern __managed__ int NPHIXSPOINTS;
extern __managed__ double NPHIXSNUINCREMENT;

extern cellhistory_struct *cellhistory;

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

extern int ncoordgrid[3];
extern int ngrid;
extern int grid_type;
extern char coordlabel[3];

extern enum model_types model_type;

extern int nprocs;
extern int rank_global;
extern int npkts;
extern int nesc;

extern double coordmax[3];
extern double mtot;
extern double vmax;
extern double rmax;
extern double totmassradionuclide[RADIONUCLIDE_COUNT];
extern double mfeg;
extern double tmax;
extern double tmin;

extern int ntstep;
extern int itstep;
extern int ftstep;
extern __managed__ int nts_global;

extern int ntbins, nnubins;
extern double nu_min_r, nu_max_r;

extern int ntlcbins;

extern double nusyn_min, nusyn_max;
extern int nfake_gam;

extern double opcase3_normal;
extern double rho_crit_para;
extern double rho_crit;

extern int debug_packet;
extern int n_middle_it;

extern int total_nlte_levels;
extern int n_super_levels;

extern mastate_t *mastate;

extern CELL cell[MGRID+1];


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