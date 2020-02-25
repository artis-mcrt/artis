#ifndef GLOBALS_H
#define GLOBALS_H

#ifndef CUDA_ENABLED
#define CUDA_ENABLED false
#endif

#if CUDA_ENABLED
#include <curand_kernel.h>
#endif

#include "types.h"
#include "sn3d.h"

extern __managed__ double syn_dir[3]; // vector pointing from origin to observer

extern __managed__ struct time time_step[MTSTEP];

#define MSYN_TIME 100
extern __managed__ int nsyn_time;
extern __managed__ double time_syn[MSYN_TIME];

#define EMISS_MAX 2
extern __managed__ int emiss_offset;
extern __managed__ int emiss_max;

#if CUDA_ENABLED
extern __managed__ curandState curandstates[MTHREADS];
#endif

extern __managed__ modelgrid_t modelgrid[MMODELGRID + 1];

extern __managed__ float compton_emiss[MMODELGRID+1][EMISS_MAX];
extern __managed__ double rpkt_emiss[MMODELGRID+1];


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
  extern __managed__ double ffheatingestimator[MMODELGRID + 1];
  extern __managed__ double colheatingestimator[MMODELGRID + 1];

  #ifdef DO_TITER
    extern double ffheatingestimator_save[MMODELGRID];
    extern double colheatingestimator_save[MMODELGRID];
    extern double gammaestimator_save[MMODELGRID * MELEMENTS * MIONS];
    extern double bfheatingestimator_save[MMODELGRID * MELEMENTS * MIONS];
  #endif
#endif

#ifdef RECORD_LINESTAT
  extern __managed__ int *ecounter;
  extern __managed__ int *acounter;
  extern __managed__ int *linestat_reduced;
#endif


extern bool file_set;

extern __managed__ bool do_comp_est;
extern __managed__ bool do_r_lc;
extern __managed__ int do_rlc_est;

extern __managed__ int n_out_it;

extern __managed__ int npts_model;
extern __managed__ double vout_model[MMODELGRID];
extern __managed__ double t_model;
extern __managed__ int ncoord1_model;
extern __managed__ int ncoord2_model;
extern __managed__ double dcoord1;
extern __managed__ double dcoord2;

extern __managed__ double CLIGHT_PROP;

extern __managed__ double gamma_grey;

extern __managed__ double min_den;

#define GREY_OP 0.1

extern __managed__ double max_path_step;

extern __managed__ int opacity_case;

extern __managed__ double dlogt;
extern __managed__ int maxion;
extern __managed__ short elements_uppermost_ion[MTHREADS][MELEMENTS];

extern __managed__ int nelements;
extern __managed__ int nlines;
extern __managed__ int includedions;
extern __managed__ elementlist_entry *elements;
extern __managed__ linelist_entry *linelist;
extern __managed__ bflist_t *bflist;

extern __managed__ rpkt_cont_opacity_struct *kappa_rpkt_cont;

extern __managed__ int ncoolingterms;
extern __managed__ int importantcoolingterms;

extern __managed__ fullphixslist_t *phixsallcont;
extern __managed__ groundphixslist_t *phixsgroundcont;

extern __managed__ int nbfcontinua;
extern __managed__ int nbfcontinua_ground;
extern __managed__ int NPHIXSPOINTS;
extern __managed__ double NPHIXSNUINCREMENT;

extern __managed__ cellhistory_struct *cellhistory;

extern __managed__ int ma_stat_activation_collexc;
extern __managed__ int ma_stat_activation_collion;
extern __managed__ int ma_stat_activation_ntcollexc;
extern __managed__ int ma_stat_activation_ntcollion;
extern __managed__ int ma_stat_activation_bb;
extern __managed__ int ma_stat_activation_bf;
extern __managed__ int ma_stat_activation_fb;
extern __managed__ int ma_stat_deactivation_colldeexc;
extern __managed__ int ma_stat_deactivation_collrecomb;
extern __managed__ int ma_stat_deactivation_bb;
extern __managed__ int ma_stat_deactivation_fb;
extern __managed__ int ma_stat_internaluphigher;
extern __managed__ int ma_stat_internaluphighernt;
extern __managed__ int ma_stat_internaldownlower;
extern __managed__ int k_stat_to_ma_collexc;
extern __managed__ int k_stat_to_ma_collion;
extern __managed__ int k_stat_to_r_ff;
extern __managed__ int k_stat_to_r_fb;
extern __managed__ int k_stat_to_r_bb;
extern __managed__ int k_stat_from_ff;
extern __managed__ int k_stat_from_bf;
extern __managed__ int nt_stat_from_gamma;
extern __managed__ int k_stat_from_earlierdecay;
extern __managed__ int escounter;
extern __managed__ int resonancescatterings;
extern __managed__ int cellcrossings;
extern __managed__ int upscatter;
extern __managed__ int downscatter;
extern __managed__ int updatecellcounter;
extern __managed__ int coolingratecalccounter;

extern __managed__ int debuglevel;

extern __managed__ int ncoordgrid[3];
extern __managed__ int ngrid;
extern __managed__ int grid_type;
extern __managed__ char coordlabel[3];

extern __managed__ enum model_types model_type;

extern __managed__ int nprocs;
extern __managed__ int rank_global;
extern __managed__ int npkts;
extern __managed__ int nesc;

extern __managed__ double coordmax[3];
extern __managed__ double mtot;
extern __managed__ double vmax;
extern __managed__ double rmax;
extern __managed__ double totmassradionuclide[RADIONUCLIDE_COUNT];
extern __managed__ double mfeg;
extern __managed__ double tmax;
extern __managed__ double tmin;

extern __managed__ int ntstep;
extern __managed__ int itstep;
extern __managed__ int ftstep;
extern __managed__ int nts_global;

__managed__ extern int ntbins;
__managed__ extern int nnubins;
__managed__ extern double nu_min_r;
__managed__ extern double nu_max_r;

extern int ntlcbins;

extern double nusyn_min, nusyn_max;
extern int nfake_gam;

extern __managed__ double opcase3_normal;
extern __managed__ double rho_crit_para;
extern __managed__ double rho_crit;

extern __managed__ int debug_packet;
extern __managed__ int n_middle_it;

extern __managed__ int total_nlte_levels;
extern __managed__ int n_super_levels;

extern __managed__ mastate_t *mastate;

extern __managed__ CELL cell[MGRID+1];


extern __managed__ bool homogeneous_abundances;

extern __managed__ bool simulation_continued_from_saved;
extern __managed__ int nthreads;
extern __managed__ double nu_rfcut;
extern __managed__ int n_lte_timesteps;
extern __managed__ double cell_is_optically_thick;
extern __managed__ int n_grey_timesteps;
extern __managed__ int n_titer;
extern __managed__ bool initial_iteration;
extern __managed__ int max_bf_continua;
extern __managed__ int n_kpktdiffusion_timesteps;
extern __managed__ float kpktdiffusion_timescale;


#endif // GLOBALS_H