#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_integration.h>

#ifdef _OPENMP
  #include "omp.h"
#endif


#define TYPES
#define MGRID 125000  //125000 //1000000 //1000000//262144 //2100000 //125000 //1000000  /* Max number of grid cells.*/
#define MMODELGRID 125000 //125000 //12800 //12800 //125 //3200 //200 //200 //200 //8192 //125 //125000 //200 //125000 //8200 //200 //8200 //200 //125000
#define MPKTS 50000 //62500 //31250 //25000 //40000 //4000 //10000 //10000 //1250 //10000 //100000 //5000 //15625 //15625 /* Maximum number of energy packets in calculation. */
#define MELEMENTS 26 //26 //27 //9
#define MIONS 7 //9
#define MTHREADS 4    /// Max number of OpenMP threads


/// Coolinglist
///============================================================================
/*typedef struct
{
  double contribution;
  short type;
  short level;
  short upperlevel;
  int lineindex;
} ionscoolinglist_t;

typedef struct
{
  short type;
  short element;
  short ion;
  short level;
  short upperlevel;
  int lineindex;
} samplegridcoolinglist_t;*/

typedef struct
{
  double contribution;
  short type;
  short element;
  short ion;
  short level;
  short upperlevel;
  int lineindex;
} cellhistorycoolinglist_t;


typedef struct
{
  double collisional;
  double collbb;
  double collbf;
  double fb;
  double ff;
  double adiabatic;
} coolingrates_t;

typedef struct
{
  double collisional;
  double collbb;
  double collbf;
  double bf;
  double ff;
  double gamma;
} heatingrates_t;


/// PHIXSLIST
///============================================================================
// typedef struct
// {
//   double tau_at_edge;
//   short level;          ///limited to 32767 levels
// } ionsphixslist_t;
//
// typedef struct
// {
//   double nu_edge;
//   short element;
//   short ion;
//   short level;          ///limited to 32767 levels
// } samplegridphixslist_t;


typedef struct
{
  double nu_edge;
  double kappa_bf_contr;
  short element;
  short ion;
  short level;          ///limited to 32767 levels
  int index_in_groundphixslist;
  //double nnlevel;
} fullphixslist_t;


typedef struct
{
  double nu_edge;
  //double photoion_contr;
  double gamma_contr;
  //double stimrecomb_contr;
  //double bfheating_contr;
  short element;
  short ion;
  short level;          ///limited to 32767 levels
} groundphixslist_t;


typedef struct
{
  fullphixslist_t *allcont;
  groundphixslist_t *groundcont;
} phixslist_t;



struct packet
{
  int number;     /// A unique number to identify which packet caused potential troubles.
  int where;      /// The grid cell that the packet is in.
  int type;       /// Identifies the type of packet (k-, r-, etc.)
  double pos[3];  /// Position of the packet (x,y,z).
  double dir[3];  /// Direction of propagation. (x,y,z). Always a unit vector.
  int last_cross; /// To avoid rounding errors on cell crossing.
  double tdecay;  /// Time at which pellet decays.
  double e_cmf;   /// The energy the packet carries in the co-moving frame.
  double e_rf;    /// The energy the packet carries in the rest frame.
  double nu_cmf;  /// The frequency in the co-moving frame.
  double nu_rf;   /// The frequency in the rest frame.

  int escape_type; /// Flag to tell us in which form it escaped from the grid.
  int escape_time; /// Time at which is passes out of the grid.
                   /// Pos, dir, where, e_rf, nu_rf should all remain set at the exit point.

  int scat_count;  /// WHAT'S THAT???

  int next_trans;  /// This keeps track of the next possible line interaction of a rpkt by storing
                   /// its linelist index (to overcome numerical problems in propagating the rpkts).

  int interactions;/// debug: number of interactions the packet undergone
  int last_event;  /// debug: stores information about the packets history

  int emissiontype;   /// records how the packet was emitted if it is a r-pkt
  double em_pos[3]; /// Position of the packet (x,y,z).
  int absorptiontype;     /// records linelistindex of the last absorption
                          /// negative values give ff-abs (-1), bf-abs (-2), compton scattering of gammas (-3),
                          /// photoelectric effect of gammas (-4), pair production of gammas (-5)
                          /// decaying pellets of the 52Fe chain (-6) and pellets which decayed before the
                          /// onset of the simulation (-7)
  double absorptionfreq;  /// records nu_cmf of packet at last absorption
  int nscatterings;   /// records number of electron scatterings a r-pkt undergone since it was emitted
  int em_time;
  double absorptiondir[3]; /// Direction of propagation (x,y,z) when a packet was last absorbed in a line. Always a unit vector.
  //short timestep;
  double stokes[3]; //I, Q and U Stokes parameters
  double pol_dir[3]; //unit vector which defines the coordinate system against which Q and U are measured; should always be perpendicular to dir

};
typedef struct packet PKT;

typedef struct
{
  int element;              /// macro atom of type element (this is an element index)
  int ion;                  /// in ionstage ion (this is an ion index)
  int level;                /// and level=level (this is a level index)
  double nnlevel;           /// population number of the active level
  float statweight;         /// statistical weight of the active level
  int activatedfromlevel;   /// Helper variable for bb-activation of macro atoms due to a rpkt event
                            /// It holds information about the lower level of the bb-transition.
  int activatingline;       /// Linelistindex of the activating line for bb activated MAs, -99 else.
  int lastaction;           /// Holds information on last action performed by do_ma
  double einstein;
} mastate_t;
//mastate_t mastate;


/// GRID
///============================================================================
typedef struct
{
  float abundance;         /// Abundance of the element (by mass!).
  float *groundlevelpop;   /// Pointer to an array of floats which contains the groundlevel populations
                           /// of all included ionisation stages for the element.
  float *partfunct;        /// Pointer to an array of floats which contains the partition functions
                           /// of all included ionisation stages for the element.
  //float *ltepartfunct;     /// Pointer to an array of floats which contains the LTE partition functions
  //                         /// of all included ionisation stages for the element.
} compositionlist_entry;

struct grid
{
  double pos_init[3]; /// Initial co-ordinates of inner most cell.
  int xyz[3];        /// Integer position of cell in grid.
  int modelgridindex;

  //float cen_init[3]; /* Initial co-ordinates of centre of cell. */
  //double wid_init[3];/*Initial width of cell in x, y and z directions.*/
  //  double pos[3]; /* Co-ordinates of inner most corner (x,y,z). */
  //  double cen[3]; /* Centre of cell (x,y,z). */
  //  double wid[3];  /*Width of cell in x, y and z directions.*/

/*  float rho_init;    /// Mass density.
  float rho;         /// Mass density.
  float nne;         /// Number density of free electrons.
  float nne_tot;     /// Number density of all electrons (bound and free ones).

  float f_ni;        /// Initial abundance (by mass) of Ni56 in cell
  float f_co;        /// Initial abundance (by mass) of Co56 in cell
  float f_fe;        /// Abundance (by mass) of Fe-group elements in cell
  float f_ni_stable; /// Initial abundance (by mass) of stable Ni in cell
  float f_co_stable; /// Initial abundance (by mass) of stable Co in cell
  float f_fe_init;   /// Initial abundance (by mass) of Fe56 in cell

  float kappa_grey;  /// Grey opacity for this cell


  float T_e;         /// Electron temperature
  float T_R;         /// Radiation temperature
  float W;           /// Dilution factor
  float nn_tot;*/
  //compositionlist_entry *composition;    /// Pointer to an array which contains the time dependent abundance of all included elements
                                         /// and all the groundlevel populations and partition functions for their ions
};
typedef struct grid CELL;


typedef struct
{
  double *contrib;
} mgicooling_t;

typedef struct
{
  int associated_cells;
  short thick;
  float initial_radial_pos;
  float rhoinit;
  float rho;
  //modelgrid nn_tot
  float nne;
  float nnetot;
  float fni;
  float fco;
  float f52fe;
  float f48cr;
  float ffe;
  float fnistable;
  float fcostable;
  float ffestable;
  float fmnstable;
  float fcrstable;
  float fvstable;
  float ftistable;
  float kappagrey;
  float Te;
  float TR;
  float TJ;
  float W;
  float grey_depth;                      /// Grey optical depth to surface of the modelgridcell
                                         /// This is only stored to print it outside the OpenMP loop in update_grid to the estimatorsfile
                                         /// so there is no need to communicate it via MPI so far!

  compositionlist_entry *composition;    /// Pointer to an array which contains the time dependent abundance of all included elements
                                         /// and all the groundlevel populations and partition functions for their ions
  double totalcooling;
  mgicooling_t *cooling;
} modelgrid_t;






/// Type definition for the sampling grid
/*
#define SAMPLEGRIDSIZE 3
#define SAMPLEGRIDSIZESQUARED 9
#define MSAMPLEGRID 27
typedef struct
{
  float rho;
  float T_R;
  float T_e;
  float W;
  //samplegridcoolinglist_t *coolinglist;
  //samplegridcoolinglist_t *heatinglist;
  samplegridphixslist_t *phixslist;
} samplegrid_t;
*/


#define MGAM_LINES 30 /* Max gamma ray lines per nucleus.*/
struct gamma_ll
{
  int type[3*MGAM_LINES]; /* is it a Ni, Co or fake line */
  int index[3*MGAM_LINES]; /* which of the lines of that element is it */
  int total;         /* the total number of lines in the list */
};
typedef struct gamma_ll LIST;

#define NSYN 1 /* number of frequency points in syn calculation */
struct syn_ray
{
  int last_cross; /* last boundary crossed */
  int where; /* the grid cell the ray is in */
  int status; /*WAITING, then ACTIVE then FINISHED*/
  int lindex[NSYN]; /* array of ray positions in the line list */
  double tstart; /* time at which the ray enters the grid */
  double rstart[3]; /* vector position at which the ray enters the grid */
  double pos[3]; /* current position of the ray */
  double nu_rf[NSYN]; /*rest frame frequencies of rays */
  double nu_rf_init[NSYN]; /*rest frame frequencies of rays at start - for a test only */
  double nu_cmf[NSYN]; /*cmf freqencies */
  double e_rf[NSYN]; /*rest frame energy of rays */
  double e_cmf[NSYN]; /* cmf energies of rays */
};
typedef struct syn_ray RAY;






/// ATOMIC DATA
///============================================================================

typedef struct
{
  short targetlevel;
  double epsilon;
  short stat_weight;
  int lineindex;
} permittedtransitionlist_entry;

typedef struct
{
  float einstein_A;
  float oscillator_strength;
  int linelistindex;
} transitionlist_entry;


typedef struct
{
  double epsilon;                            /// Excitation energy of this level relative to the neutral ground level.
  short stat_weight;                         /// Statistical weight of this level.
  int cont_index;                            /// Index of the continuum associated to this level. Negative number.
  short metastable;                          /// 1 if the level is metastable, else 0

  //double photoion_xs_nu_edge;             /// Number of grid points in the photoion_xs lookup-table.
  float *photoion_xs;         /// Pointer to a lookup-table providing photoionisation cross-sections for this level.

  double *spontrecombcoeff;
//   double *spontrecombcoeff_E;
//   double *photoioncoeff_below;
//   double *photoioncoeff_above;
  double *corrphotoioncoeff;
//  double *corrphotoioncoeff_above;
  double *bfheating_coeff;
//  double *bfheating_coeff_above;
  double *bfcooling_coeff;
  //double *stimulated_bfcooling_coeff;
  //double *stimulated_recomb_coeff;

  //float *modified_spontrecombcoeff;
  //float *stimrecombcoeff;
  //float *modified_stimrecombcoeff;

  ///float population;
  ///float sahafact;
  ///float spontaneousrecomb_ratecoeff;
  ///float modifiedspontaneousrecomb_ratecoeff;
  ///float corrphotoion_ratecoeff;
  ///float modifiedcorrphotoion_ratecoeff;
  ///marates_struct marates;

  /// time dependent macroatom event rates
  //double rad_deexc;                          /// Radiative deexcitation rate from this level.
  //double internal_down_same;                 /// Rate for internal downward transitions within the same ionisation stage.
  //double internal_up_same;                   /// Rate for internal upward transitions within the same ionisation stage.

  //transitionlist_entry *transitions;
  permittedtransitionlist_entry *uptrans;    /// Allowed upward transitions from this level
  permittedtransitionlist_entry *downtrans;  /// Allowed downward transitions from this level
  int closestgroundlevelcont;
} levellist_entry;

typedef struct
{
  short ionstage;                            /// Which ionisation stage: XI=0, XII=1, XIII=2, ...
  int nlevels;                             /// Number of levels for this ionisation stage
  int ionisinglevels;                             /// Number of levels which have a bf-continuum
  int coolingoffset;
  int ncoolingterms;
  double ionpot;                             /// Ionisation threshold to the next ionstage
  //int nbfcontinua;
  //ionsphixslist_t *phixslist;
  float *Alpha_sp;
//  float *zeta;
  levellist_entry *levels;                   /// Carries information for each level: 0,1,...,nlevels-1
} ionlist_entry;

typedef struct
{
  short anumber;                             /// Atomic number
  short nions;                               /// Number of ions for the current element
//  short uppermost_ion;                       /// Highest ionisation stage which has a decent population for a given cell
//                                             /// Be aware that this must not be used outside of the update_grid routine
//                                             /// and their doughters. Neither it will work with OpenMP threads.
  float abundance;                           ///
  float mass;                                /// Atomic mass number in multiple of MH
  ionlist_entry *ions;                       /// Carries information for each ion: 0,1,...,nions-1
} elementlist_entry;



typedef struct
{
  double nu;                                 /// Frequency of the line transition
  float einstein_A;
  float osc_strength;
  short elementindex;                        /// It's a transition of element (not its atomic number,
                                             /// but the (x-1)th element included in the simulation.
  short ionindex;                            /// The same for the elements ion
  short upperlevelindex;                     /// And the participating upper
  short lowerlevelindex;                     /// and lower levels

} linelist_entry;
/// Global pointer to beginning of linelist


typedef struct
{
  short elementindex;
  short ionindex;
  short levelindex;
} bflist_t;


typedef struct
{
  short lower;
  short upper;
  double A;
} transitiontable_entry;  ///only used temporary during input


typedef struct
{
  int cellnumber;
} nne_solution_paras;

typedef struct
{
  int cellnumber;
  double t_current;
} Te_solution_paras;


typedef struct
{
  float T;
  float T2;
  double nu_edge;
} gslintegration_paras;

typedef struct
{
  float T_e;
  int cellnumber;
} gslintegration_ffheatingparas;

typedef struct
{
  double nu_edge;
  int cellnumber;
} gslintegration_bfheatingparas;





typedef struct
{
  double total;
  double es;
  double ff;
  double bf;
  double bf_inrest;
  double ffheating;
  //double bfheating;
} rpkt_cont_opacity_struct;








/// Cell history
///============================================================================
typedef struct
{
  double contribution;
} coolinglist_contributions;

typedef struct
{
  //float population;                      /// The levels population
  double population;                      /// The levels population
  double sahafact;                        ///

  double rad_deexc;                       /// Radiative deexcitation rate from this level.
  double col_deexc;                       /// Collisional deexcitation rate from this level.
  double rad_recomb;                      /// Radiative recombination from this level.
  double col_recomb;                      /// Collisional recombination rate from this level.
  double internal_down_same;              /// Rate for internal downward transitions to same ionisation stage.
  double internal_up_same;                /// Rate for internal upward transitions to same ionisation stage.
  double internal_down_lower;             /// Rate for internal downward transitions to lower ionisation stage.
  double internal_up_higher;              /// Rate for internal upward transitions to higher ionisation stage.

  double *individ_rad_deexc;
  double *individ_internal_down_same;
  double *individ_internal_up_same;

  double spontaneousrecombrate;
  double bfcooling;
  double corrphotoioncoeff;
} chlevels_struct;

typedef struct
{
  chlevels_struct *chlevels;              /// Pointer to the ions levellist.
} chions_struct;

typedef struct
{
  chions_struct *chions;                  /// Pointer to the elements ionlist.
} chelements_struct;


typedef struct
{
  int cellnumber;                         /// Identifies the cell the data is valid for.
  //double totalcooling;                    /// Total cooling rate in this cell.
  //double bfcooling;                    /// Total cooling rate in this cell.
//  coolinglist_contributions *coolinglist; /// Cooling contributions by the different processes.
  cellhistorycoolinglist_t *coolinglist; /// Cooling contributions by the different processes.
  chelements_struct *chelements;          /// Pointer to a nested list which helds compositional
                                          /// information for all the elements=0,1,...,nelements-1
} cellhistory_struct;



typedef struct
{
  short *to;
} transitions_t;

