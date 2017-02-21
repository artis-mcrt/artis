#ifndef TYPES_H
#define TYPES_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
//#include <gsl/gsl_sf_expint.h>

#define MGRID  125000 //125000 //1000000 //1000000//262144 //2100000 //125000 //1000000  /* Max number of grid cells.*/
#define MMODELGRID 178 //125000 //12800 //12800 //125 //3200 //200 //200 //200 //8192 //125 //125000 //200 //125000 //8200 //200 //8200 //200 //125000
#define MPKTS 600000 //25000 //40000 //4000 //10000 //10000 //1250 //10000 //100000 //5000 //15625 //15625 /* Maximum number of energy packets in calculation. */
#define MELEMENTS 26 //26 //27 //9
#define MIONS 5 //9
#define MTHREADS 8    /// Max number of OpenMP threads


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

enum coolingtype {
  COOLINGTYPE_FF         = 880,
  COOLINGTYPE_FB         = 881,
  COOLINGTYPE_COLLEXC    = 882,
  COOLINGTYPE_COLLION    = 883,
};

typedef struct cellhistorycoolinglist_t
{
  double contribution;
  enum coolingtype type;
  int element;
  int ion;
  int level;
  int upperlevel;
  int lineindex;
} cellhistorycoolinglist_t;

/*enum heatingtype {
  HEATINGTYPE_FF         = 884,
  HEATINGTYPE_BF         = 885,
  HEATINGTYPE_COLLDEEXC  = 886,
  HEATINGTYPE_COLLRECOMB = 887,
};*/

typedef struct coolingrates_t
{
  double collisional;
  //double collbb;
  //double collbf;
  double fb;
  double ff;
  double adiabatic;
} coolingrates_t;

typedef struct heatingrates_t
{
  double collisional;
  //double collbb;
  //double collbf;
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


typedef struct fullphixslist_t
{
  double nu_edge;
  double kappa_bf_contr;
  int element;
  int ion;
  int level;          ///limited to 32767 levels
  int index_in_groundphixslist;
  //double nnlevel;
} fullphixslist_t;

typedef struct groundphixslist_t
{
  double nu_edge;
  //double photoion_contr;
  double gamma_contr;
  //double stimrecomb_contr;
  //double bfheating_contr;
  int element;
  int ion;
  int level;          ///limited to 32767 levels
} groundphixslist_t;


typedef struct phixslist_t
{
  fullphixslist_t *restrict allcont;
  groundphixslist_t *restrict groundcont;
} phixslist_t;

enum packet_type {
  TYPE_ESCAPE = 32,
  TYPE_NICKEL_PELLET = 100,
  TYPE_COBALT_PELLET = 101,
  TYPE_48CR_PELLET = 102,
  TYPE_48V_PELLET = 103,
  TYPE_52FE_PELLET = 104,
  TYPE_52MN_PELLET = 105,
  TYPE_COBALT_POSITRON_PELLET = 106,
  TYPE_GAMMA = 10,
  TYPE_RPKT = 11,
  TYPE_KPKT = 12,
  TYPE_MA = 13,
  TYPE_EMINUS = 20,
  TYPE_PRE_KPKT = 120,
  TYPE_GAMMA_KPKT = 121,
};

enum cell_boundary {
  NEG_X = 101,
  POS_X = 102,
  NEG_Y = 103,
  POS_Y = 104,
  NEG_Z = 105,
  POS_Z = 106,
  NONE = 107,
};

typedef struct packet
{
  int where;      /// The grid cell that the packet is in.
  enum packet_type type;       /// Identifies the type of packet (k-, r-, etc.)
  enum cell_boundary last_cross; /// To avoid rounding errors on cell crossing.
  int interactions;/// debug: number of interactions the packet undergone
  int nscatterings;   /// records number of electron scatterings a r-pkt undergone since it was emitted
  int last_event;  /// debug: stores information about the packets history
  double pos[3];  /// Position of the packet (x,y,z).
  double dir[3];  /// Direction of propagation. (x,y,z). Always a unit vector.
  double e_cmf;   /// The energy the packet carries in the co-moving frame.
  double e_rf;    /// The energy the packet carries in the rest frame.
  double nu_cmf;  /// The frequency in the co-moving frame.
  double nu_rf;   /// The frequency in the rest frame.
  int next_trans;  /// This keeps track of the next possible line interaction of a rpkt by storing
                   /// its linelist index (to overcome numerical problems in propagating the rpkts).
  int emissiontype;   /// records how the packet was emitted if it is a r-pkt
  double em_pos[3];   /// Position of the packet (x,y,z).
  int em_time;
  int absorptiontype;     /// records linelistindex of the last absorption
                          /// negative values give ff-abs (-1), bf-abs (-2), compton scattering of gammas (-3),
                          /// photoelectric effect of gammas (-4), pair production of gammas (-5)
                          /// decaying pellets of the 52Fe chain (-6) and pellets which decayed before the
                          /// onset of the simulation (-7)
  double absorptionfreq;  /// records nu_cmf of packet at last absorption
  double absorptiondir[3]; /// Direction of propagation (x,y,z) when a packet was last absorbed in a line. Always a unit vector.
  //short timestep;
  double stokes[3]; //I, Q and U Stokes parameters
  double stokes_qu[2]; //Q and U Stokes parameters
  double pol_dir[3]; //unit vector which defines the coordinate system against which Q and U are measured; should always be perpendicular to dir
  enum packet_type escape_type; /// Flag to tell us in which form it escaped from the grid.
  double tdecay;  /// Time at which pellet decays.
  int escape_time; /// Time at which is passes out of the grid.
                   /// Pos, dir, where, e_rf, nu_rf should all remain set at the exit point.
  int scat_count;  /// WHAT'S THAT???
  int number;     /// A unique number to identify which packet caused potential troubles.
} PKT;

enum ma_action {
  MA_ACTION_NONE = 0,
  MA_ACTION_RADDEEXC = 1,
  MA_ACTION_COLDEEXC = 2,
  MA_ACTION_RADRECOMB = 3,
  MA_ACTION_COLRECOMB = 4,
  MA_ACTION_INTERNALDOWNSAME = 5,
  MA_ACTION_INTERNALDOWNLOWER = 6,
  MA_ACTION_INTERNALUPSAME = 7,
  MA_ACTION_INTERNALUPHIGHER = 8,
};

typedef struct mastate_t
{
  double nnlevel;           /// population number of the active level
  double einstein;
  // double statweight;         /// statistical weight of the active level
  int element;              /// macro atom of type element (this is an element index)
  int ion;                  /// in ionstage ion (this is an ion index)
  int level;                /// and level=level (this is a level index)
  int activatedfromlevel;   /// Helper variable for bb-activation of macro atoms due to a rpkt event
                            /// It holds information about the lower level of the bb-transition.
  int activatingline;       /// Linelistindex of the activating line for bb activated MAs, -99 else.
  enum ma_action lastaction;           /// Holds information on last action performed by do_ma
} mastate_t;


/// GRID
///============================================================================
typedef struct compositionlist_entry
{
  float abundance;         /// Abundance of the element (by mass!).
  float *groundlevelpop;   /// Pointer to an array of floats which contains the groundlevel populations
                           /// of all included ionisation stages for the element.
  float *partfunct;        /// Pointer to an array of floats which contains the partition functions
                           /// of all included ionisation stages for the element.
  //float *ltepartfunct;     /// Pointer to an array of floats which contains the LTE partition functions
  //                         /// of all included ionisation stages for the element.
} compositionlist_entry;

typedef struct grid
{
  double pos_init[3]; /// Initial co-ordinates of inner most corner of cell.
  int xyz[3];         /// Integer position of cell in grid.
  int modelgridindex;
} CELL;

typedef struct mgicooling_t
{
  double *contrib;
} mgicooling_t;


typedef struct modelgrid_t
{
  double Te;
  double TR;
  double TJ;
  double W;
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
  float grey_depth;                      /// Grey optical depth to surface of the modelgridcell
                                         /// This is only stored to print it outside the OpenMP loop in update_grid to the estimatorsfile
                                         /// so there is no need to communicate it via MPI so far!

  compositionlist_entry *composition;    /// Pointer to an array which contains the time dependent abundance of all included elements
                                         /// and all the groundlevel
                                         /// populations and partition
                                         /// functions for their ions
  double *nlte_pops;                      /// Pointer to an array that
					 /// contains the nlte-level
					 /// populations for this cell

  double totalcooling;
  mgicooling_t *cooling;
  short thick;
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
typedef struct gamma_ll
{
  int type[3 * MGAM_LINES]; /* is it a Ni, Co or fake line */
  int index[3 * MGAM_LINES]; /* which of the lines of that element is it */
  int total;         /* the total number of lines in the list */
} LIST;

#define NSYN 1 /* number of frequency points in syn calculation */

enum ray_status {
  ACTIVE = 1,
  WAITING = 2,
  FINISHED = 3,
};

typedef struct syn_ray
{
  double tstart; /* time at which the ray enters the grid */
  double rstart[3]; /* vector position at which the ray enters the grid */
  double pos[3]; /* current position of the ray */
  double nu_rf[NSYN]; /*rest frame frequencies of rays */
  double nu_rf_init[NSYN]; /*rest frame frequencies of rays at start - for a test only */
  double nu_cmf[NSYN]; /*cmf freqencies */
  double e_rf[NSYN]; /*rest frame energy of rays */
  double e_cmf[NSYN]; /* cmf energies of rays */
  int where; /* the grid cell the ray is in */
  int lindex[NSYN]; /* array of ray positions in the line list */
  enum cell_boundary last_cross; /* last boundary crossed */
  enum ray_status status; /*WAITING, then ACTIVE then FINISHED*/
} RAY;



/// ATOMIC DATA
///============================================================================

typedef struct transitionlist_entry
{
  double epsilon;
  int targetlevel;
  int lineindex;
  int stat_weight;
} transitionlist_entry;

/*
typedef struct
{
  float einstein_A;
  float oscillator_strength;
  int linelistindex;
} transitionlist_entry;
*/

typedef struct phixstarget_entry
{
  double *restrict spontrecombcoeff;
  double *restrict corrphotoioncoeff;
  double *restrict bfheating_coeff;
  double *restrict bfcooling_coeff;

  double probability;        // fraction of phixs cross section leading to this final level
  int levelindex;         // index of upper ion level after photoionisation
} phixstarget_entry;

typedef struct levellist_entry
{
  transitionlist_entry *restrict uptrans;    /// Allowed upward transitions from this level
  transitionlist_entry *restrict downtrans;  /// Allowed downward transitions from this level

  double epsilon;                            /// Excitation energy of this level relative to the neutral ground level.
  phixstarget_entry *restrict phixstargets;  /// pointer to table of target states and probabilities
  float *restrict photoion_xs;               /// Pointer to a lookup-table providing photoionisation cross-sections for this level.
  int nphixstargets;                         /// length of phixstargets array:
  int stat_weight;                           /// Statistical weight of this level.

  int cont_index;                            /// Index of the continuum associated to this level. Negative number.
  int closestgroundlevelcont;

  bool is_nlte;                              /// 1 if the level is to
                                             /// be treated in nlte
  bool metastable;                            ///

//  double photoion_xs_nu_edge;              /// nu of the first grid point in the photoion_xs lookup-table.

//double *spontrecombcoeff_E;
//double *photoioncoeff_below;
//double *photoioncoeff_above;
//double *corrphotoioncoeff_above;
//double *bfheating_coeff_above;
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

/// time dependent macroatom event rates
//double rad_deexc;                          /// Radiative deexcitation rate from this level.
//double internal_down_same;                 /// Rate for internal downward transitions within the same ionisation stage.
//double internal_up_same;                   /// Rate for internal upward transitions within the same ionisation stage.

} levellist_entry;

typedef struct ionlist_entry
{
  levellist_entry *restrict levels;                   /// Carries information for each level: 0,1,...,nlevels-1
  int ionstage;                            /// Which ionisation stage: XI=0, XII=1, XIII=2, ...
  int nlevels;                               /// Number of levels for this ionisation stage
  int nlevels_nlte;                          /// number of nlte levels for this ion
  int first_nlte;                            /// reference index for counting of nlte levels
  int ionisinglevels;                        /// Number of levels which have a bf-continuum
  int coolingoffset;
  int ncoolingterms;
  double ionpot;                             /// Ionisation threshold to the next ionstage
  double *Alpha_sp;
  //int nbfcontinua;
  //ionsphixslist_t *phixslist;
} ionlist_entry;

typedef struct elementlist_entry
{
  ionlist_entry *restrict ions;                       /// Carries information for each ion: 0,1,...,nions-1
  int nions;                               /// Number of ions for the current element
  int anumber;                             /// Atomic number
//  int uppermost_ion;                       /// Highest ionisation stage which has a decent population for a given cell
//                                             /// Be aware that this must not be used outside of the update_grid routine
//                                             /// and their daughters. Neither it will work with OpenMP threads.
  float abundance;                           ///
  float mass;                                /// Atomic mass number in multiple of MH
} elementlist_entry;

typedef struct linelist_entry
{
  double nu;                                 /// Frequency of the line transition
  float einstein_A;
  float osc_strength;
  float coll_str;
  int elementindex;                        /// It's a transition of element (not its atomic number,
                                           /// but the (x-1)th element included in the simulation.
  int ionindex;                            /// The same for the elements ion
  int upperlevelindex;                     /// And the participating upper
  int lowerlevelindex;                     /// and lower levels
  bool forbidden;
} linelist_entry;

typedef struct bflist_t
{
  int elementindex;
  int ionindex;
  int levelindex;
} bflist_t;

typedef struct nne_solution_paras
{
  int cellnumber;
} nne_solution_paras;

typedef struct gslintegration_paras
{
  double nu_edge;
  double T;
} gslintegration_paras;

typedef struct rpkt_cont_opacity_struct
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
typedef struct coolinglist_contributions
{
  double contribution;
} coolinglist_contributions;

typedef struct
{
  double spontaneousrecombrate;
  double bfcooling;
  double bfheatingcoeff;
  double corrphotoioncoeff;
  double sahafact;
} chphixstargets_struct;


typedef struct chlevels_struct
{
  double population;                      /// The level's population

  double rad_deexc;                       /// Radiative deexcitation rate from this level.
  double col_deexc;                       /// Collisional deexcitation rate from this level.
  double rad_recomb;                      /// Radiative recombination from this level.
  double col_recomb;                      /// Collisional recombination rate from this level.
  double internal_down_same;              /// Rate for internal downward transitions to same ionisation stage.
  double internal_up_same;                /// Rate for internal upward transitions to same ionisation stage.
  double internal_down_lower;             /// Rate for internal downward transitions to lower ionisation stage.
  double internal_up_higher;              /// Rate for internal upward transitions to higher ionisation stage.

  double *restrict individ_rad_deexc;
  double *restrict individ_internal_down_same;
  double *restrict individ_internal_up_same;

  chphixstargets_struct *restrict chphixstargets;
} chlevels_struct;

typedef struct chions_struct
{
  chlevels_struct *restrict chlevels;              /// Pointer to the ions levellist.
} chions_struct;

typedef struct chelements_struct
{
  chions_struct *restrict chions;                  /// Pointer to the elements ionlist.
} chelements_struct;

typedef struct cellhistory_struct
{
//  double totalcooling;                    /// Total cooling rate in this cell.
//  double bfcooling;                       /// Total cooling rate in this cell.
//  coolinglist_contributions *coolinglist; /// Cooling contributions by the different processes.
  cellhistorycoolinglist_t *restrict coolinglist;    /// Cooling contributions by the different processes.
  chelements_struct *restrict chelements;            /// Pointer to a nested list which helds compositional
                                            /// information for all the elements=0,1,...,nelements-1
  int cellnumber;                           /// Identifies the cell the data is valid for.
} cellhistory_struct;


typedef struct transitions_t
{
  int *restrict to;
} transitions_t;

#endif //TYPES_H
