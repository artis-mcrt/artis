#ifndef TYPES_H
#define TYPES_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
//#include <gsl/gsl_sf_expint.h>

// Max number of propagation grid cells
//125000 //1000000 //1000000//262144 //2100000 //125000 //1000000
#define MGRID  125000

// Max number of input grid cells
//125000 //12800 //12800 //125 //3200 //200 //200 //200 //8192 //125 //125000 //200 //125000 //8200 //200 //8200 //200 //125000
#define MMODELGRID 8200

// Maximum number of energy packets in calculation.
//25000 //40000 //4000 //10000 //10000 //1250 //10000 //100000 //5000 //15625 //15625
#define MPKTS 2000000

//26 //27 //9
#define MELEMENTS 17

// 9
#define MIONS 5

// Max number of OpenMP threads
#define MTHREADS 8


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
  double nt_frac_heating; // = heatingrates.gamma / total_gamma_deposition, = get_nt_frac_heating(modelgridindex) when T_e solver runs
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
  TYPE_56NI_PELLET = 100,
  TYPE_56CO_PELLET = 101,
  TYPE_48CR_PELLET = 102,
  TYPE_48V_PELLET = 103,
  TYPE_52FE_PELLET = 104,
  TYPE_52MN_PELLET = 105,
  TYPE_56CO_POSITRON_PELLET = 106,
  TYPE_57NI_PELLET = 107,
  TYPE_57NI_POSITRON_PELLET = 108,
  TYPE_57CO_PELLET = 109,
  TYPE_GAMMA = 10,
  TYPE_RPKT = 11,
  TYPE_KPKT = 12,
  TYPE_MA = 13,
  TYPE_NTLEPTON = 20,
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
                          /// decay of a positron pellet (-10)
  int trueemissiontype;  // emission type coming from a kpkt to rpkt (last thermal emission)
  double absorptionfreq;  /// records nu_cmf of packet at last absorption
  double absorptiondir[3]; /// Direction of propagation (x,y,z) when a packet was last absorbed in a line. Always a unit vector.
  //short timestep;
  double stokes[3]; //I, Q and U Stokes parameters
  double pol_dir[3]; //unit vector which defines the coordinate system against which Q and U are measured; should always be perpendicular to dir
  double tdecay;  /// Time at which pellet decays.
  enum packet_type escape_type; /// Flag to tell us in which form it escaped from the grid.
  int escape_time; /// Time at which is passes out of the grid.
                   /// Pos, dir, where, e_rf, nu_rf should all remain set at the exit point.
  int scat_count;  /// WHAT'S THAT???
  int number;     /// A unique number to identify which packet caused potential troubles.
  bool originated_from_positron; // first-non-pellet packet type was positron
  float trueemissionvelocity;
} PKT;

enum ma_action {
  /// Radiative deexcitation rate from this level.
  MA_ACTION_RADDEEXC = 0,
  /// Collisional deexcitation rate from this level.
  MA_ACTION_COLDEEXC = 1,
  /// Radiative recombination from this level.
  MA_ACTION_RADRECOMB = 2,
  /// Collisional recombination rate from this level.
  MA_ACTION_COLRECOMB = 3,
  /// Rate for internal downward transitions to same ionisation stage.
  MA_ACTION_INTERNALDOWNSAME = 4,
  /// Rate for internal upward transitions to same ionisation stage.
  MA_ACTION_INTERNALDOWNLOWER = 5,
  /// Rate for internal downward transitions to lower ionisation stage.
  MA_ACTION_INTERNALUPSAME = 6,
  /// Rate for internal upward transitions to higher ionisation stage.
  MA_ACTION_INTERNALUPHIGHER = 7,
  /// Rate for internal upward transitions to higher ionisation stage due to non-thermal collisions.
  MA_ACTION_INTERNALUPHIGHERNT = 8,
  MA_ACTION_COUNT = 9,
};

typedef struct mastate_t
{
  int element;              /// macro atom of type element (this is an element index)
  int ion;                  /// in ionstage ion (this is an ion index)
  int level;                /// and level=level (this is a level index)
  int activatingline;       /// Linelistindex of the activating line for bb activated MAs, -99 else.
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
  // int xyz[3];         /// Integer position of cell in grid.
  int modelgridindex;
} CELL;

typedef struct mgicooling_t
{
  double *contrib;
} mgicooling_t;

enum radionuclides {
  NUCLIDE_NI57 = 0,
  NUCLIDE_NI56 = 1,
  NUCLIDE_CO56 = 2,
  FAKE_GAM_LINE_ID = 3,
  NUCLIDE_CR48 = 4,
  NUCLIDE_V48 = 5,
  NUCLIDE_CO57 = 6,
  NUCLIDE_FE52 = 7,
  RADIONUCLIDE_COUNT = 8,
};

typedef struct gamma_ll
{
  enum radionuclides *nuclidetype; // is it a Ni56, Co56, a fake line, etc
  int *index;               // which of the lines of that element is it */
  int total;                // the total number of lines in the list */
} LIST;

typedef struct modelgrid_t
{
  float Te;
  float TR;
  float TJ;
  float W;
  float nne;
  float initial_radial_pos;
  float rhoinit;
  float rho;
  //modelgrid nn_tot
  float nnetot;           // total electron density (free + bound).
  float fradionuclides[RADIONUCLIDE_COUNT];
  float ffegrp;
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

  compositionlist_entry *composition;    /// Pointer to an array which contains the time dependent abundances
                                        /// of all included elements and all the groundlevel
                                         /// populations and partition functions for their ions
  double *nlte_pops;                     /// Pointer to an array that contains the nlte-level
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

typedef struct phixstarget_entry
{
  double *restrict spontrecombcoeff;
  #if (!NO_LUT_PHOTOION)
    double *restrict corrphotoioncoeff;
  #endif
  #if (!NO_LUT_BFHEATING)
  double *restrict bfheating_coeff;
    #endif
  double *restrict bfcooling_coeff;

  double probability;        // fraction of phixs cross section leading to this final level
  int levelindex;         // index of upper ion level after photoionisation
} phixstarget_entry;


typedef struct levellist_entry
{
  double epsilon;                            /// Excitation energy of this level relative to the neutral ground level.
  int *uptrans_lineindicies;    /// Allowed upward transitions from this level
  int *downtrans_lineindicies;  /// Allowed downward transitions from this level
  int nuptrans;
  int ndowntrans;
  double phixs_threshold;                    /// Energy of first point in the photion_xs table
  phixstarget_entry *restrict phixstargets;  /// pointer to table of target states and probabilities
  float *restrict photoion_xs;               /// Pointer to a lookup-table providing photoionisation cross-sections for this level.
  int nphixstargets;                         /// length of phixstargets array:
  int stat_weight;                           /// Statistical weight of this level.

  int cont_index;                            /// Index of the continuum associated to this level. Negative number.
#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
  int closestgroundlevelcont;
#endif

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
  levellist_entry *levels;                   /// Carries information for each level: 0,1,...,nlevels-1
  int ionstage;                              /// Which ionisation stage: XI=0, XII=1, XIII=2, ...
  int nlevels;                               /// Number of levels for this ionisation stage
  int nlevels_nlte;                          /// number of nlte levels for this ion
  int first_nlte;                            /// index into nlte_pops array of a grid cell
  int ionisinglevels;                        /// Number of levels which have a bf-continuum
  int maxrecombininglevel;                   /// level index of the highest level with a non-zero recombination rate
  int coolingoffset;
  int ncoolingterms;
  float *Alpha_sp;
  double ionpot;                             /// Ionisation threshold to the next ionstage
  //int nbfcontinua;
  //ionsphixslist_t *phixslist;
} ionlist_entry;

typedef struct elementlist_entry
{
  ionlist_entry *restrict ions;              /// Carries information for each ion: 0,1,...,nions-1
  int nions;                                 /// Number of ions for the current element
  int anumber;                               /// Atomic number
//  int uppermost_ion;                       /// Highest ionisation stage which has a decent population for a given cell
                                             /// Be aware that this must not be used outside of the update_grid routine
                                             /// and their daughters. Neither it will work with OpenMP threads.
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
// typedef struct coolinglist_contributions
// {
//   double contribution;
// } coolinglist_contributions;

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
  double processrates[MA_ACTION_COUNT];

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


#endif //TYPES_H
