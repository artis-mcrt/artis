#include "constants.h"
#include "artisoptions.h"

#ifndef TYPES_H
#define TYPES_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#ifndef __CUDA_ARCH__
#include <gsl/gsl_rng.h>
#endif
//#include <gsl/gsl_sf_expint.h>

struct time
{
  double start; // time at start of this timestep. [s]
  double width; // Width of timestep. [s]
  double mid; // Mid time in step - computed logarithmically. [s]
  double gamma_dep; // cmf gamma ray energy deposition from absorption events [erg]
  double gamma_dep_pathint; // cmf gamma ray energy deposition from packet trajectories [erg]
  double positron_dep; // cmf positron energy deposition [erg]
  double eps_positron_ana_power; // cmf positron KE energy generation rate analytical [erg/s]
  double electron_dep; // cmf electron energy deposition [erg]
  double electron_emission; // cmf electron KE energy generation [erg]
  double eps_electron_ana_power; // cmf electron KE energy generation rate analytical [erg/s]
  double alpha_dep; // cmf alpha energy deposition [erg]
  double alpha_emission; // cmf alpha KE energy generation [erg]
  double eps_alpha_ana_power; // cmf alpha KE energy generation rate analytical [erg/s]
  double gamma_emission; // gamma decay energy generation in this timestep [erg]
  double qdot_betaminus; // energy generation from beta-minus decays (including neutrinos) [erg/s/g]
  double qdot_alpha;  // energy generation from alpha decays (including neutrinos) [erg/s/g]
  double qdot_total;  // energy generation from all decays (including neutrinos) [erg/s/g]
  double cmf_lum; // cmf luminosity light curve [erg]
  int pellet_decays; // Number of pellets that decay in this time step.
};


/// PHIXSLIST

struct fullphixslist
{
  double nu_edge;
  int element;
  int ion;
  int level;
  int phixstargetindex;
  int index_in_groundphixslist;
};

struct groundphixslist
{
  double nu_edge;
  int element;
  int ion;
  int level;
  int phixstargetindex;
};


struct phixslist
{
  #if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
  double *groundcont_gamma_contr;
  #endif
  double *kappa_bf_sum;
#if (DETAILED_BF_ESTIMATORS_ON)
  double *gamma_contr;
#endif
};




/// GRID
///============================================================================
struct compositionlist_entry
{
  float abundance;         /// Abundance of the element (by mass!).
  float *groundlevelpop;   /// Pointer to an array of floats which contains the groundlevel populations
                           /// of all included ionisation stages for the element.
  float *partfunct;        /// Pointer to an array of floats which contains the partition functions
                           /// of all included ionisation stages for the element.
  //float *ltepartfunct;     /// Pointer to an array of floats which contains the LTE partition functions
  //                         /// of all included ionisation stages for the element.
};

struct gridcell
{
  double pos_min[3];   // Initial co-ordinates of inner most corner of cell.
  // int xyz[3];       // Integer position of cell in grid.
  int modelgridindex;
};


#define NSYN 1 /* number of frequency points in syn calculation */

enum ray_status {
  ACTIVE = 1,
  WAITING = 2,
  FINISHED = 3,
};


#include "boundary.h"

struct syn_ray
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
};



/// ATOMIC DATA
///============================================================================

struct phixstarget_entry
{
  double *spontrecombcoeff;
  #if (!NO_LUT_PHOTOION)
    double *corrphotoioncoeff;
  #endif
  #if (!NO_LUT_BFHEATING)
  double *bfheating_coeff;
  #endif
  double *bfcooling_coeff;

  double probability;        // fraction of phixs cross section leading to this final level
  int levelindex;         // index of upper ion level after photoionisation
};


struct levellist_entry
{
  double epsilon;                            /// Excitation energy of this level relative to the neutral ground level.
  int *uptrans_lineindicies;    /// Allowed upward transitions from this level
  int *downtrans_lineindicies;  /// Allowed downward transitions from this level
  int nuptrans;
  int ndowntrans;
  double phixs_threshold;                    /// Energy of first point in the photion_xs table
  struct phixstarget_entry *phixstargets;  /// pointer to table of target states and probabilities
  float *photoion_xs;               /// Pointer to a lookup-table providing photoionisation cross-sections for this level.
  int nphixstargets;                         /// length of phixstargets array:
  float stat_weight;                           /// Statistical weight of this level.

  int cont_index;                            /// Index of the continuum associated to this level. Negative number.
#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
  int closestgroundlevelcont;
#endif

  int uniquelevelindex;
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

};

struct ionlist_entry
{
  struct levellist_entry *levels;                   /// Carries information for each level: 0,1,...,nlevels-1
  int ionstage;                              /// Which ionisation stage: XI=0, XII=1, XIII=2, ...
  int nlevels;                               /// Number of levels for this ionisation stage
  int nlevels_nlte;                          /// number of nlte levels for this ion
  int first_nlte;                            /// index into nlte_pops array of a grid cell
  int ionisinglevels;                        /// Number of levels which have a bf-continuum
  int maxrecombininglevel;                   /// level index of the highest level with a non-zero recombination rate
  int nlevels_groundterm;
  int coolingoffset;
  int ncoolingterms;
  int uniqueionindex;
  float *Alpha_sp;
  double ionpot;                             /// Ionisation threshold to the next ionstage
  //int nbfcontinua;
  //ionsphixslist *phixslist;
};

struct elementlist_entry
{
  ionlist_entry *ions;              /// Carries information for each ion: 0,1,...,nions-1
  int nions;                                 /// Number of ions for the current element
  int anumber;                               /// Atomic number
//  int uppermost_ion;                       /// Highest ionisation stage which has a decent population for a given cell
                                             /// Be aware that this must not be used outside of the update_grid routine
                                             /// and their daughters. Neither it will work with OpenMP threads.
  float abundance;                           ///
  float initstablemeannucmass;                   /// Atomic mass number in multiple of MH
};

struct linelist_entry
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
};

struct bflist_t
{
  int elementindex;
  int ionindex;
  int levelindex;
  int phixstargetindex;
};

struct nne_solution_paras
{
  int cellnumber;
};

struct gslintegration_paras
{
  double nu_edge;
  double T;
  float *photoion_xs;
};

struct rpkt_cont_opacity
{
  double nu; // frequency at which opacity was calculated
  double total;
  double es;
  double ff;
  double bf;
  double ffheating;
  //double bfheating;
  int modelgridindex;
  bool recalculate_required; // e.g. when cell or timestep has changed
};

/// Cell history
///============================================================================
// typedef struct coolinglist_contributions
// {
//   double contribution;
// } coolinglist_contributions;

struct chphixstargets
{
  double corrphotoioncoeff;
#if (SEPARATE_STIMRECOMB)
  double stimrecombcoeff;
#endif
};

#include "macroatom.h"

struct chlevels
{
  double processrates[MA_ACTION_COUNT];
  struct chphixstargets *chphixstargets;
  double bfheatingcoeff;
  double population;
  double *individ_rad_deexc;
  double *individ_internal_down_same;
  double *individ_internal_up_same;
};

struct chions
{
  struct chlevels *chlevels;              /// Pointer to the ions levellist.
};

struct chelements
{
  struct chions *chions;                  /// Pointer to the elements ionlist.
};

struct cellhistory
{
  double *cooling_contrib;    /// Cooling contributions by the different processes.
  struct chelements *chelements;            /// Pointer to a nested list which helds compositional
                                            /// information for all the elements=0,1,...,nelements-1
  int cellnumber;                           /// Identifies the cell the data is valid for.
  int bfheating_mgi;
};


#endif //TYPES_H
