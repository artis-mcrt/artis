#ifndef ARTISOPTIONS_H
#define ARTISOPTIONS_H

#include <cstdlib>
#include "constants.h"

// Number of energy packets per process (MPI rank). OpenMP threads share these packets
#define MPKTS 1000000

// Max number of propagation grid cells
//125000 //1000000 //1000000//262144 //2100000 //125000 //1000000
#define MGRID  125000

// Max number of input model grid cells
//125000 //12800 //12800 //125 //3200 //200 //200 //200 //8192 //125 //125000 //200 //125000 //8200 //200 //8200 //200 //125000
#define MMODELGRID 177

#define GRID_TYPE GRID_UNIFORM
#define CUBOID_NCOORDGRID_X 50
#define CUBOID_NCOORDGRID_Y 50
#define CUBOID_NCOORDGRID_Z 50
// #define GRID_TYPE GRID_SPHERICAL1D

// Max number of elements
#define MELEMENTS 17

// Max number of ion stages for any element
#define MIONS 5

// non-LTE population solver
static const bool NLTE_POPS_ON = true;

// solve the NLTE population matrix equation simultaneously for levels of all ions of an element
static const bool NLTE_POPS_ALL_IONS_SIMULTANEOUS = true;

// maximum number of NLTE/Te/Spencer-Fano iterations
static const int NLTEITER = 30;

// this macro function determines which levels of which ions will be treated in full NLTE (if NLTE_POPS_ON is true)
// for now, all NLTE levels should be contiguous and include the ground state
// (i.e. level indices < X should return true for some X)
#define LEVEL_IS_NLTE(element, ion, level) \
if (get_element(element) == 26 && get_ionstage(element, ion) == 2) \
return (level <= 197); \
else \
return (level <= 80);

// if uniform pellet energies are not used, a uniform decay time distribution is used with scaled packet energies
#define UNIFORM_PELLET_ENERGIES true

#define DIRECT_COL_HEAT
#define NO_INITIAL_PACKETS
#define RECORD_LINESTAT

/// Rate coefficients
#define TABLESIZE 100 //200 //100
#define MINTEMP 1000.
#define MAXTEMP 30000. //1000000.

// temperature for which total ion recombination rate are calibrated to input data (recombrates.txt)
#define RECOMBCALIBRATION_T_ELEC 6000.

// Polarisation for real packets
// #define DIPOLE
// #define POL_ON

// Polarisation for virtual packets
// #define VPKT_ON

// GSL integration workspace size
static const size_t GSLWSIZE = 16384;


#define TRACK_ION_STATS false
#define TRACK_ION_MASTATS false

#define MLINES 500000    // Increase linelist by this blocksize

#define MINDENSITY 1e-40         /// Minimum cell density. Below cells are treated as empty.
#define MINPOP 1e-40

#define NU_MIN_R 1e13   /// lower frequency boundary for UVOIR spectra and BB sampling
#define NU_MAX_R 5e15   /// upper frequency boundary for UVOIR spectra and BB sampling


// ****
// Start of radiation field model options
//

// if using this, avoid look up tables and switch on the direct integration options below
// (since LUTs created with Planck function J_nu)
static const bool MULTIBIN_RADFIELD_MODEL_ON = true;

#define RADFIELDBINCOUNT 256

static const int FIRST_NLTE_RADFIELD_TIMESTEP = 12;

static const double nu_lower_first_initial = (CLIGHT / (40000e-8)); // CLIGHT / ([lambda Angstroms]e-8)
static const double nu_upper_last_initial = (CLIGHT /  (1085e-8));  // not including the very top super bin
static const double nu_upper_superbin = (CLIGHT /  (10e-8)); // very top end super bin

static const double T_R_min = 500;
static const double T_R_max = 250000;

// store Jb_lu estimators for particular lines chosen in radfield::init()
static const bool DETAILED_LINE_ESTIMATORS_ON = false;

// store detailed bound-free rate estimators
#define DETAILED_BF_ESTIMATORS_ON true

// if DETAILED_BF_ESTIMATORS_ON, then use BF estimators at the following timestep and later
#define DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP 13

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

// if SEPARATE_STIMRECOMB is false, then stimulated recombination is treated as negative photoionisation
#define SEPARATE_STIMRECOMB false

//
// End of radiation field model options
// ****


// ****
// Start of non-thermal solution options
//

/// non-thermal ionisation
static const bool NT_ON = true;

/// use the detailed Spencer-Fano solver instead of the work function approximation
static const bool NT_SOLVE_SPENCERFANO = true;

// number of energy points in the Spencer-Fano solution vector
#define SFPTS 4096

// eV
#define SF_EMAX 16000.

// eV
#define SF_EMIN .1

// use a grid of energy points with constant spacing in log energy
#define SF_USE_LOG_E_INCREMENT false

// trigger a Spencer-Fano solution at least once every n timesteps
// 0 can only use solutions from previous NLTE iterations on the current timestep
// <=-1 will always solve the SF equation for every iteration of every timestep
static const int SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS = 0;

// a change in the electron fraction (e.g. 0.5 is a 50% change) since the previous solution will also trigger a solution
static const double NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS = 0.05;

// just consider excitation from the first N levels and to the first M upper levels,
// because these transitions really slow down the solver
static const int NTEXCITATION_MAXNLEVELS_LOWER = 5;  // set to zero for none
static const int NTEXCITATION_MAXNLEVELS_UPPER = 250; // maximum number of upper levels included

// limit the number of stored non-thermal excitation transition rates to reduce memory cost.
// if this is higher than SFPTS, then you might as well just store
// the full NT degradation spectrum and calculate the rates as needed (although CPU costs)
static const int MAX_NT_EXCITATIONS_STORED = 25000;

// set to true to keep a list of non-thermal excitation rates for use
// in the NLTE pop solver, macroatom, and NTLEPTON packets.
// Even with this off, excitations will be included in the solution
// and their combined deposition fraction is calculated
#define NT_EXCITATION_ON false

// increase the excitation and ionization lists by this blocksize when reallocating
#define NT_BLOCKSIZEEXCITATION 5192

// calculate eff_ionpot and ionisation rates by always dividing by the valence shell potential for the ion
// instead of the specific shell potentials
#define NT_USE_VALENCE_IONPOTENTIAL false

// allow ions to lose more than one electron per impact ionisation using Auger effect probabilities
// associate with electron shells
// if this is greater than zero, make sure NT_USE_VALENCE_IONPOTENTIAL is false!
#define NT_MAX_AUGER_ELECTRONS 2

// add the Auger electron term to the Spencer-Fano equation
#define SF_AUGER_CONTRIBUTION_ON true

// set true to divide up the mean Auger energy by the number of electrons that come out
#define SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN false

//
// End of non-thermal solution options
// ****


#endif //ARTISOPTIONS_H