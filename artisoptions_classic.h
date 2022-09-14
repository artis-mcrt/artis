#ifndef ARTISOPTIONS_H
#define ARTISOPTIONS_H

#include <cstdlib>

#include "constants.h"

// Number of energy packets per process (MPI rank). OpenMP threads share these packets
#define MPKTS 100000

#define GRID_TYPE GRID_UNIFORM
#define CUBOID_NCOORDGRID_X 100
#define CUBOID_NCOORDGRID_Y 100
#define CUBOID_NCOORDGRID_Z 100
// #define GRID_TYPE GRID_SPHERICAL1D

// non-LTE population solver
constexpr bool NLTE_POPS_ON = false;

// solve the NLTE population matrix equation simultaneously for levels of all ions of an element
constexpr bool NLTE_POPS_ALL_IONS_SIMULTANEOUS = false;

// maximum number of NLTE/Te/Spencer-Fano iterations
constexpr int NLTEITER = 30;

// this macro function determines which levels of which ions will be treated in full NLTE (if NLTE_POPS_ON is true)
// for now, all NLTE levels should be contiguous and include the ground state
// (i.e. level indices < X should return true for some X)
#define LEVEL_IS_NLTE(element, ion, level) return false;

// atomic data and LTE
#define LTEPOP_EXCITATIONTEMPERATURE grid::get_TJ(modelgridindex)

// Only include a single level for the highest ion stage
constexpr bool single_level_top_ion = true;

// if false, read from file or autodetect
constexpr bool single_ground_level = true;

// option to enforce connecting the lower n levels to all other levels with collisions
// disable by returning zero
#define NLEVELS_REQUIRETRANSITIONS(Z, ionstage) 0

// if uniform pellet energies are not used, a uniform decay time distribution is used with scaled packet energies
#define UNIFORM_PELLET_ENERGIES true

// #define DIRECT_COL_HEAT
// #define NO_INITIAL_PACKETS
#define RECORD_LINESTAT

/// Rate coefficients
#define TABLESIZE 100  // 200 //100
#define MINTEMP 3500.
#define MAXTEMP 140000.  // 1000000.

// temperature for which total ion recombination rate are calibrated to input data (recombrates.txt)
#define RECOMBCALIBRATION_T_ELEC 6000.

// Polarisation for real packets
#define DIPOLE
#define POL_ON

// Polarisation for virtual packets
// #define VPKT_ON

// GSL integration workspace size
constexpr size_t GSLWSIZE = 16384;

#define TRACK_ION_STATS false
#define TRACK_ION_MASTATS false

// Increase linelist by this blocksize
#define MLINES 500000

// Minimum cell density. Below cells are treated as empty.
#define MINDENSITY 1e-40
#define MINPOP 1e-30

// lower frequency boundary for UVOIR spectra and BB sampling
#define NU_MIN_R 1e14
// upper frequency boundary for UVOIR spectra and BB sampling
#define NU_MAX_R 5e15

// ** Start of radiation field model options **

// if using this, avoid look up tables and switch on the direct integration options below
// (since LUTs created with Planck function J_nu)
constexpr bool MULTIBIN_RADFIELD_MODEL_ON = false;

// number of bins (including the T=T_e superbin)
#define RADFIELDBINCOUNT 256

// from this timestep number (0-indexed), radfield switches from T_R,W dilute blackbody to binned radfield
constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 12;

// nu_lower of lowest-frequency bin = CLIGHT / ([lambda Angstroms]e-8)
constexpr double nu_lower_first_initial = (CLIGHT / (40000e-8));
// nu_upper of highest-frequency bin just below the very top super bin
constexpr double nu_upper_last_initial = (CLIGHT / (1085e-8));
// nu_upper of top super bin with fixed T=T_e (nu_lower will be nu_upper_last_initial)
constexpr double nu_upper_superbin = (CLIGHT / (10e-8));

// range of 'temperatures' for the scaled Planck function fit to each bin
constexpr double T_R_min = 500;
constexpr double T_R_max = 250000;

// store Jb_lu estimators for particular lines (which are chosen in radfield::init())
constexpr bool DETAILED_LINE_ESTIMATORS_ON = false;

// track detailed bound-free rate estimators instead of doing radiation field model integrals or lookup table
#define DETAILED_BF_ESTIMATORS_ON false

// if DETAILED_BF_ESTIMATORS_ON, then use BF estimators at the following timestep (0-indexed) and later
#define DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP 13

// extremely slow and memory consuming - for debugging only
// not safe for MPI or OpenMP - single process and single thread only!
// this will output a list of contributions to each bound-free rate estimator
// with each packet emission type ranked by their contribution to the rate
#define DETAILED_BF_ESTIMATORS_BYTYPE false

// dynamically calculate photoionization rates for the current radiation field
// instead of interpolating values from a lookup table for a blackbody radiation field
#define NO_LUT_PHOTOION false

// as above for bound-free heating
#define NO_LUT_BFHEATING false

// if SEPARATE_STIMRECOMB is false, then stimulated recombination is treated as negative photoionisation
#define SEPARATE_STIMRECOMB false

// ** End of radiation field model options **

// ** Start of non-thermal solution options **

/// non-thermal ionisation
constexpr bool NT_ON = false;

/// use the detailed Spencer-Fano solver instead of the work function approximation (only works if NT_ON)
constexpr bool NT_SOLVE_SPENCERFANO = false;

// number of energy points in the Spencer-Fano solution vector
#define SFPTS 4096

// eV
#define SF_EMAX 16000.

// eV
#define SF_EMIN 0.1

// use a grid of energy points with constant spacing in log energy
#define SF_USE_LOG_E_INCREMENT false

// trigger a Spencer-Fano solution at least once every n timesteps
// 0 can only use solutions from previous NLTE iterations on the current timestep
// <=-1 will always solve the SF equation for every iteration of every timestep
constexpr int SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS = 0;

// a change in the electron fraction (e.g. 0.5 is a 50% change) since the previous solution will also trigger a solution
constexpr double NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS = 0.05;

// just consider excitation from the first N levels and to the first M upper levels,
// because these transitions really slow down the solver
constexpr int NTEXCITATION_MAXNLEVELS_LOWER = 5;    // set to zero for none
constexpr int NTEXCITATION_MAXNLEVELS_UPPER = 250;  // maximum number of upper levels included

// limit the number of stored non-thermal excitation transition rates to reduce memory cost.
// if this is higher than SFPTS, then you might as well just store
// the full NT degradation spectrum and calculate the rates as needed (although CPU costs)
constexpr int MAX_NT_EXCITATIONS_STORED = 25000;

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

// ** End of non-thermal solution options **

#define TEMPERATURE_SOLVER_ACCURACY 1e-2

#define CONTINUUM_NU_INTEGRAL_ACCURACY 1e-2

#define RATECOEFF_INTEGRAL_ACCURACY 1e-2

// when calculating ion ionisation rate coefficient (for estimator files only, not the simulation), contribute the
// lowest n levels that make up at least IONGAMMA_POPFRAC_LEVELS_INCLUDED fraction of the ion population
#define IONGAMMA_POPFRAC_LEVELS_INCLUDED 1.

// incomplete work in progress
constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT = false;

// when converting mass fraction to a number density, use a mean atomic mass
// calcuated from the nuclear composition (plus stable component),
// rather than just from the compositiondata.txt values
constexpr bool USE_CALCULATED_MEANATOMICWEIGHT = false;

// output emission.out, absorption.out during the sn3d simulation (without running exspec)
// this can occupy a lot of memory with large atomic data sets
constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = false;

// setting to false is highly experimental
constexpr bool INSTANT_PARTICLE_DEPOSITION = true;

// Options for different types of timestep set_ups, only one of these can be true at one time. The hybrid timestep
// schemes that switch between log and fixed require a transition time from one scheme to the other as well as the
// fixed timestep width to be set. These values need to be consistent with the number of timesteps i.e. don't give
// values that would give the same number or more more fixed timesteps than the total number of timesteps in the
// simulation. The times are set in days.

enum timestepsizemethods {
  TIMESTEP_SIZES_LOGARITHMIC = 0,
  TIMESTEP_SIZES_CONSTANT = 1,
  TIMESTEP_SIZES_LOGARITHMIC_THEN_CONSTANT = 2,
  TIMESTEP_SIZES_CONSTANT_THEN_LOGARITHMIC = 3,
};

constexpr enum timestepsizemethods TIMESTEP_SIZE_METHOD = TIMESTEP_SIZES_LOGARITHMIC;

constexpr double FIXED_TIMESTEP_WIDTH = 0.1;

constexpr double TIMESTEP_TRANSITION_TIME = 5;

#endif  // ARTISOPTIONS_H