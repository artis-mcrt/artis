```
// Number of energy packets per process (MPI rank). OpenMP threads share these packets
constexpr int MPKTS;

constexpr int GRID_TYPE = {GRID_UNIFORM, GRID_SPHERICAL1D}
constexpr int CUBOID_NCOORDGRID_X;
constexpr int CUBOID_NCOORDGRID_Y;
constexpr int CUBOID_NCOORDGRID_Z;

// non-LTE population solver
constexpr bool NLTE_POPS_ON;

// solve the NLTE population matrix equation simultaneously for levels of all ions of an element
constexpr bool NLTE_POPS_ALL_IONS_SIMULTANEOUS;

// maximum number of NLTE/Te/Spencer-Fano iterations
constexpr int NLTEITER;

// this macro function determines which levels of which ions will be treated in full NLTE (if NLTE_POPS_ON is true)
// for now, all NLTE levels should be contiguous and include the ground state
// (i.e. level indices < X should return true for some X)
constexpr bool LEVEL_IS_NLTE(int element_z, int ionstage, int level) { return false; }

// atomic data and LTE
#define LTEPOP_EXCITATIONTEMPERATURE grid::get_Te(modelgridindex)

// Only include a single level for the highest ion stage
constexpr bool single_level_top_ion;

// if false, read from file or autodetect
constexpr bool single_ground_level;

// option to enforce connecting the lower n levels to all other levels with collisions
// disable by returning zero
constexpr int NLEVELS_REQUIRETRANSITIONS(int Z, int ionstage) {
  return ((Z == 26 || Z == 28) && ionstage >= 1) ? 80 : 0;
}

// if uniform pellet energies are not used, a uniform decay time distribution is used with scaled packet energies
constexpr bool UNIFORM_PELLET_ENERGIES;

constexpr bool DIRECT_COL_HEAT;
constexpr bool NO_INITIAL_PACKETS;
constexpr bool RECORD_LINESTAT;

// allows non-zero energy density at time t_model using q column in model.txt
// NO_INITIAL_PACKETS must be disabled to make use of this
constexpr bool USE_MODEL_INITIAL_ENERGY;

/// Rate coefficients
constexpr int TABLESIZE;
constexpr double MINTEMP;
constexpr double MAXTEMP;

// temperature for which total ion recombination rate are calibrated to input data (recombrates.txt)
constexpr double RECOMBCALIBRATION_T_ELEC;

// Polarisation for real packets
constexpr bool DIPOLE;
constexpr bool POL_ON;

// Polarisation for virtual packets
// #define VPKT_ON


constexpr bool TRACK_ION_STATS;
constexpr bool TRACK_ION_MASTATS;

constexpr double MINPOP;

constexpr double NU_MIN_R;  /// lower frequency boundary for UVOIR spectra and BB sampling
constexpr double NU_MAX_R;  /// upper frequency boundary for UVOIR spectra and BB sampling

// ** Start of radiation field model options **

// if using this, avoid look up tables and switch on the direct integration options below
// (since LUTs created with Planck function J_nu)
constexpr bool MULTIBIN_RADFIELD_MODEL_ON;

constexpr int RADFIELDBINCOUNT;

constexpr int FIRST_NLTE_RADFIELD_TIMESTEP;

constexpr double nu_lower_first_initial;  // CLIGHT / ([lambda Angstroms]e-8)
constexpr double nu_upper_last_initial;    // not including the very top super bin
constexpr double nu_upper_superbin;          // very top end super bin

constexpr double T_R_min;
constexpr double T_R_max;

// store Jb_lu estimators for particular lines chosen in radfield::init()
constexpr bool DETAILED_LINE_ESTIMATORS_ON;

// store detailed bound-free rate estimators
constexpr bool DETAILED_BF_ESTIMATORS_ON;

// if DETAILED_BF_ESTIMATORS_ON, then use BF estimators at the following timestep and later
constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP;

// dynamically calculate photoionization rates for the current radiation field
// instead of interpolating values from a lookup table for a blackbody radiation field
constexpr bool NO_LUT_PHOTOION;

// as above for bound-free heating
constexpr bool NO_LUT_BFHEATING;

// if SEPARATE_STIMRECOMB is false, then stimulated recombination is treated as negative photoionisation
#define SEPARATE_STIMRECOMB false

// ** End of radiation field model options **

// ** Start of non-thermal solution options **

// non-thermal ionisation
constexpr bool NT_ON;

// use the detailed Spencer-Fano solver instead of the work function approximation (only works if NT_ON)
constexpr bool NT_SOLVE_SPENCERFANO;

// number of energy points in the Spencer-Fano solution vector
constexpr int SFPTS;

// eV
constexpr double SF_EMAX;

// eV
constexpr double SF_EMIN;

// trigger a Spencer-Fano solution at least once every n timesteps
// 0 can only use solutions from previous NLTE iterations on the current timestep
// <=-1 will always solve the SF equation for every iteration of every timestep
constexpr int SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS;

// a change in the electron fraction (e.g. 0.5 is a 50% change) since the previous solution will also trigger a solution
constexpr double NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS;

// just consider excitation from the first N levels and to the first M upper levels,
// because these transitions really slow down the solver
constexpr int NTEXCITATION_MAXNLEVELS_LOWER;    // set to zero for none
constexpr int NTEXCITATION_MAXNLEVELS_UPPER;  // maximum number of upper levels included

// limit the number of stored non-thermal excitation transition rates to reduce memory cost.
// if this is higher than SFPTS, then you might as well just store
// the full NT degradation spectrum and calculate the rates as needed (although CPU costs)
constexpr int MAX_NT_EXCITATIONS_STORED;

// set to true to keep a list of non-thermal excitation rates for use
// in the NLTE pop solver, macroatom, and NTLEPTON packets.
// Even with this off, excitations will be included in the solution
// and their combined deposition fraction is calculated
constexpr bool NT_EXCITATION_ON = false;

// increase the excitation and ionization lists by this blocksize when reallocating
constexpr int NT_BLOCKSIZEEXCITATION;

// calculate eff_ionpot and ionisation rates by always dividing by the valence shell potential for the ion
// instead of the specific shell potentials
constexpr bool NT_USE_VALENCE_IONPOTENTIAL;

// allow ions to lose more than one electron per impact ionisation using Auger effect probabilities
// associate with electron shells
// if this is greater than zero, make sure NT_USE_VALENCE_IONPOTENTIAL is false!
constexpr int NT_MAX_AUGER_ELECTRONS;

// add the Auger electron term to the Spencer-Fano equation
constexpr bool SF_AUGER_CONTRIBUTION_ON;

// set true to divide up the mean Auger energy by the number of electrons that come out
constexpr bool SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN;

// ** End of non-thermal solution options **

constexpr double TEMPERATURE_SOLVER_ACCURACY;

constexpr double CONTINUUM_NU_INTEGRAL_ACCURACY;

constexpr double RATECOEFF_INTEGRAL_ACCURACY;

// when calculating ion ionisation rate coefficient (for estimator files), contribute the lowest n levels that
// make up at least IONGAMMA_POPFRAC_LEVELS_INCLUDED fraction of the ion population
constexpr double IONGAMMA_POPFRAC_LEVELS_INCLUDED;

constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT;

// when converting mass fraction to a number density, use a mean atomic mass
// calcuated from the nuclear composition (plus stable component),
// rather than just from the compositiondata.txt values
constexpr bool USE_CALCULATED_MEANATOMICWEIGHT;

constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC;

constexpr bool INSTANT_PARTICLE_DEPOSITION;

// Options for different types of timestep set-ups, only one of these can be true at one time. The hybrid timestep
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

constexpr enum timestepsizemethods TIMESTEP_SIZE_METHOD;

constexpr double FIXED_TIMESTEP_WIDTH;

constexpr double TIMESTEP_TRANSITION_TIME;

constexpr bool USE_GSL_RANDOM;

// once a new gridsave and packets*.tmp have been written, don't delete the previous set
constexpr bool KEEP_ALL_RESTART_FILES;

// multiply bound-free cooling coefficient by upper level population instead of the upper ion target level population
constexpr bool BFCOOLING_USELEVELPOPNOTIONPOP;

// classic artis used a comparison function with side effects and the resulting line list
// was not correctly sorted by frequency
constexpr bool CLASSIC_LINELIST_SORT;
```