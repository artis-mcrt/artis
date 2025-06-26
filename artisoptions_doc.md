```
// Number of energy packets per process (MPI rank). OpenMP threads share these packets
constexpr int MPKTS;

constexpr auto GRID_TYPE = {GridType::CARTESIAN3D, GridType::CYLINDRICAL2D, GridType::SPHERICAL1D}

// for GridType::CARTESIAN3D, set the dimensions. This will have no effect with a 3D model.txt since they will be set to match the input
constexpr int CUBOID_NCOORDGRID_X;
constexpr int CUBOID_NCOORDGRID_Y;
constexpr int CUBOID_NCOORDGRID_Z;

// for 2D cylindrical and 3D Cartesian, remove the corners (v > vmax) to force a spherical escape surface
constexpr bool FORCE_SPHERICAL_ESCAPE_SURFACE;

// maximum number of NLTE/Te/Spencer-Fano iterations
constexpr int NLTEITER;

// The return of this function is passed to LEVEL_IS_NLTE and specifies how many levels will be treated in full NLTE.
// It is only active when LEVEL_IS_NLTE does not return false. The additiopn of this function enables
// NLEVELS_REQUIRETRANSITIONS to have a value specified for each ion relative to the number of NLTE levels for that ion:
// It is advisable that if NLEVELS_REQUIRETRANSITIONS is active at least a subset of the levels treated in the superlevel
// are alweays collisionally connected to lower levels even when transition data doesn't exist to avoid population versions
// occuring in the superlevel
constexpr int ION_NLEVELS_NLTE(int element_z, int ionstage) {
  if (element_z < 20) {
    return 100;
  }
  return 200;
}

// this macro function determines which levels of which ions will be treated in full NLTE
// for now, all NLTE levels should be contiguous and include the ground state
// (i.e. level indices < X should return true for some X)
constexpr bool LEVEL_IS_NLTE(int element_z, int ionstage, int level) { return false; }

// Use TJ radiation density temperature for Boltzmann excitation formula instead of electron temperature Te
// This is default on for classic, and off for nebularnlte, where it affects the super-level
constexpr bool LTEPOP_EXCITATION_USE_TJ = false;

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

// INITIAL PACKETS will seed the cells on the first timestep at tmin with K-packets
// representing decay energy from t_model to tmin, and,
// if USE_MODEL_INITIAL_ENERGY is true, also the snapshot energy at t_model
constexpr bool INITIAL_PACKETS_ON;

// allows non-zero energy density at time t_model using q column in model.txt
// INITIAL_PACKETS_ON must be true to make use of this
constexpr bool USE_MODEL_INITIAL_ENERGY;

// record counts of emissions and absorptions in each line
constexpr bool RECORD_LINESTAT;

// Rate coefficients
constexpr int TABLESIZE;
constexpr double MINTEMP;
constexpr double MAXTEMP;

// temperature for which total ion recombination rate are calibrated to input data (recombrates.txt)
constexpr double RECOMBCALIBRATION_T_ELEC;

// Polarisation for real packets
constexpr bool DIPOLE;

// Only affects exspec and enables writing specpol.out, emissionpol.out, absorptionpol.out
constexpr bool POL_ON;

// Polarisation for virtual packets
constexpr bool VPKT_ON;

constexpr bool TRACK_ION_STATS;

constexpr double MINPOP;

constexpr double NU_MIN_R;  // lower frequency boundary for UVOIR spectra and BB sampling
constexpr double NU_MAX_R;  // upper frequency boundary for UVOIR spectra and BB sampling

// use nearest-neighbour instead of linear interpolation of photoionisation cross sections
// to match classic artis
constexpr bool PHIXS_CLASSIC_NO_INTERPOLATION;

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

// interpolate values from a lookup table for a blackbody radiation field
// instead of dynamically integrating photoionization rates for the exact radiation field
constexpr bool USE_LUT_PHOTOION;

// as above for bound-free heating
constexpr bool USE_LUT_BFHEATING;

// Previously the NLTE solver only checked if level populations were negative and replaced these populations
// with the LTE population. However this can cause numerical problems (e.g. when the ground populations is very
// small and and the negative population is replaced with a significantly larger population ratios taken e.g.
// when calculating partition functions can over the float limit.) This option provides additional checks
// on the populations calculated by the NLTE solver (ground population > MINPOP, checks for population inversions)
// and now returns a solver fail for certain cases.
constexpr bool STRICT_POPULATION_CHECKING = false;

// Previously when the NLTE solution failed the populations for the entire element were set to LTE values.
// However it is often the uppermost/lowermost ions which have problems in their solution due to the
// dynamic range of the simulation (e.g. high ion stages which have small populations in cells which
// are in a low ionisation state and vice versa) while a good NLTE solution can be obtained for the rest
// of the ions in the element. This option provides functionality to strip the top/bottom ions progressively
// from elements (provided they have small populations) when the NLTE fails before retrying the solution
// to determine if a successful solution can be obtained with this reduced range of ions.
constexpr bool NLTE_LIMIT_ION_STAGES_AFTER_FAILURE = false;

// Controls by what factor the populations of a level have to be inverted relative to the ground to result
// in a NLTE solver fail being returned
constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL = 1000.;

// Controls by what factor the populations of a level have to be inverted relative to print out a warning
// that the level is inverted
constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING = 10.;

// Controls the highest ratio the population a level can have relative to the total population of the element
// for the ion to still be removed from the NLTE solution when using the NLTE_LIMIT_ION_STAGES_AFTER_FAILURE
// functionality. The ratio of the nlte ground populations, each nlte excited population and the superlevel
// population are all individually checked.
constexpr double NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION = 1e-9;

// If set_element_pops_lte is called when NLTE solver fails the uppermost_ion used in set_groundlevelpops
// is the one in set in memory based on the NLTE phi factors. In rare cases this can result in partition
// function overflows. This option recalls find_uppermost_ion so the uppermost_ion used is based on the
// correct LTE phi factors instead.
constexpr bool NLTEPOP_FAILURE_USE_FIND_UPPERMOST_ION_FOR_LTE_RESET = false;

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
// 0 can only re-use solutions from previous NLTE iterations of the current timestep
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

// load shells.txt containing shell occupancy data instead of simple algorithmic guesses
constexpr bool NT_WORKFUNCTION_USE_SHELL_OCCUPANCY_FILE = false;

// ** End of non-thermal solution options **

constexpr double TEMPERATURE_SOLVER_ACCURACY;

constexpr double CONTINUUM_NU_INTEGRAL_ACCURACY;

constexpr double RATECOEFF_INTEGRAL_ACCURACY;

// when calculating ion ionisation rate coefficient (for estimator files), contribute the lowest n levels that
// make up at least IONGAMMA_POPFRAC_LEVELS_INCLUDED fraction of the ion population
constexpr double IONGAMMA_POPFRAC_LEVELS_INCLUDED;

constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT;

// when converting mass fraction to a number density, use a mean atomic mass
// calculated from the nuclear composition (plus stable component),
// rather than just from the compositiondata.txt values
constexpr bool USE_CALCULATED_MEANATOMICWEIGHT;

constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC;

constexpr bool INSTANT_PARTICLE_DEPOSITION;

// Options for different types of timestep set-ups, only one of these can be true at one time. The hybrid timestep
// schemes that switch between log and fixed require a transition time from one scheme to the other as well as the
// fixed timestep width to be set. These values need to be consistent with the number of timesteps i.e. don't give
// values that would give the same number or more more fixed timesteps than the total number of timesteps in the
// simulation. The times are set in days.


constexpr enum class timestepsizemethods TIMESTEP_SIZE_METHOD;

constexpr double FIXED_TIMESTEP_WIDTH;

constexpr double TIMESTEP_TRANSITION_TIME;

// once a new gridsave and packets*.tmp have been written, don't delete the previous set
constexpr bool KEEP_ALL_RESTART_FILES;

// multiply bound-free cooling coefficient by upper level population instead of the upper ion target level population
constexpr bool BFCOOLING_USELEVELPOPNOTIONPOP;

// set true to calculate and use expansion opacities instead of line-by-line
constexpr bool EXPANSIONOPACITIES_ON;

// thermalisation probability (1 - P is probabiltiy of scattering). EXPANSIONOPACITIES_ON must be true for this to work.
// set this to < 0 to use the macroatom
constexpr float RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY;

// Use XCOM data for gamma photoionisation instead of Si+Fe Equation 2 of Ambwani & Sutherland (1988), Veigele (1973)
constexpr bool USE_XCOM_GAMMAPHOTOION;

```
