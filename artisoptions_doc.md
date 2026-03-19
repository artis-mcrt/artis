```c++
// Number of energy packets per process (MPI rank). OpenMP threads share these packets
constexpr int MPKTS;

// set to GridType:: CARTESIAN3D, CYLINDRICAL2D, or SPHERICAL1D
constexpr GridType GRID_TYPE;

// for GridType::CARTESIAN3D, set the dimensions. This will have no effect with a 3D model.txt since they will be set to match the input
constexpr int CUBOID_NCOORDGRID_X;
constexpr int CUBOID_NCOORDGRID_Y;
constexpr int CUBOID_NCOORDGRID_Z;

// for 2D cylindrical and 3D Cartesian, remove the corners (v > vmax) to force a spherical escape surface
constexpr bool FORCE_SPHERICAL_ESCAPE_SURFACE;

// maximum number of NLTE/Te/Spencer-Fano iterations
constexpr int NLTEITER;

// Specify how many levels will be treated in full NLTE, not including the ground state or the superlevel.
constexpr int ION_NLEVELS_EXCITED_NLTE(int element_z, int ionstage);

// Use TJ radiation density temperature for Boltzmann excitation formula instead of electron temperature Te
// This is default on for classic (with Boltzmann factor level pops), and off for nebularnlte, where it affects the superlevel sublevel populations
constexpr bool LTEPOP_EXCITATION_USE_TJ;

// force Saha ionisation balance for a given element (contraint applied to NLTE population solver and classic phi function)
constexpr bool FORCE_SAHA_ION_BALANCE(int element_z);

// Only include a single level for the highest ion stage
constexpr bool SINGLE_LEVEL_TOP_ION;

// if false, read from file or autodetect
// this only affects the recombrates scaling, since rates are given per ground multiplet population.
constexpr bool SINGLE_GROUND_LEVEL;

// Add any missing collisional transitions between the lower n levels and all other levels (or disable by returning zero)
// This can prevent fully disconnected levels, whose NLTE populations cannot be determined
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

// write virtual packet per-direction contributions to vpkt_contrib.out files
constexpr bool VPKT_WRITE_CONTRIBS;

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

constexpr double RADFIELDBINS_NU_MIN;  // CLIGHT / ([lambda Angstroms]e-8)
constexpr double RADFIELDBINS_NU_MAX;    // not including the very top super bin
constexpr double RADFIELDBINS_T_E_SUPERBIN_NU_MAX;          // very top end super bin

// store Jb_lu estimators for particular lines chosen in radfield::init()
constexpr bool DETAILED_LINE_ESTIMATORS_ON;

// store detailed bound-free rate estimators
constexpr bool DETAILED_BF_ESTIMATORS_ON;

// select which bf-continua are tracked in the detailed estimators (only used when DETAILED_BF_ESTIMATORS_ON is true)
constexpr bool LEVEL_HAS_BFEST(int element_z, int ionstage, int level);

// if DETAILED_BF_ESTIMATORS_ON, then use BF estimators at the following timestep and later
constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP;

// interpolate values from a lookup table for a blackbody radiation field
// instead of dynamically integrating photoionisation rates for the exact radiation field
constexpr bool USE_LUT_PHOTOION;

// enable per-ion bound-free heating estimators and associated renormalisation
constexpr bool USE_ION_BFHEATING_ESTIMATORS;

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
constexpr bool NLTE_LIMIT_ION_STAGES_AFTER_FAILURE;

// Controls by what factor the populations of a level have to be inverted relative to the ground to result
// in a NLTE solver fail being returned
constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL;

// Controls by what factor the populations of a level have to be inverted relative to print out a warning
// that the level is inverted
constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING;

// Controls the highest ratio the population a level can have relative to the total population of the element
// for the ion to still be removed from the NLTE solution when using the NLTE_LIMIT_ION_STAGES_AFTER_FAILURE
// functionality. The ratio of the nlte ground populations, each nlte excited population and the superlevel
// population are all individually checked.
constexpr double NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION;

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

constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT;

// when converting mass fraction to a number density, use a mean atomic mass
// calculated from the nuclear composition (plus stable component),
// rather than just from the compositiondata.txt values
constexpr bool USE_CALCULATED_MEANATOMICWEIGHT;

constexpr bool WRITE_EMISSIONABSORPTION_SPEC_AT_END;

// thermalisation scheme for non-thermal particles (positrons, electrons, alphas). ThermalisationScheme::INSTANT
// instantly deposits all particle energy. ThermalisationScheme::DETAILED uses time-dependent Monte Carlo transport.
// ThermalisationScheme::BARNES, and WOLLAEGER use analytic thermalisation efficiency functions.
constexpr ThermalisationScheme PARTICLE_THERMALISATION_SCHEME;

// thermalisation scheme for gamma-ray photons. ThermalisationScheme::DETAILED uses full gamma-ray transport.
constexpr ThermalisationScheme GAMMA_THERMALISATION_SCHEME;

// Options for different types of timestep set-ups, only one of these can be true at one time. The hybrid timestep
// schemes that switch between log and fixed require a transition time from one scheme to the other as well as the
// fixed timestep width to be set. These values need to be consistent with the number of timesteps i.e. don't give
// values that would give the same number or more more fixed timesteps than the total number of timesteps in the
// simulation. The times are set in days.

constexpr TimeStepSizeMethod TIMESTEP_SIZE_METHOD;

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
constexpr std::optional<float> RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY;

// Use XCOM data for gamma photoionisation instead of Si+Fe Equation 2 of Ambwani & Sutherland (1988), Veigele (1973)
constexpr bool USE_XCOM_GAMMAPHOTOION;

// use fissiondecays.txt and fissionproducts_GEF_100keV.txt to handle spontaneous fission decays
constexpr bool DECAY_SPONTFISSION_ON = false;

// TODO: explain where clumping factors are actually used
constexpr bool USE_MICROCLUMPING;

// TODO: figure out the arguments to this function
constexpr float clumping_factor(double tmid, double rad_vel) { return 1.; }
```
