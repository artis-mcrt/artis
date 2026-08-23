```c++
// Number of energy packets per process (MPI rank). OpenMP threads share these packets
constexpr int MPKTS;

// override to GridType::CARTESIAN3D, GridType::CYLINDRICAL2D, GridType::SPHERICAL1D, or leave with no value to autodetect from model.txt
constexpr std::optional<GridType> GRID_TYPE_OVERRIDE;

// for GridType::CARTESIAN3D, set the grid size from 1D and 2D input models. For a 3D model.txt, these will be ignored and set to match the input grid.
constexpr int CUBOID_NCOORDGRID_X;
constexpr int CUBOID_NCOORDGRID_Y;
constexpr int CUBOID_NCOORDGRID_Z;

// for 2D cylindrical and 3D Cartesian, remove the corners (v > vmax) to force a spherical escape surface
constexpr bool FORCE_SPHERICAL_ESCAPE_SURFACE;

// maximum number of NLTE/Te/Spencer-Fano iterations
constexpr int NLTEITER;

// Anderson acceleration of the NLTE/Te/Spencer-Fano iteration. The iteration maps the electron temperature
// and the electron density of one pass to the values of the next pass. With the acceleration on, the
// accelerator combines the last iterates so that the residual of the combination is minimal. The next
// pass then starts from that combination. The state has two components, so the accelerator keeps two
// iterates, the largest useful depth. False keeps the plain successive substitution.
//
// The accelerator forgets its history when an element falls back to LTE or changes its solved ion range,
// because the map is then discontinuous. A new non-thermal solution also clears the history. The loop
// rejects four kinds of step:
// - a step outside [MINTEMP, MAXTEMP];
// - a step with an electron density below MINPOP;
// - a step with an electron density above the total electron density of the cell;
// - a step away from the map output that is larger than twice the residual.
// The convergence test comes before the injection, so a converged cell keeps the output of a plain pass. A
// cell that reaches NLTEITER also ends with a plain pass.
constexpr bool NLTE_TE_NNE_USE_ANDERSON_ACCEL;

// Relative tolerance of the NLTE/Te/Spencer-Fano iteration for nne and T_e. Without the acceleration the
// loop stops when the change between two passes is below the tolerance. With acceleration the loop tests an
// estimate of the error that remains instead. The estimate is the change times rho / (1 - rho). Here rho is
// the ratio of the last two changes, limited to the range 0.5 to 0.95. nne and T_e each get their own ratio,
// and the larger estimate counts. The estimate is the error of a linear contraction with that ratio, so the
// tolerance then means an error and not a change. The estimate is never below the change. With the charge
// transfer reactions or the acceleration, the same value is the tolerance of the ion population test.
constexpr double NLTE_OUTER_RELTOL;

// Specify how many levels will be treated in full NLTE, not including the ground state or the superlevel.
constexpr int ION_NLEVELS_EXCITED_NLTE(int element_z, int ionstage);

// Use TJ radiation density temperature for Boltzmann excitation formula instead of electron temperature Te
// This is default on for classic (with Boltzmann factor level pops), and off for nebularnlte, where it affects the superlevel sublevel populations
constexpr bool LTEPOP_EXCITATION_USE_TJ;

// force Saha ionisation balance for a given element (constraint applied to NLTE population solver and classic phi function)
constexpr bool FORCE_SAHA_ION_BALANCE(int element_z);

// Only include a single level for the highest ion stage
constexpr bool SINGLE_LEVEL_TOP_ION;

// Add any missing collisional transitions between the lower n levels and all other levels (or disable by returning zero)
// This can prevent fully disconnected levels, whose NLTE populations cannot be determined
constexpr int NLEVELS_REQUIRETRANSITIONS(int Z, int ionstage) {
  return ((Z == 26 || Z == 28) && ionstage >= 1) ? 80 : 0;
}

// if uniform pellet energies are not used, a uniform decay time distribution is used with scaled packet energies
constexpr bool UNIFORM_PELLET_ENERGIES;

// directly calculate collisional heating from ion level populations instead of using estimators
constexpr bool DIRECT_COL_HEAT;

// INITIAL PACKETS will seed the cells on the first timestep at tmin with K-packets
// representing decay energy from t_model to tmin, and,
// if USE_MODEL_INITIAL_ENERGY is true, also the snapshot energy at t_model
constexpr bool INITIAL_PACKETS_ON;

// allows non-zero energy density at time t_model using q column in model.txt
// INITIAL_PACKETS_ON must be true to make use of this
constexpr bool USE_MODEL_INITIAL_ENERGY;

// Rate coefficients
constexpr int TABLESIZE;
constexpr double MINTEMP;
constexpr double MAXTEMP;

// temperature for which total ion recombination rate are calibrated to input data (recombrates.txt)
constexpr double RECOMBCALIBRATION_T_ELEC;

// Polarisation for real packets
constexpr bool DIPOLE;

// Compute polarisation and write specpol.out, emissionpol.out, absorptionpol.out
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

// Solve the NLTE level population equations with the Grassmann-Taksar-Heyman (GTH) state-elimination
// algorithm instead of the default LU decomposition of the rate matrix with a normalisation row.
// GTH exploits the fact that the assembled rate matrix is the transpose of a Markov-chain generator:
// it never reads the diagonal and performs no subtractions, giving populations with componentwise
// relative accuracy even when they span many orders of magnitude, with no matrix equilibration,
// balance vector, or iterative refinement. Elements with FORCE_SAHA_ION_BALANCE always use the LU
// solver, because the Saha constraint rows break the generator structure that GTH requires.
constexpr bool NLTE_USE_GTH_SOLVER;

// Solve the ionisation balance and the thermal balance time-dependently from this timestep on. This
// follows the SUMO code (Pognan, Jerkstrand & Grumer 2022, MNRAS, 510, 3806, eqs. 8 and 17). No value
// keeps the statistical equilibrium and the steady-state thermal balance for the whole run. With a value,
// the ion populations of the NLTE elements get a backward Euler time term from the previous grid update.
// The electron temperature of every cell gets the same time term. The excitation inside each ion stays
// in statistical equilibrium.
//
// The populations and the temperature are fractions of the element and of the total electron number.
// The expansion and the radioactive decay then add no terms. A decay daughter atom takes the current
// ionisation distribution of its element. A cell without a previous solution uses the steady-state
// equations for one timestep. This applies to the first NLTE timestep of the cell, to a LTE or thick
// timestep, and to an element that fell back to LTE.
//
// The time-dependent timesteps always use the LU solver, because the time term breaks the generator
// structure that NLTE_USE_GTH_SOLVER needs. The steady-state timesteps can still use GTH. Elements with
// FORCE_SAHA_ION_BALANCE and the elements without NLTE levels keep their equilibrium ionisation balance.
// The error of the backward Euler step is first order in width/mid. Keep width/mid at 0.1 or less.
//
// The k-packets carry the same energy budget. Each time a k-packet selects a cooling process, its energy gets
// the factor 1 - (c_adiabatic + c_heatcapacity) / heating of its cell (see kpkt.h). This factor applies in
// every timestep with a thermal balance, also without this option. With this option, the stored thermal
// energy also stays out of the radiation field. A gas that cools releases its stored energy into the
// packets, and the factor is then above 1. The code removes no k-packet.
constexpr std::optional<int> NLTE_TIME_DEPENDENT_FIRST_TIMESTEP;

// non-thermal ionisation
constexpr bool NT_ON;

// use the detailed Spencer-Fano solver instead of the work function approximation (only works if NT_ON)
constexpr bool NT_SOLVE_SPENCERFANO;

// NB: the energy grid of the Spencer-Fano solution vector is not an artisoptions.h option. SFPTS (the number of
// energy points) and SF_EMIN / SF_EMAX (the grid limits in eV) are set at the top of nonthermal.cc and apply to
// every preset.

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
// associated with electron shells.
constexpr int NT_MAX_AUGER_ELECTRONS;

// add the Auger electron term to the Spencer-Fano equation
constexpr bool SF_AUGER_CONTRIBUTION_ON;

constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT;

// when converting mass fraction to a number density, use a mean atomic mass
// calculated from the nuclear composition (plus stable component),
// rather than just from the compositiondata.txt values
constexpr bool USE_CALCULATED_MEANATOMICWEIGHT;

constexpr bool WRITE_EMISSIONABSORPTION_SPEC_AT_END;

// track escaped gamma-ray packets and write gamma_light_curve.out
constexpr bool KEEP_ESCAPED_GAMMAS;

// thermalisation scheme for non-thermal particles (positrons, electrons, alphas). INSTANTFULLDEPOSITION
// instantly deposits all particle energy. TIMEDEPENDENT uses time-dependent Monte Carlo transport.
// TIMEDEPENDENT_WITH_ADIABATIC_LOSS is the same as TIMEDEPENDENT but also includes adiabatic (expansion)
// losses: particle energy decreases as E/t in addition to collisional losses, and the packet comoving
// energy e_cmf is reduced by the expansion factor ts/t_new.
// BARNES, and WOLLAEGER use analytic thermalisation efficiency functions.
// TIMEDEPENDENTWITHGAMMAPRODUCTS also replaces the instant "gamma" deposition of Compton electrons and pair-production particles
// with separate handling as particle deposition
constexpr ParticleThermalisationScheme PARTICLE_THERMALISATION_SCHEME;

// thermalisation scheme for gamma-ray photons.
// FREQUENCYDEPENDENT is time-dependent gamma-ray transport (with frequency-dependent opacities if GAMMA_USE_KAPPA_GREY is not set)
// BARNES, WOLLAEGER, and GUTTMAN use analytic thermalisation efficiencies to instantly deposit some fraction of gamma energy.
constexpr GammaThermalisationScheme GAMMA_THERMALISATION_SCHEME;

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

// The bound-free cooling coefficient of each (level, target) continuum is normalised per population of the upper-ion
// target level (like the spontaneous recombination coefficient alpha_sp). Set true to multiply it by that target
// level's population. Set false (classic ARTIS behaviour) to multiply by the whole upper ion population instead,
// which treats the ion as if it were entirely in the continuum's target levels: a level with several targets shares
// the ion population among them by their LTE fractions, and a single-target level gets the whole ion population.
constexpr bool BFCOOLING_USELEVELPOPNOTIONPOP;

// set true to calculate and use expansion opacities instead of line-by-line in non-grey mode
constexpr bool RPKT_USE_EXPANSION_OPACITIES;

// set true to calculate and use expansion opacities instead of line-by-line for virtual packets
constexpr bool VPKT_USE_EXPANSION_OPACITIES;

// Optionally replace macroatom with a thermalisation probability P (where 1 - P is probability of scattering).
constexpr std::optional<float> RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY;

// For cells in grey mode, select a method of calculating opacity:
//   - FEGROUP_APPROX: Fe-group line expansion opacity scaled to local composition (default).
//   - TANAKA2020_ELECTRONFRAC: opacity parametrised using the Tanaka et al. (2020) fit to electron fraction (Ye).
//   - JUST2022_TEMP_LANTHANIDEFRAC: opacity parametrised using Just et al. (2022) fit to temperature and lanthanide fraction.
constexpr RpktGreyType RPKT_GREY_TYPE = RpktGreyType::FEGROUP_APPROX;

// Use XCOM data for gamma photoionisation instead of Si+Fe Equation 2 of Ambwani & Sutherland (1988), Veigele (1973)
constexpr bool USE_XCOM_GAMMAPHOTOION;

// Override frequency-dependent gamma-ray opacity with a grey opacity [cm^2/g]
constexpr std::optional<double> GAMMA_USE_KAPPA_GREY;

// Include charge transfer reactions in the NLTE population solver. The published fits come from
// data/chargetransfer.txt, which holds reactions with hydrogen and helium. The code generates
// estimates at startup for the other electron captures from a neutral donor (see chargetransfer.cc).
// A singly charged ion gets a flat rate: 1e-12 cm3/s, the median of the tabulated rates, for an
// energy release up to 4 eV, and the radiative floor of 1e-14 cm3/s above it. An ion with a charge
// of two or more gets a multichannel Landau-Zener estimate, with the levels of the lower ion as the
// capture channels. The reverse rates come from detailed balance. The rates enter the NLTE rate matrix as
// per-ion coefficients between neighbouring ion stages. A reaction is active only when both
// elements have NLTE levels and a free ionisation balance, so that both sides of the reaction get
// their transition and the total ionic charge stays constant. The solver does not add the reaction
// heat to the thermal balance.
constexpr bool ENABLE_CHARGE_TRANSFER_REACTIONS;

// Use microclumping, which enhances collisional (de)excitation, collisional ionisation, collisional recombination,
// radiative recombination, collisional capture, stimulated recombination, free-free heating, free-free cooling by the clumping factor
// (including in the photoionisation-equilibrium ionisation balance for elements without NLTE levels; the Saha LTE
// ion balance is unaffected). Note: with USE_LUT_PHOTOION, the stimulated recombination correction inside the
// tabulated photoionisation coefficients cannot include the per-cell clumping factor; the enhancement applies to the
// direct bound-free integrals and the packet opacity departure ratios.
constexpr bool USE_MICROCLUMPING;

// Calculate clumping factors based on time and radial velocity
// Will be passed globals::timesteps[nts].mid and the cell's mean radial velocity
// (grid::get_modelcell_mean_radial_pos_tmin(mgi) / globals::tmin) [cm/s]
constexpr float clumping_factor([[maybe_unused]] double tmid, [[maybe_unused]] double rad_vel) { return 1.; }
```
