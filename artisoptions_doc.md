```c++
// The number of energy packets of all MPI ranks together. Each rank gets an equal share, and the first ranks also
// get one packet each of the remainder. The OpenMP threads of a rank share the packets of that rank.
constexpr std::int64_t NUM_PACKETS;

// Set GridType::CARTESIAN3D to map a 1D or a 2D model onto a 3D Cartesian grid. No value keeps the grid type of
// model.txt. GridType::SPHERICAL1D and GridType::CYLINDRICAL2D are valid only when they equal the grid type of
// model.txt. Any other combination stops the run.
constexpr std::optional<GridType> GRID_TYPE_OVERRIDE;

// The size of the 3D Cartesian grid for a 1D or a 2D model. A 3D model.txt sets the size, and the code ignores
// these values.
constexpr int CUBOID_NCOORDGRID_X;
constexpr int CUBOID_NCOORDGRID_Y;
constexpr int CUBOID_NCOORDGRID_Z;

// On a 2D cylindrical or a 3D Cartesian grid, empty every propagation cell with a centre outside the sphere of
// radius vmax * tmin. A packet that enters a cell with an inner corner outside the sphere escapes at once.
constexpr bool FORCE_SPHERICAL_ESCAPE_SURFACE;

// The largest iteration number of the NLTE/Te/Spencer-Fano loop of a cell. The first pass is iteration 0, so a
// cell makes at most NLTE_TE_NNE_MAXITER + 1 passes.
constexpr int NLTE_TE_NNE_MAXITER;

// Anderson acceleration of the NLTE/Te/Spencer-Fano iteration. The iteration maps the state of one pass to the
// state of the next pass. The state holds the logarithms of the electron temperature, of the electron density,
// and of the population of each significant ion (see NLTE_SIGNIFICANT_ION_FRACTION in nltepop.h). The charge
// transfer reactions move population between two elements and hold nne, so the ion populations make that mode
// visible to the accelerator. update_grid.cc sets the maximum number of ions in the state and the maximum depth.
//
// With the acceleration on, the accelerator combines the last iterates so that the residual of the combination
// is minimal. The next pass then starts from that combination. False keeps the plain successive substitution.
//
// The accelerator forgets its history when an element falls back to LTE or changes its solved ion range, because
// the map is then discontinuous. It also forgets its history when the set of significant ions changes, because
// the state is then different. A new non-thermal solution also clears the history. The loop rejects five kinds
// of step:
// - a step outside [MINTEMP, MAXTEMP];
// - a step with an electron density below MINPOP;
// - a step with an electron density above the total electron density of the cell;
// - a step with an ion population below MINPOP or above the population of its element;
// - a step away from the map output that is larger than twice the residual.
// The convergence test comes before the injection, so a converged cell keeps the output of a plain pass. A cell
// that reaches NLTE_TE_NNE_MAXITER also ends with a plain pass.
constexpr bool NLTE_TE_NNE_USE_ANDERSON_ACCEL;

// The relative tolerance of the NLTE/Te/Spencer-Fano iteration for T_e and nne. Without the acceleration, the
// loop stops when the change of each value between two passes is at or below the tolerance. With the
// acceleration, the loop also needs an estimate of the remaining error at or below the tolerance.
//
// The estimate is change * rho / (1 - rho), where rho is the ratio of the last two changes of that value,
// limited to [0.5, 0.95]. The first sample after a reset uses rho = 0.5. T_e and nne each get an estimate, and
// the larger one counts. A pass with a changed map skips the estimate. With the charge transfer reactions or
// the acceleration, the same value is the tolerance of the ion population test.
constexpr double NLTE_TE_NNE_RELTOL;

// The number of excited levels of the ion in full NLTE. The ground state and the superlevel do not count.
constexpr int ION_NLEVELS_EXCITED_NLTE(int element_z, int ionstage);

// Use the radiation temperature T_J instead of T_e in the Boltzmann factor of the excitation, also for the
// sublevels of a superlevel. With false, the T_e finder also recomputes the photoionisation rates of the ground
// continua of the elements without NLTE levels when T_e moves by more than 10 percent.
constexpr bool LTEPOP_EXCITATION_USE_TJ;

// Force the Saha ionisation balance of an element. The NLTE solver replaces the ground level row of each ion
// above the lowest solved ion with the Saha constraint, and the elements without NLTE levels use the Saha phi
// function. The element then gets no charge transfer reaction, no time-dependent ionisation term, and always
// the LU solver.
constexpr bool FORCE_SAHA_ION_BALANCE(int element_z);

// Keep only one level, and no transition, for the highest ion stage of each element.
constexpr bool SINGLE_LEVEL_TOP_ION;

// Add a collisional transition with no radiative rate between each of the lowest n levels and every other level
// of the ion where the data has none. A level without a transition has no NLTE solution. Return zero for no
// added transitions.
constexpr int NLEVELS_REQUIRETRANSITIONS(int element_z, int ionstage);

// True: every pellet gets the same energy, and the decay time of each pellet is a sample of the decay curve.
// False: the decay times are uniform between the first decay time and tmax, and the energy of each pellet
// follows the decay power at its time.
constexpr bool UNIFORM_PELLET_ENERGIES;

// Sum the collisional de-excitation heating from the level populations instead of the Monte Carlo estimator.
// The estimator also holds the collisional recombination heating. The sum does not.
constexpr bool COL_HEAT_FROM_LEVELPOPS;

// Seed the cells at tmin with k-packets that carry the decay energy from t_model to tmin and the snapshot
// energy at t_model (the q column of model.txt). The expansion from the decay time to tmin reduces each energy.
constexpr bool INITIAL_PACKETS_ON;

// The number of temperature points of the rate coefficient tables, spaced in log T between MINTEMP and MAXTEMP.
constexpr int RATECOEFF_TABLESIZE;

// The limits [K] of T_e, T_J, the initial temperature, and the T_R of the whole spectrum. The code clamps them
// into this range. The rate coefficient tables cover this range, and the accelerator rejects a step outside it.
// The bins of the radiation field have their own T_R limits in radfield.cc.
constexpr double MINTEMP;
constexpr double MAXTEMP;

// The temperature [K] at which the code calibrates the total recombination rate of each ion to
// recombrates.txt, if that file exists.
constexpr double RECOMBCALIBRATION_T_ELEC;

// Sample the electron scattering direction of a real packet from the dipole phase function. False gives
// isotropic scattering. A virtual packet always uses the dipole function.
constexpr bool DIPOLE;

// Track the Stokes parameters and write specpol.out, emissionpol.out, and absorptionpol.out.
constexpr bool POL_ON;

// Enable the virtual packets that vpkt.txt sets up. This needs POL_ON.
constexpr bool VPKT_ON;

// Write a line to a vpackets_<rank>.out file for each emission of a real packet in a thin cell whose virtual
// packets escape in at least one observer direction of vpkt.txt. The line holds the arrival time, the frequency,
// and the energy of the contribution to each direction. This needs VPKT_ON.
constexpr bool VPKT_WRITE_CONTRIBS;

// The lower bound of the level populations, the ion populations, and the electron density nne [cm^-3]. A level
// population from the NLTE solver has no lower bound, and an absent element has zero population.
constexpr double MINPOP;

// The frequency limits of the UVOIR spectra and of the blackbody sampling of a k-packet [Hz]
constexpr double NU_MIN_R;
constexpr double NU_MAX_R;

// Take the photoionisation cross section from the table point at or below the frequency, as classic ARTIS did,
// instead of a linear interpolation.
constexpr bool PHIXS_CLASSIC_NO_INTERPOLATION;

// Fit a dilute blackbody to each frequency bin of the radiation field, in addition to the fit of the whole
// spectrum. The fit of the whole spectrum stays the fallback for a bin without a fit. Set USE_LUT_PHOTOION to
// false with this option, because the tables assume a Planck function. Nothing checks this.
constexpr bool MULTIBIN_RADFIELD_MODEL_ON;

// The number of bins, including the T_e superbin
constexpr int RADFIELDBINCOUNT;

// The first timestep at which radfield() reads the binned radiation field. It must be at or after the first
// NLTE timestep.
constexpr int FIRST_NLTE_RADFIELD_TIMESTEP;

// The frequency range of the regular bins [Hz], e.g. CLIGHT / (lambda[Angstrom] * 1e-8)
constexpr double RADFIELDBINS_NU_MIN;
constexpr double RADFIELDBINS_NU_MAX;

// The upper frequency of the T_e superbin, the last bin above RADFIELDBINS_NU_MAX. The temperature of this bin
// is the electron temperature of the cell. Only its dilution factor comes from the fit.
constexpr double RADFIELDBINS_T_E_SUPERBIN_NU_MAX;

// Store the Jb_lu estimators of the lines that radfield::init() selects.
constexpr bool DETAILED_LINE_ESTIMATORS_ON;

// Store the detailed bound-free rate estimators. This needs USE_LUT_PHOTOION false.
constexpr bool DETAILED_BF_ESTIMATORS_ON;

// Select the continua that the detailed bound-free estimators track. Only used with DETAILED_BF_ESTIMATORS_ON.
constexpr bool LEVEL_HAS_BFEST(int element_z, int ionstage, int level);

// Use the detailed bound-free estimators from this timestep on, inclusive. Only used with
// DETAILED_BF_ESTIMATORS_ON.
constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP;

// Take the photoionisation rate coefficient from a table for a Planck radiation field, instead of an integral
// over the radiation field of the cell. The table value gets the dilution factor W of the cell and the ratio of
// the Monte Carlo estimator to the blackbody rate of the nearest ground continuum.
constexpr bool USE_LUT_PHOTOION;

// Store a bound-free heating estimator for each ground continuum. Multiply the analytic bound-free heating of
// each level by the ratio of the estimator to the analytic rate of the nearest ground continuum.
constexpr bool USE_ION_BFHEATING_ESTIMATORS;

// Write the heating and the cooling rate of each ion to the estimators file, in addition to the totals. A
// bound-free term belongs to the lower ion of the continuum, and a collisional term belongs to the ion that
// holds the level. The free-free heating estimator holds no per-ion information, so each ion gets the share
// that it has in the free-free opacity. The per-ion collisional heating needs COL_HEAT_FROM_LEVELPOPS. A cell
// without a thermal balance gets no per-ion values.
constexpr bool WRITE_ION_HEATING_COOLING_RATES;

// Reject an NLTE solution with one of these faults:
// - a population that is not finite;
// - a ground population below MINPOP;
// - a population below -MINPOP;
// - a population inversion above STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL.
// Without this option, the solver replaces a negative population with the Boltzmann population. That
// population can be much larger, and a partition function can then overflow.
constexpr bool STRICT_POPULATION_CHECKING;

// After a failed NLTE solution of an element, remove the highest ion, or the lowest ion when the highest ion
// cannot go, and solve again. An ion can go only when its population is small (see
// NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION). Without this option, or when no ion can go,
// the whole element falls back to LTE.
constexpr bool NLTE_LIMIT_ION_STAGES_AFTER_FAILURE;

// The population of a level per statistical weight, relative to the same ratio of the ground level, above
// which STRICT_POPULATION_CHECKING rejects the solution
constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL;

// The same ratio above which the solver logs a warning
constexpr float STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING;

// NLTE_LIMIT_ION_STAGES_AFTER_FAILURE removes an ion only when its ground population, each NLTE excited
// population, and its superlevel population are each below this fraction of the element population.
constexpr double NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION;

// Solve the NLTE rate matrix with the Grassmann-Taksar-Heyman (GTH) state elimination instead of the LU
// decomposition with a normalisation row. The rate matrix is the transpose of a Markov chain generator. GTH
// never reads the diagonal and makes no subtraction, so each population gets a relative accuracy, also when the
// populations span many orders of magnitude. GTH needs no equilibration, no balance vector, and no iterative
// refinement. An element with FORCE_SAHA_ION_BALANCE or a time-dependent timestep (see
// NLTE_TIME_DEPENDENT_FIRST_TIMESTEP) uses the LU solver, because the Saha rows and the time term break the
// generator structure. The steady-state timesteps still use GTH.
constexpr bool NLTE_USE_GTH_SOLVER;

// Solve the ionisation balance and the thermal balance with a time term from this timestep on, as the SUMO code
// does (Pognan, Jerkstrand & Grumer 2022, MNRAS, 510, 3806-3837, doi:10.1093/mnras/stab3674, eqs. 8 and 17).
// No value keeps the statistical equilibrium and the steady-state thermal balance for the whole run. With a
// value, the ion populations of each NLTE element and the electron temperature of each cell with a thermal
// balance get a backward Euler term from the previous grid update. The excitation inside each ion stays in
// statistical equilibrium.
//
// The solver stores the ion populations as fractions of the element population and rebuilds the previous nne
// from these fractions and the current element densities. The expansion and the radioactive decay then add no
// terms, and a decay daughter atom takes the current ionisation distribution of its element. A cell without a
// previous solution uses the steady-state equations for one timestep. This applies to the first NLTE timestep
// of the cell, to the timestep after an LTE or thick timestep, and to an element that fell back to LTE.
//
// A time-dependent timestep always uses the LU solver, because the time term breaks the generator structure
// that NLTE_USE_GTH_SOLVER needs. Elements with FORCE_SAHA_ION_BALANCE and elements without NLTE levels keep
// their equilibrium ionisation balance. The error of the backward Euler step is first order in width/mid. Keep
// width/mid at 0.1 or less.
//
// The k-packets carry the same energy budget. Each time a k-packet selects a cooling process, its energy gets
// the factor 1 - (c_adiabatic + c_heatcapacity) / heating of its cell, limited to [0, 100] (see kpkt.cc). A
// cell without heating uses the factor 1. The factor applies in every timestep with a thermal balance, also
// without this option. With this option, the
// stored thermal energy also stays out of the radiation field. A gas that cools releases its stored energy
// into the packets, and the factor is then above 1. The code removes no k-packet.
constexpr std::optional<int> NLTE_TIME_DEPENDENT_FIRST_TIMESTEP;

// How the code deposits the energy of the non-thermal leptons.
// NT_OFF: no non-thermal ionisation.
// NT_SPENCERFANO: the Spencer-Fano solution. It also gives the non-thermal excitation rates for the NLTE
// population solver, the macroatom, and the NTLEPTON packets.
// NT_AXELRODAPPROX: the work function approximation of Axelrod (1980, PhD thesis, University of California,
// Santa Cruz). The energy fractions are then 0.03 for the ionisation and 0.97 for the heating, with no
// excitation rates.
constexpr NonThermalScheme NT_SCHEME;

// The energy grid of the Spencer-Fano solution is not an option of artisoptions.h. SFPTS (the number of energy
// points) is in nonthermal.h, and SF_EMIN and SF_EMAX (the grid limits in eV) are at the top of nonthermal.cc.
// They apply to every preset.

// Reuse a Spencer-Fano solution for at most this many timesteps after the timestep of the solution. 0 reuses a
// solution only within the NLTE iterations of the same timestep. A negative value solves at every iteration of
// every timestep.
constexpr int SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS;

// A change of nne per ion (nne divided by the total ion density) since the last solution at or above this
// fraction, e.g. 0.5 for 50 percent, also triggers a solution.
constexpr double NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS;

// Include non-thermal excitation only from the lowest NTEXCITATION_MAXNLEVELS_LOWER levels of an ion and to
// its lowest NTEXCITATION_MAXNLEVELS_UPPER levels, because these transitions slow the solver. A zero in either
// value includes no transition.
constexpr int NTEXCITATION_MAXNLEVELS_LOWER;
constexpr int NTEXCITATION_MAXNLEVELS_UPPER;

// The number of stored non-thermal excitation rates. The solver keeps the transitions with the largest
// deposition fractions. A transition outside the list gets no excitation rate, and an NTLEPTON packet that
// selects one becomes a k-packet.
constexpr int MAX_NT_EXCITATIONS_STORED;

// Divide by the valence shell potential of the ion, instead of the potential of each shell, in the effective
// ionisation potential and the ionisation rates. Not compatible with NT_MAX_AUGER_ELECTRONS above zero.
constexpr bool NT_USE_VALENCE_IONPOTENTIAL;

// The maximum number of Auger electrons that one impact ionisation releases, from the shell probabilities.
// Zero gives one electron per ionisation.
constexpr int NT_MAX_AUGER_ELECTRONS;

// Add the source term of the Auger electrons to the Spencer-Fano equation.
constexpr bool SF_AUGER_CONTRIBUTION_ON;

// Use the full relativistic Doppler factor instead of the first-order 1 - n.v/c. The line resonance distance
// and the walk over the expansion opacity bins then use a linear interpolation of the frequency along the path.
constexpr bool USE_RELATIVISTIC_DOPPLER_SHIFT;

// Convert a mass fraction to a number density with the mean atomic mass of the nuclear composition of the cell,
// including the stable component, instead of the mass in compositiondata.txt.
constexpr bool USE_CALCULATED_MEANATOMICWEIGHT;

// Keep the escaped gamma-ray packets in the packet files, and write gamma_light_curve.out and gamma_spec.out.
constexpr bool KEEP_ESCAPED_GAMMAS;

// The thermalisation of the non-thermal particles (positrons, electrons, and alpha particles):
// - INSTANTFULLDEPOSITION deposits the particle energy at once;
// - TIMEDEPENDENT transports the particles with the Monte Carlo method;
// - TIMEDEPENDENT_WITH_ADIABATIC_LOSS adds the adiabatic loss rate E/t to the collisional loss rate. Only the
//   collisional share of the lost energy heats the gas;
// - TIMEDEPENDENTWITHGAMMAPRODUCTS also transports the electrons and positrons from Compton scattering,
//   photoelectric absorption, and pair production, instead of an instant deposition;
// - BARNES and WOLLAEGER use analytic thermalisation efficiencies (Barnes, Kasen, Wu & Martínez-Pinedo 2016,
//   ApJ, 829, 110, doi:10.3847/0004-637X/829/2/110; Wollaeger, Korobkin, Fontes, Rosswog, Even & Fryer 2018,
//   MNRAS, 478, 3298-3334, doi:10.1093/mnras/sty1018).
constexpr ParticleThermalisationScheme PARTICLE_THERMALISATION_SCHEME;

// The thermalisation of the gamma-ray photons. FREQUENCYDEPENDENT transports the gamma rays with the Monte
// Carlo method, with frequency-dependent opacities unless GAMMA_USE_KAPPA_GREY has a value. BARNES, WOLLAEGER,
// and GUTTMAN deposit a fraction of the gamma energy at once, from an analytic thermalisation efficiency. See
// PARTICLE_THERMALISATION_SCHEME for the first two references. GUTTMAN: Guttman, Shenhar, Sarkar & Waxman 2024,
// MNRAS, 533, 994-1011, doi:10.1093/mnras/stae1795.
constexpr GammaThermalisationScheme GAMMA_THERMALISATION_SCHEME;

// The timestep scheme: LOGARITHMIC, CONSTANT, LOGARITHMIC_THEN_CONSTANT, or CONSTANT_THEN_LOGARITHMIC. The two
// hybrid schemes need FIXED_TIMESTEP_WIDTH and TIMESTEP_TRANSITION_TIME, both in days. The number of constant
// timesteps must be below the total number of timesteps, or the run stops.
constexpr TimeStepSizeMethod TIMESTEP_SIZE_METHOD;

// The maximum width of a constant timestep [days], for the hybrid schemes. The code divides the constant
// interval into equal timesteps of at most this width.
constexpr double FIXED_TIMESTEP_WIDTH;

// The time of the change between the two schemes [days], for the hybrid schemes
constexpr double TIMESTEP_TRANSITION_TIME;

// The bound-free cooling coefficient of each (level, target) continuum is per population of the target level of
// the upper ion, like the spontaneous recombination coefficient alpha_sp. True multiplies it by the population
// of that target level. False (classic ARTIS) multiplies it by the population of the whole upper ion. A level
// with several targets then shares the ion population among them by their LTE fractions, and a level with one
// target gets the whole ion population.
constexpr bool BFCOOLING_USELEVELPOPNOTIONPOP;

// Use expansion opacities instead of line-by-line opacities for the real packets in the cells that are not
// thick. Not compatible with VPKT_ON.
constexpr bool RPKT_USE_EXPANSION_OPACITIES;

// Use expansion opacities instead of line-by-line opacities for the virtual packets.
constexpr bool VPKT_USE_EXPANSION_OPACITIES;

// The line weight in the expansion opacity of each wavelength bin of RPKT_USE_EXPANSION_OPACITIES,
// VPKT_USE_EXPANSION_OPACITIES, and RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY. A line with the Sobolev optical
// depth tau adds (lambda / delta_lambda) * weight to the sum of its bin:
// - EXPANSION: 1 - exp(-tau) (Eastman & Pinto 1993, ApJ, 412, 731-751, doi:10.1086/172957).
// - LINEBINNEDCAPPED: min(1, tau), the line-binned opacity with a limit of 1 for each line.
// - LINEBINNED: tau, the line-binned opacity (Fontes, Fryer, Hungerford, Wollaeger & Korobkin 2020, MNRAS, 493,
//   4143-4171, doi:10.1093/mnras/staa485).
constexpr ExpansionOpacityMethod EXPANSION_OPACITY_METHOD;

// Replace the macroatom with a thermalisation probability P for each bound-bound absorption, and a scattering
// in the absorbing line with probability 1 - P. Every k-packet in a cell that is not thick then emits a blackbody
// spectrum weighted with the sum of the expansion opacity and the free-free opacity. By Kirchhoff's law, the
// emission type is free-free or a line of the wavelength bin, in proportion to their opacities. The code therefore
// computes the expansion opacities also without RPKT_USE_EXPANSION_OPACITIES. A thick cell samples a plain
// Planck function. No value keeps the macroatom.
constexpr std::optional<float> RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY;

// The grey opacity of a thick cell:
// - FEGROUP_APPROX: 0.1 cm^2/g times (0.9 X + 0.1) / (0.9 <X> + 0.1), where X is the initial Fe-group mass
//   fraction of the cell and <X> is the mean of the model;
// - TANAKA2020_ELECTRONFRAC: a fit to the electron fraction Y_e (Tanaka, Kato, Gaigalas & Kawaguchi 2020,
//   MNRAS, 496, 1369-1392, doi:10.1093/mnras/staa1576);
// - JUST2022_TEMP_LANTHANIDEFRAC: a fit to the temperature and the lanthanide fraction (Just, Kullmann,
//   Goriely, Bauswein, Janka & Collins 2022, MNRAS, 510, 2820-2840, doi:10.1093/mnras/stab3327). The code
//   recomputes this opacity at each timestep.
constexpr RpktGreyType RPKT_GREY_TYPE;

// Use the XCOM table for the gamma-ray photoelectric absorption instead of the Si and Fe power laws of
// Ambwani & Sutherland 1988, ApJ, 325, 820-827, doi:10.1086/166052, eq. 2, after Veigele 1973, Atomic Data and
// Nuclear Data Tables, 5, 51-111, doi:10.1016/S0092-640X(73)80015-4.
constexpr bool USE_XCOM_GAMMAPHOTOION;

// Replace the frequency-dependent gamma-ray opacity with this grey opacity [cm^2/g]. No value keeps the
// frequency-dependent opacity.
constexpr std::optional<double> GAMMA_USE_KAPPA_GREY;

// Include the charge transfer reactions in the NLTE population solver. The published fits come from
// data/chargetransfer.txt, which holds reactions with hydrogen and helium. The code estimates the other electron
// captures from a neutral donor at startup (see chargetransfer.cc). A singly charged ion gets a flat rate of
// 1e-12 cm3/s, the median of the tabulated rates, for an energy release up to 4 eV. Above 4 eV it gets the
// radiative floor of 1e-14 cm3/s. An ion with a charge of two or more gets a multichannel Landau-Zener
// estimate, with the levels of the lower ion as the capture channels.
//
// The reverse rates come from detailed balance. The rates enter the NLTE rate matrix as per-ion coefficients
// between neighbouring ion stages. A reaction is active only when both elements have NLTE levels and a free
// ionisation balance. Both sides of the reaction then get their transition, and the total ionic charge stays
// constant. The solver adds no reaction heat to the thermal balance.
constexpr bool ENABLE_CHARGE_TRANSFER_REACTIONS;

// Multiply these rates by the clumping factor of the cell:
// - collisional excitation, de-excitation, ionisation, and recombination;
// - radiative and stimulated recombination, and collisional capture;
// - free-free heating and cooling, and bound-free cooling;
// - the charge transfer reactions;
// - the photoionisation equilibrium of the elements without NLTE levels. The Saha balance is unaffected.
// With USE_LUT_PHOTOION, the stimulated recombination correction inside the tabulated photoionisation
// coefficients cannot include the clumping factor of the cell. The stimulated recombination then gets the factor
// through the direct bound-free integrals and through the departure ratios of the packet opacity.
constexpr bool USE_MICROCLUMPING;

// The clumping factor of a cell from the time and the radial velocity. The code passes
// globals::timesteps[nts].mid and grid::get_modelcell_mean_radial_pos_tmin(mgi) / globals::tmin [cm/s]. The
// result must be finite and at least 1. A value of 1 means no clumping.
constexpr float clumping_factor(double tmid, double rad_vel);
```
