# Release notes

This page lists the changes that a user sees, for each release of about the last year. The newest release comes
first. The page leaves out the changes of the internal structure, of the continuous integration, and of the
dependencies. The [releases page](https://github.com/artis-mcrt/artis/releases) on GitHub holds the full list of
each release, and the run times of two reference models.

A release is a git tag, e.g. `v2026.5.15`. The log file of rank 0 gives the version of a run.

Each release has these sections, when it has such changes:

- **The numerical results**: the changes that alter the output of a simulation. Each entry names the end-to-end
  tests whose reference checksums changed (see [Testing](testing.md)).
- **Input files** and **Output files**: the changes of a file format or of the files that a run reads and writes.
- **Compile-time options**: the options that the release added, renamed, or removed (see
  [Compile-time options](options.md)).
- **Build, programs, and scripts**: the other changes that a user sees.

The entries use these names for the groups of tests:

- classic: `classicmode_1d_3dgrid` and `classicmode_3d`;
- kilonova: the tests whose name starts with `kilonova`, without the time-dependent NLTE test;
- nebular: `nebular_1d_3dgrid` and `nebular_1d_3dgrid_limitbfest`;
- nltephotospheric: `nltephotospheric_dynamic_ion_range_1d_1dgrid`;
- timedepnlte: `kilonova_1d_timedepnlte`;
- all: every end-to-end test of that time.

Each entry gives the current name of a test. The names of the 2D kilonova tests held the word `2dgrid` before
[#346](https://github.com/artis-mcrt/artis/pull/346), e.g. `kilonova_2d_2dgrid_expansionopac`.

:::{warning}
A restart file must only be consistent with the binary that wrote it. A simulation in progress cannot continue with
a newer version of the code. Complete the simulation with the version that started it.
:::

## Changes after v2026.5.15

The `main` branch holds these changes. No release holds them at this time.

### The numerical results

#### Atomic masses, decay data, and deposition

- The code converts a mass number to grams with the atomic mass unit. Before, it used the mass of a hydrogen atom,
  so each number density was 0.78 percent too small. Ni56 and Co56 get their current mean lives
  ([#670](https://github.com/artis-mcrt/artis/pull/670)). *Tests: all.*
- The atomic mass unit has the CODATA 2022 value ([#679](https://github.com/artis-mcrt/artis/pull/679)).
  *Tests: all.*
- The energy columns of the decay data files come from the ENDF/B-VII.1 decay data, the actinide chains are
  connected, and the Bateman solution is correct for these chains
  ([#484](https://github.com/artis-mcrt/artis/pull/484)).
  *Tests: kilonova, nltephotospheric.*
- The code reads the values of the fission product table as yields in percent for each fission
  ([#480](https://github.com/artis-mcrt/artis/pull/480)). *Tests: kilonova.*
- These standard nuclides get their total decay Q-values ([#448](https://github.com/artis-mcrt/artis/pull/448)).
  *Tests: all.*
    - Ni57;
    - Ni56;
    - Co56;
    - Cr48;
    - V48;
    - Co57;
    - Fe52;
    - Mn52.
- The mean life of Ni57 has the correct unit conversion ([#447](https://github.com/artis-mcrt/artis/pull/447)).
  *Tests: `classicmode_3d`, kilonova.*
- Fe52 and Mn52 get one gamma-ray line each. Before, their pellets deposited the energy locally
  ([#668](https://github.com/artis-mcrt/artis/pull/668)). *Tests: all.*
- The deposition estimators count a pellet of a nuclide that has no gamma-ray spectrum
  ([#592](https://github.com/artis-mcrt/artis/pull/592)). *Tests: classic, kilonova.*
- The deposition of the spontaneous fission is the mean of all ranks, not the value of rank 0
  ([#465](https://github.com/artis-mcrt/artis/pull/465)). *Tests: kilonova.*
- The gamma-ray deposition estimator uses opacities of one frame
  ([#433](https://github.com/artis-mcrt/artis/pull/433)). *Tests: all except
  `kilonova_2d_barnesthermalisation`.*
- The pair-production cross-section below 1.5 MeV has the coefficient 0.10063 in place of 1.0063
  ([#478](https://github.com/artis-mcrt/artis/pull/478)). *Tests: all except `kilonova_2d_barnesthermalisation`.*
- The decay energy of the initial temperature has the weight of the 1/t law
  ([#631](https://github.com/artis-mcrt/artis/pull/631)). *Tests: all.*
- The kinetic energy of the ejecta uses the mean of r squared of each cell. Only the Barnes thermalisation schemes
  read this energy ([#555](https://github.com/artis-mcrt/artis/pull/555)).
  *Tests: `kilonova_2d_barnesthermalisation`.*

#### Non-thermal particles

- The Axelrod work function approximation uses the number-weighted mean atomic number
  ([#396](https://github.com/artis-mcrt/artis/pull/396)). *Tests: nebular, nltephotospheric.*
- The Lotz ionisation cross-section uses the natural logarithm in the high-energy correction
  ([#401](https://github.com/artis-mcrt/artis/pull/401)) and the relativistic beta squared
  ([#479](https://github.com/artis-mcrt/artis/pull/479)). *Tests: nltephotospheric.*
- The Spencer-Fano solver has a correct source term of the Auger electrons
  ([#450](https://github.com/artis-mcrt/artis/pull/450)). *Tests: `nebular_1d_3dgrid`, nltephotospheric.*
- The excitation cross-section from a collision strength is linear in the ratio of the hydrogen ionisation
  potential to the energy, and the Coulomb logarithm is correct at low energy
  ([#452](https://github.com/artis-mcrt/artis/pull/452)). *Tests: nebular, nltephotospheric.*
- The source integral of each row of the Spencer-Fano matrix includes the energy bin of that row
  ([#454](https://github.com/artis-mcrt/artis/pull/454)). *Tests: nebular, nltephotospheric.*
- The non-thermal excitation rates leave out the transitions below the minimum energy of the Spencer-Fano grid
  ([#459](https://github.com/artis-mcrt/artis/pull/459)). *Tests: nebular, nltephotospheric.*
- A non-thermal lepton packet selects heating, ionisation, or excitation with the fractions of the Spencer-Fano
  solution ([#453](https://github.com/artis-mcrt/artis/pull/453)). It selects the excited transition from the stored
  list with the sum of that list ([#671](https://github.com/artis-mcrt/artis/pull/671)). *Tests: nebular.*
- A non-thermal ionisation selects its channel with the energy that the Spencer-Fano solution gives to that channel
  ([#709](https://github.com/artis-mcrt/artis/pull/709)). *Tests: nebular, nltephotospheric.*
- The construction and the solution of the Spencer-Fano matrix are faster
  ([#510](https://github.com/artis-mcrt/artis/pull/510), [#520](https://github.com/artis-mcrt/artis/pull/520),
  [#522](https://github.com/artis-mcrt/artis/pull/522)). *Tests: nebular.*

#### Populations, thermal balance, and radiation field

- The NLTE population matrix gets its normalisation from an iterative Osborne balance
  ([#375](https://github.com/artis-mcrt/artis/pull/375)). *Tests: nebular, nltephotospheric.*
- The interpolation of `recombrates.txt` has the correct weights
  ([#424](https://github.com/artis-mcrt/artis/pull/424)). *Tests: nebular, nltephotospheric.*
- The detailed bound-free rate estimators have no Doppler factor
  ([#436](https://github.com/artis-mcrt/artis/pull/436)). *Tests: nebular, nltephotospheric.*
- The radiative de-excitation rate coefficient is the A-value in the optically thin limit
  ([#451](https://github.com/artis-mcrt/artis/pull/451)). *Tests: nebular, nltephotospheric.*
- The Planck mean frequency of a radiation field bin does not underflow
  ([#457](https://github.com/artis-mcrt/artis/pull/457)). *Tests: nebular, nltephotospheric.*
- The sample of the free-bound emission frequency does not underflow
  ([#461](https://github.com/artis-mcrt/artis/pull/461)). *Tests: `classicmode_3d`, nebular, nltephotospheric.*
- The NLTE matrix leaves out the transitions between two levels of the same superlevel
  ([#523](https://github.com/artis-mcrt/artis/pull/523)). *Tests: nebular, nltephotospheric.*
- With `BFCOOLING_USELEVELPOPNOTIONPOP` false, the bound-free cooling counts a level with more than one target
  only once ([#530](https://github.com/artis-mcrt/artis/pull/530)). *Tests: `kilonova_1d`, `kilonova_2d`, nebular,
  nltephotospheric.*
- Each k-packet loses the share of its energy that the thermal balance gives to the adiabatic cooling
  ([#538](https://github.com/artis-mcrt/artis/pull/538)). *Tests: classic, nebular, nltephotospheric.*
- The adiabatic cooling of the thermal balance uses the mid time of the current timestep
  ([#653](https://github.com/artis-mcrt/artis/pull/653)). *Tests: classic, timedepnlte, nebular,
  nltephotospheric.*
- Every element gets its bound-free heating. Before, an element with a mass fraction of 1 percent or less got none
  without `USE_ION_BFHEATING_ESTIMATORS` ([#593](https://github.com/artis-mcrt/artis/pull/593)).
  *Tests: timedepnlte, nebular, nltephotospheric.*
- A cell with J = 0 gets T_J = `MINTEMP` on both paths of the radiation field
  ([#594](https://github.com/artis-mcrt/artis/pull/594)). *Tests: `classicmode_3d`.*
- The renormalisation factor of the photoionisation tables stays 1 for an element that the cell does not contain
  ([#627](https://github.com/artis-mcrt/artis/pull/627)). *Tests: `classicmode_3d`.*
- The T_e solver calculates the ionisation rates again at its final T_e
  ([#623](https://github.com/artis-mcrt/artis/pull/623)), and it calculates the partition functions of the elements
  without NLTE levels at each trial T_e ([#636](https://github.com/artis-mcrt/artis/pull/636)). *Tests: timedepnlte.*
- The outer iteration of the NLTE populations, T_e, and nne is correct
  ([#542](https://github.com/artis-mcrt/artis/pull/542), [#549](https://github.com/artis-mcrt/artis/pull/549)).
  *Tests: timedepnlte.*
- Three code reviews corrected defects of the parallel code, of the NLTE solver, and of the bound-free continuum
  selection ([#441](https://github.com/artis-mcrt/artis/pull/441),
  [#469](https://github.com/artis-mcrt/artis/pull/469), [#560](https://github.com/artis-mcrt/artis/pull/560)).
  *Tests: all.*

#### Packets, grid, and propagation

- The k-packet diffusion includes the adiabatic losses ([#390](https://github.com/artis-mcrt/artis/pull/390)).
  *Tests: all except `kilonova_2d_expansionopac`.*
- The distance to a cell boundary and the distance to a line have a better numerical accuracy
  ([#398](https://github.com/artis-mcrt/artis/pull/398), [#399](https://github.com/artis-mcrt/artis/pull/399),
  [#404](https://github.com/artis-mcrt/artis/pull/404)). *Tests: classic, nebular, nltephotospheric.*
- The line list stores the Einstein B coefficients as single-precision values
  ([#421](https://github.com/artis-mcrt/artis/pull/421)). *Tests: all.*
- Every build uses the xoshiro random number generator ([#426](https://github.com/artis-mcrt/artis/pull/426)).
  *Tests: all.*
- The random isotropic direction has a better floating-point precision
  ([#434](https://github.com/artis-mcrt/artis/pull/434)). *Tests: classic, nebular, nltephotospheric.*
- The maximum path step of a packet increases with time, as the grid does
  ([#471](https://github.com/artis-mcrt/artis/pull/471)). *Tests: classic, kilonova, nltephotospheric.*
- The grid update uses the decay coefficients of each timestep again for all cells, which makes it faster
  ([#494](https://github.com/artis-mcrt/artis/pull/494), [#495](https://github.com/artis-mcrt/artis/pull/495)).
  *Tests: kilonova, nebular, nltephotospheric.*
- For a model on a grid of a different geometry, the code rescales the nuclide abundances before it calculates
  other quantities from them ([#500](https://github.com/artis-mcrt/artis/pull/500)).
  *Tests: `classicmode_1d_3dgrid`, nebular.*
- The Sobolev optical depth of a line uses the time at which the packet arrives at the line
  ([#586](https://github.com/artis-mcrt/artis/pull/586)). *Tests: all except `kilonova_2d_expansionopac`.*
- A simulation in more than one job gives the results of a simulation in one job. The restart files hold the state
  of the random number generator and the partition functions
  ([#634](https://github.com/artis-mcrt/artis/pull/634)). *Tests: all, in the checksums after the resume run.*
- The radius of a propagation cell is its midpoint radius for all grid types
  ([#681](https://github.com/artis-mcrt/artis/pull/681)). *Tests: all.*
- A cell that has no mass in the elements of `compositiondata.txt` becomes an empty cell
  ([#644](https://github.com/artis-mcrt/artis/pull/644)). *Tests: `classicmode_3d`.*
- The scatter directions near a pole are correct ([#617](https://github.com/artis-mcrt/artis/pull/617)).
  *Tests: timedepnlte, nebular, nltephotospheric.*
- Three changes of the arithmetic, with no change of the method
  ([#382](https://github.com/artis-mcrt/artis/pull/382), [#402](https://github.com/artis-mcrt/artis/pull/402),
  [#466](https://github.com/artis-mcrt/artis/pull/466)). *Tests: classic, nebular, nltephotospheric,
  `kilonova_2d_barnesthermalisation`.*

#### Opacities and virtual packets

- The range of the expansion opacity bins changed ([#389](https://github.com/artis-mcrt/artis/pull/389)).
  *Tests: `classicmode_3d`, `kilonova_2d_expansionopac`.*
- A virtual packet with binned opacities goes line by line through its first bin
  ([#394](https://github.com/artis-mcrt/artis/pull/394)), and that step leaves out the line that emitted the packet
  ([#589](https://github.com/artis-mcrt/artis/pull/589)). *Tests: `classicmode_3d`.*
- The continuum opacities of a virtual packet have the Doppler factor
  ([#438](https://github.com/artis-mcrt/artis/pull/438)), its emission probability has the relativistic beaming
  ([#439](https://github.com/artis-mcrt/artis/pull/439)), and its opacities have the correct time scaling
  ([#468](https://github.com/artis-mcrt/artis/pull/468)). *Tests: classic.*
- A virtual packet observer at a pole gets the exact direction and the exact Stokes frame
  ([#618](https://github.com/artis-mcrt/artis/pull/618)). *Tests: classic, timedepnlte, nebular.*
- The virtual packets can use the `LINEBINNED` line weight
  ([#646](https://github.com/artis-mcrt/artis/pull/646)). *Tests: `classicmode_3d`.*
- With the parameterised thermalisation, the code resets the scattering count of a packet
  ([#440](https://github.com/artis-mcrt/artis/pull/440)). *Tests: `kilonova_2d_expansionopac`.*
- The sum of the binned opacities keeps the Sobolev optical depth in double precision
  ([#491](https://github.com/artis-mcrt/artis/pull/491)). *Tests: `classicmode_3d`, `kilonova_2d_expansionopac`.*
- A second row of `transitiondata.txt` for the same pair of levels adds its A-value to the line, so the Sobolev
  opacity uses the sum ([#591](https://github.com/artis-mcrt/artis/pull/591)). *Tests: timedepnlte.*
- With `USE_RELATIVISTIC_DOPPLER_SHIFT`, a binned opacity always gets the path correction
  ([#607](https://github.com/artis-mcrt/artis/pull/607)). *Tests: `kilonova_2d_expansionopac`.*

#### The presets

A change of a preset value alters the results of each run that uses the preset without an edit.

- `kilonova_lte` uses `PARTICLE_THERMALISATION_SCHEME = TIMEDEPENDENT_WITH_ADIABATIC_LOSS`
  ([#562](https://github.com/artis-mcrt/artis/pull/562)). *Tests: `kilonova_1d`, `kilonova_2d`,
  `kilonova_2d_expansionopac`, `kilonova_2d_xcomgammaphotoion`.* The preset now has the value
  `TIMEDEPENDENT_WITH_ADIABATIC_LOSS_ULTRARELATIVISTICLIMIT`, which keeps the adiabatic loss rate of that change
  ([#669](https://github.com/artis-mcrt/artis/pull/669)).
- `nltenebular` uses `NLTE_USE_GTH_SOLVER = true` ([#703](https://github.com/artis-mcrt/artis/pull/703)).
  *Tests: nebular.*
- `nltenebular` and `nltephotospheric` use `BFCOOLING_USELEVELPOPNOTIONPOP = true`
  ([#704](https://github.com/artis-mcrt/artis/pull/704)). *Tests: nebular, nltephotospheric.*
- The T_e solver has a fractional accuracy of 1e-3 in every preset
  ([#543](https://github.com/artis-mcrt/artis/pull/543)). *Tests: classic.*
- Every preset has `KEEP_ESCAPED_GAMMAS = false` ([#629](https://github.com/artis-mcrt/artis/pull/629)). A run
  then keeps no escaped gamma packets, and it writes no gamma-ray light curve or spectrum.
- `kilonova_lte` has the values of the kilonova NLTE preset for the options that an LTE run does not use
  ([#544](https://github.com/artis-mcrt/artis/pull/544), [#629](https://github.com/artis-mcrt/artis/pull/629)).
  This change altered no test.

#### Corrections that changed no test

These corrections alter the results, or let a run complete, only for some models or some option values. The list
is not complete:

- The sample of a free-bound emission frequency stays inside the cross-section table
  ([#460](https://github.com/artis-mcrt/artis/pull/460)).
- A model with alpha decays and no `X_He4` column keeps the helium of the decays
  ([#551](https://github.com/artis-mcrt/artis/pull/551)).
- `FORCE_SPHERICAL_ESCAPE_SURFACE` removes the mass outside the sphere for a model on a grid of a different
  geometry also ([#552](https://github.com/artis-mcrt/artis/pull/552)), and it works for each value of the maximum
  velocity ([#561](https://github.com/artis-mcrt/artis/pull/561)).
- A cell that changes from thick to thin gets valid photoionisation rates
  ([#570](https://github.com/artis-mcrt/artis/pull/570)) and a valid radiation field
  ([#590](https://github.com/artis-mcrt/artis/pull/590)) in its first thin timestep.
- A bound-free absorption divides the energy between the ionisation and the thermal pool at the frequency of the
  packet, with a limit of one for the ratio ([#572](https://github.com/artis-mcrt/artis/pull/572)).
- With `PARTICLE_THERMALISATION_SCHEME = TIMEDEPENDENT_WITH_ADIABATIC_LOSS`, the heating rate and the energy that a
  particle packet deposits are correct for the adiabatic loss
  ([#482](https://github.com/artis-mcrt/artis/pull/482)). The adiabatic loss rate of that value is now relativistic,
  with the mass of the particle ([#669](https://github.com/artis-mcrt/artis/pull/669)).
- With `NT_SCHEME = NT_AXELRODAPPROX`, the thermal balance counts all of the lepton deposition as heat
  ([#654](https://github.com/artis-mcrt/artis/pull/654)).
- With `RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY`, the thermal emission includes the bound-free heating opacity
  ([#673](https://github.com/artis-mcrt/artis/pull/673)).
- With `RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY`, a blackbody emission, a thermal redistribution, and a line
  scattering each start virtual packets. A pre-k-packet in a thin cell also starts them ([#619](https://github.com/artis-mcrt/artis/pull/619)).
- The virtual packet spectra hold only the flux inside the wavelength ranges of `vpkt.txt`. The restart files of the
  virtual packets have the full precision, so the spectra do not depend on the division of a run into jobs
  ([#556](https://github.com/artis-mcrt/artis/pull/556)).
- With `NLTE_TIME_DEPENDENT_FIRST_TIMESTEP`, a k-packet that becomes radiation gets the radiative share of the
  thermal balance of its cell ([#537](https://github.com/artis-mcrt/artis/pull/537)).
- The recombination coefficient of an element without NLTE levels is correct for a photoionisation with more than
  one target level, and for calibrated data with one target level ([#529](https://github.com/artis-mcrt/artis/pull/529)).
- The Ni ions with the ion stages 23 to 28 get the collisional ionisation fits of `collion.txt`. Before, the code
  did not read the last ten rows of that file ([#649](https://github.com/artis-mcrt/artis/pull/649)).
- The Compton deposition of a gamma packet below approximately 5 keV keeps its precision ([#467](https://github.com/artis-mcrt/artis/pull/467)).
- A GPU build with more than one rank gives each rank its own packet seeds. Before, neighbouring ranks had most of
  their seeds in common ([#481](https://github.com/artis-mcrt/artis/pull/481)).
- A photoionisation or a bound-free heating rate of zero does not stop the run. The renormalisation factor of the
  estimator is then 1, with a warning in the log ([#496](https://github.com/artis-mcrt/artis/pull/496)).
- The sample of a decay time completes for a decay chain with a very long mean lifetime, e.g. U238 ([#579](https://github.com/artis-mcrt/artis/pull/579)).
- A virtual packet takes only the continuum opacity over a path that is too short to change its comoving frequency.
  Before, it added every remaining line ([#575](https://github.com/artis-mcrt/artis/pull/575)).
- With `RPKT_USE_BINNED_OPACITIES` and no `RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY`, the search for a line
  inside an opacity bin uses the time of the packet. Before, it used the midpoint of the timestep ([#571](https://github.com/artis-mcrt/artis/pull/571)).
- With `USE_MICROCLUMPING`, the ionisation balance of an element without NLTE levels uses the electron density
  inside the clumps ([#476](https://github.com/artis-mcrt/artis/pull/476)).
- With atomic data in which two ions have the same ground threshold, each ion gets its own photoionisation and
  bound-free heating estimators. Such a dataset now stops the run at the start ([#554](https://github.com/artis-mcrt/artis/pull/554)).
- A cell in which every element has only its lowest ion stage keeps the NLTE populations of an element with an NLTE
  solution ([#651](https://github.com/artis-mcrt/artis/pull/651)).
- The ion list of a cell keeps the ions that the non-thermal ionisation can make when the cell has lepton
  deposition and no gamma-ray deposition, e.g. in a kilonova model ([#503](https://github.com/artis-mcrt/artis/pull/503)).
- The number of levels in the ground term of an ion comes from the rules of LS coupling. Before, an inverted 3P
  ground term counted as one level, e.g. for O I, Ne III, S I, and Ar III ([#533](https://github.com/artis-mcrt/artis/pull/533)).
- The Grassmann-Taksar-Heyman solver stops its rescue loop after eight passes. Before, a rate that was not finite
  made the loop run without end ([#498](https://github.com/artis-mcrt/artis/pull/498)).
- A sample of a gamma-ray line, of a non-thermal ionisation target, of an Auger outcome, or of a collisional
  transition does not select an outcome with zero weight ([#463](https://github.com/artis-mcrt/artis/pull/463)).
- The Lotz cross-section gets the energy in the correct unit for the heating below the Spencer-Fano energy grid
  ([#405](https://github.com/artis-mcrt/artis/pull/405)).
- `GAMMA_THERMALISATION_SCHEME = GUTTMAN` has a corrected calculation ([#446](https://github.com/artis-mcrt/artis/pull/446)).
- With `RPKT_USE_BINNED_OPACITIES` and no `RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY`, a packet finds its position
  in the line list again from its comoving frequency after a step through the opacity bins ([#597](https://github.com/artis-mcrt/artis/pull/597)).
- A cell with only neutral ions keeps an electron density above zero in a build with the default fast-math option.
  Before, the next k-packet of such a cell stopped the run ([#628](https://github.com/artis-mcrt/artis/pull/628)).
- A model with Cu XX or Zn VIII reads its row of the Auger table. Before, such a model stopped the run ([#588](https://github.com/artis-mcrt/artis/pull/588)).
- The stored list of the non-thermal excitations holds no transition of an ion whose excitation fraction is above 1
  or not finite ([#458](https://github.com/artis-mcrt/artis/pull/458)).
- A packet that moves exactly along the z axis of a 2D cylindrical grid finds each cell boundary ([#477](https://github.com/artis-mcrt/artis/pull/477)).
- The transition table of an ion uses a stable sort, so duplicate rows of `transitiondata.txt` keep the sequence
  of the file ([#622](https://github.com/artis-mcrt/artis/pull/622)).
- For a 2D model on a Cartesian grid, the code also rescales the nuclide abundances to keep the mass of each
  nuclide ([#502](https://github.com/artis-mcrt/artis/pull/502)).
- With `ENABLE_CHARGE_TRANSFER_REACTIONS` and no Anderson acceleration, the iteration of a cell does not stop on a
  pass that changed the ion range of an element ([#558](https://github.com/artis-mcrt/artis/pull/558)).
- The NLTE solver does not accept a solution that is not finite ([#488](https://github.com/artis-mcrt/artis/pull/488)).
- The index of the NLTE population array has 64 bits. Before, a model whose number of cells with matter times its
  number of NLTE levels reached 2^31 wrote populations into the wrong cell ([#557](https://github.com/artis-mcrt/artis/pull/557)).

### Input files

- A build with libzstd reads a zstd compressed input file, e.g. `model.txt.zst`, when the plain file is absent
  ([#632](https://github.com/artis-mcrt/artis/pull/632)).
- The level indices and the cell indices of the input files can start at 0 or at 1. `sn3d` finds the first index
  from `adata.txt` and from `model.txt`, and the log messages give the 0-based index
  ([#492](https://github.com/artis-mcrt/artis/pull/492), [#687](https://github.com/artis-mcrt/artis/pull/687)).
  [Input files](inputfiles.md) gives the rule of each file.
- `Ye.txt` uses the cell numbering of `model.txt`. Before, a 0-based file moved each value by one cell
  ([#553](https://github.com/artis-mcrt/artis/pull/553)).
- `compositiondata.txt` can hold comments after a `#` character
  ([#616](https://github.com/artis-mcrt/artis/pull/616)).
- `INITIAL_PACKETS_ON` also uses the `q` column of `model.txt`, which holds the trapped radiation energy per mass
  ([#544](https://github.com/artis-mcrt/artis/pull/544), [#631](https://github.com/artis-mcrt/artis/pull/631)).
- `sn3d` sets the `nprocs_exspec` line of `input.txt` at the start of a new simulation
  ([#611](https://github.com/artis-mcrt/artis/pull/611)). Line 20 (`max_bf_continua`) has no effect
  ([#526](https://github.com/artis-mcrt/artis/pull/526)).
- A `recombrates.txt` that is not complete, or that has an unreadable row, stops the run. Before, the reader kept
  the values of the row before ([#483](https://github.com/artis-mcrt/artis/pull/483)).
- An invalid input stops the run with a message that names the file and the values
  ([#701](https://github.com/artis-mcrt/artis/pull/701)). Some examples of the new checks:
    - a cell count or a snapshot time of `model.txt` that is not positive, and an electron fraction outside 0 to 1
      ([#526](https://github.com/artis-mcrt/artis/pull/526));
    - a row of a packet file or of `abundances.txt` that is not complete
      ([#559](https://github.com/artis-mcrt/artis/pull/559));
    - a photoionisation target level that occurs twice ([#548](https://github.com/artis-mcrt/artis/pull/548)) or
      that the ion does not have ([#569](https://github.com/artis-mcrt/artis/pull/569));
    - a change of `tmin_days` or `tmax_days` before a resume run
      ([#569](https://github.com/artis-mcrt/artis/pull/569)).
- A velocity grid map of `vpkt.txt` with a time range or a wavelength range outside the ranges of the spectra stops
  the run. Before, such a map was empty or biased ([#667](https://github.com/artis-mcrt/artis/pull/667)).
  *Tests: `classicmode_1d_3dgrid`.*
- `sn3d` finds an optional input file in `data/` and in `artis/data/`, as it finds a required file, e.g. `Ye.txt`,
  `recombrates.txt`, and `autoion.txt` ([#662](https://github.com/artis-mcrt/artis/pull/662)).
- `autoion.txt` can give more than one transition for one level. The rows of one level must come together, or the
  run stops ([#663](https://github.com/artis-mcrt/artis/pull/663)).
- A level of `adata.txt` with an energy below the energy of the previous level gives a warning in the log. The
  check came with [#567](https://github.com/artis-mcrt/artis/pull/567).

### Output files

- `sn3d` writes the log and the estimator files of each job into a job folder `job_from_ts<n>`, where n is the
  first timestep of the job. A new simulation first removes the output files of the previous simulation
  ([#609](https://github.com/artis-mcrt/artis/pull/609)).
- A build with libzstd writes each `.out` file zstd compressed, e.g. `spec.out.zst`. `sn3d` writes the packet files
  into `packets/`, and the virtual packet files into `vpackets/`, `vspecpol/`, and `vpkt_grid/`
  ([#632](https://github.com/artis-mcrt/artis/pull/632)).
- `sn3d` writes all of the spectrum files that `exspec` writes, e.g. `specpol.out` and `gamma_spec.out`
  ([#598](https://github.com/artis-mcrt/artis/pull/598)). The emission and absorption files hold every timestep of
  the model ([#625](https://github.com/artis-mcrt/artis/pull/625)).
- The wavelengths of `linestat.out` have nine significant digits in place of six
  ([#682](https://github.com/artis-mcrt/artis/pull/682)). *Tests: all.*
- `linelist.out` replaces `linestat.out`. It has one row for each line, a header row, and 0-based level indices
  ([#683](https://github.com/artis-mcrt/artis/pull/683)). The wavelengths keep the nine significant digits.
- `sn3d` replaces `deposition.out` in one step, and a failure of that step stops the run ([#534](https://github.com/artis-mcrt/artis/pull/534)).
- A file error during the write of the restart files stops the run. `sn3d` then keeps `input.txt` and the restart
  files of the timestep before ([#578](https://github.com/artis-mcrt/artis/pull/578)).
- With [WRITE_COMBINED_ALLRANK_OUT_FILES](options.md#write_combined_allrank_out_files), the estimator, nlte, and
  radfield files of all ranks go into one file each, e.g. `estimators_allranks.out`. Every preset sets the option
  ([#642](https://github.com/artis-mcrt/artis/pull/642), [#678](https://github.com/artis-mcrt/artis/pull/678),
  [#703](https://github.com/artis-mcrt/artis/pull/703)).
- The spectrum files have no space at the end of a line ([#387](https://github.com/artis-mcrt/artis/pull/387),
  [#599](https://github.com/artis-mcrt/artis/pull/599)).
- The packet files:
    - with `RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY` below 1, `nscatterings` does not count a line scattering
      ([#657](https://github.com/artis-mcrt/artis/pull/657));
    - `exspec` reads each column of a row that holds `nan`, and a row with too many columns stops the run
      ([#580](https://github.com/artis-mcrt/artis/pull/580));
    - they have no `stokes_i` column ([#386](https://github.com/artis-mcrt/artis/pull/386));
    - the absorption frequency also belongs to a free-free or a bound-free absorption
      ([#656](https://github.com/artis-mcrt/artis/pull/656)), and to a thermal redistribution
      ([#383](https://github.com/artis-mcrt/artis/pull/383));
    - a gamma packet has its emission position and time ([#637](https://github.com/artis-mcrt/artis/pull/637));
    - with [SAMPLE_RPKT_EMISSION](options.md#sample_rpkt_emission), each packet keeps one sampled r-packet
      emission, with its position and time ([#661](https://github.com/artis-mcrt/artis/pull/661),
      [#677](https://github.com/artis-mcrt/artis/pull/677)).
- The estimator files:
    - with [WRITE_ION_HEATING_COOLING_RATES](options.md#write_ion_heating_cooling_rates), each cell has the
      heating and cooling rates of each ion ([#581](https://github.com/artis-mcrt/artis/pull/581));
    - the `gammaestimator` line of an LTE or a grey cell has the unit of the other cells
      ([#630](https://github.com/artis-mcrt/artis/pull/630)).
- The nlte files hold no populations for an element without an NLTE solution
  ([#624](https://github.com/artis-mcrt/artis/pull/624)).
- `deposition.out`: the `Qdot_ana_erg/s/g` column includes the annihilation energy of the positron of a beta-plus
  decay ([#658](https://github.com/artis-mcrt/artis/pull/658)).
- `light_curve.out`: the comoving luminosity has no extra Lorentz factor
  ([#574](https://github.com/artis-mcrt/artis/pull/574)). *Tests: all.*
- The log of rank 0 holds the compile-time options of the run. After each timestep, it gives the progress and an
  estimate of the remaining time ([#700](https://github.com/artis-mcrt/artis/pull/700)).
- The restart file `gridsave_ts<n>.tmp` is a binary file ([#633](https://github.com/artis-mcrt/artis/pull/633)).

### Compile-time options

New options:

| Option | Purpose | Pull request |
| --- | --- | --- |
| [USE_MICROCLUMPING](options.md#use_microclumping), `MICROCLUMPING_FACTOR` | microclumping of the ejecta | [#275](https://github.com/artis-mcrt/artis/pull/275) |
| [NLTE_USE_GTH_SOLVER](options.md#nlte_use_gth_solver) | the Grassmann-Taksar-Heyman solver for the NLTE populations | [#490](https://github.com/artis-mcrt/artis/pull/490) |
| [ENABLE_CHARGE_TRANSFER_REACTIONS](options.md#enable_charge_transfer_reactions) | charge transfer reactions in the NLTE solver | [#535](https://github.com/artis-mcrt/artis/pull/535) |
| [NLTE_TIME_DEPENDENT_FIRST_TIMESTEP](options.md#nlte_time_dependent_first_timestep) | time-dependent ionisation and thermal balance | [#536](https://github.com/artis-mcrt/artis/pull/536) |
| [NLTE_TE_NNE_USE_ANDERSON_ACCEL](options.md#nlte_te_nne_use_anderson_accel), `NLTE_TE_NNE_RELTOL` | Anderson acceleration and tolerance of the iteration of the NLTE populations, T_e, and nne | [#539](https://github.com/artis-mcrt/artis/pull/539) |
| [WRITE_ION_HEATING_COOLING_RATES](options.md#write_ion_heating_cooling_rates) | heating and cooling rates of each ion in the estimator files | [#581](https://github.com/artis-mcrt/artis/pull/581) |
| [BINNED_OPACITY_LINE_WEIGHT](options.md#binned_opacity_line_weight), `LINE_BINNED_CAPPED_TAU_MAX` | the line weight of the binned opacities | [#621](https://github.com/artis-mcrt/artis/pull/621), [#705](https://github.com/artis-mcrt/artis/pull/705) |
| [WRITE_COMBINED_ALLRANK_OUT_FILES](options.md#write_combined_allrank_out_files) | one output file for all ranks | [#642](https://github.com/artis-mcrt/artis/pull/642), [#678](https://github.com/artis-mcrt/artis/pull/678) |
| [SAMPLE_RPKT_EMISSION](options.md#sample_rpkt_emission) | one sampled r-packet emission of each packet | [#661](https://github.com/artis-mcrt/artis/pull/661) |

[PARTICLE_THERMALISATION_SCHEME](options.md#particle_thermalisation_scheme) has the new values
`TIMEDEPENDENT_WITH_ADIABATIC_LOSS` ([#393](https://github.com/artis-mcrt/artis/pull/393)) and
`TIMEDEPENDENT_WITH_ADIABATIC_LOSS_ULTRARELATIVISTICLIMIT` ([#669](https://github.com/artis-mcrt/artis/pull/669)).

Options with a new name or a new form:

| Old option | New option | Pull request |
| --- | --- | --- |
| `MPKTS`, the packets of each rank | [NUM_PACKETS](options.md#num_packets), the packets of all ranks together | [#612](https://github.com/artis-mcrt/artis/pull/612) |
| `NT_ON`, `NT_SOLVE_SPENCERFANO`, `NT_EXCITATION_ON` | [NT_SCHEME](options.md#nt_scheme) | [#544](https://github.com/artis-mcrt/artis/pull/544) |
| `NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS` | `SF_RE_SOLVE_WITHIN_TIMESTEP_MIN_NNEPERION_FRACDIFF` | [#674](https://github.com/artis-mcrt/artis/pull/674) |
| `NLTEITER` | `NLTE_TE_NNE_MAXITER` | [#544](https://github.com/artis-mcrt/artis/pull/544) |
| `TABLESIZE` | `RATECOEFF_TABLESIZE` | [#544](https://github.com/artis-mcrt/artis/pull/544) |
| `DIRECT_COL_HEAT` | `COLLISIONAL_HEATING_FROM_LEVELPOPS` | [#544](https://github.com/artis-mcrt/artis/pull/544), [#705](https://github.com/artis-mcrt/artis/pull/705) |
| `MULTIBIN_RADFIELD_MODEL_ON`, `FIRST_NLTE_RADFIELD_TIMESTEP` | [MULTIBIN_RADFIELD_FIRST_TIMESTEP](options.md#multibin_radfield_first_timestep) | [#705](https://github.com/artis-mcrt/artis/pull/705) |
| `DETAILED_BF_ESTIMATORS_ON`, `DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP` | [DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP](options.md#detailed_bf_estimators_first_timestep) | [#705](https://github.com/artis-mcrt/artis/pull/705) |
| `LEVEL_HAS_BFEST` | `DETAILED_BF_ESTIMATORS_LEVEL_SELECTED` | [#705](https://github.com/artis-mcrt/artis/pull/705) |
| `NLTE_LIMIT_ION_STAGES_AFTER_FAILURE` | `NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION`, with no value for off | [#705](https://github.com/artis-mcrt/artis/pull/705) |
| `NLEVELS_REQUIRETRANSITIONS` | `NLEVELS_ADD_MISSING_COLLISIONAL_TRANSITIONS` | [#705](https://github.com/artis-mcrt/artis/pull/705) |
| `USE_XCOM_GAMMAPHOTOION` | `GAMMA_USE_XCOM_PHOTOION` | [#705](https://github.com/artis-mcrt/artis/pull/705) |
| `RPKT_USE_EXPANSION_OPACITIES`, `VPKT_USE_EXPANSION_OPACITIES` | `RPKT_USE_BINNED_OPACITIES`, `VPKT_USE_BINNED_OPACITIES` | [#705](https://github.com/artis-mcrt/artis/pull/705) |
| `FIXED_TIMESTEP_WIDTH`, `TIMESTEP_TRANSITION_TIME` | `TIMESTEP_SIZE_CONSTANT_MAX_WIDTH_DAYS`, `TIMESTEP_SIZE_TRANSITION_TIME_DAYS` | [#705](https://github.com/artis-mcrt/artis/pull/705) |

Removed options:

| Option | Fixed behaviour | Pull request |
| --- | --- | --- |
| `SFPTS`, `SF_EMIN`, `SF_EMAX`, `SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN` | constants in the non-thermal source file | [#385](https://github.com/artis-mcrt/artis/pull/385) |
| `NT_WORKFUNCTION_USE_SHELL_OCCUPANCY_FILE` | the code always reads the non-thermal data files of `data/` | [#380](https://github.com/artis-mcrt/artis/pull/380) |
| `TEMPERATURE_SOLVER_ACCURACY` | 1e-3 | [#544](https://github.com/artis-mcrt/artis/pull/544) |
| `KEEP_ALL_RESTART_FILES` | `sn3d` removes the restart files of the previous timestep | [#544](https://github.com/artis-mcrt/artis/pull/544) |
| `WRITE_EMISSIONABSORPTION_SPEC_AT_END` | `sn3d` always writes these files | [#544](https://github.com/artis-mcrt/artis/pull/544) |
| `USE_MODEL_INITIAL_ENERGY` | `INITIAL_PACKETS_ON` controls the initial energy of the model | [#544](https://github.com/artis-mcrt/artis/pull/544) |
| `SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS` | each timestep gets a new Spencer-Fano solution | [#674](https://github.com/artis-mcrt/artis/pull/674) |
| `SINGLE_GROUND_LEVEL` | no replacement | commit [fe33ca21](https://github.com/artis-mcrt/artis/commit/fe33ca21) |

The presets:

- `artisoptions_kilonova_nlte.h` is new ([#539](https://github.com/artis-mcrt/artis/pull/539),
  [#629](https://github.com/artis-mcrt/artis/pull/629)).
- `artisoptions_nltephotospheric.h` is the new name of `artisoptions_nltephotospheric_dynamic_ion_range.h`
  ([#629](https://github.com/artis-mcrt/artis/pull/629)).
- `artisoptions_nltewithoutnonthermal.h` and `artisoptions_christinenonthermal.h` are gone
  ([#544](https://github.com/artis-mcrt/artis/pull/544)).
- The build stops for some invalid combinations of options, e.g. the multibin radiation field model together with
  `USE_LUT_PHOTOION` ([#526](https://github.com/artis-mcrt/artis/pull/526),
  [#650](https://github.com/artis-mcrt/artis/pull/650), [#705](https://github.com/artis-mcrt/artis/pull/705)).

### Build, programs, and scripts

- The default branch has the name `main`. Its old name was `develop`.
- A build for an Nvidia GPU completes a run. Three corrections removed crashes in the list of the photoionisation
  cross-sections, in the state of each packet, and in the memory of the packet update ([#412](https://github.com/artis-mcrt/artis/pull/412), [#413](https://github.com/artis-mcrt/artis/pull/413),
  [#416](https://github.com/artis-mcrt/artis/pull/416)).
- `scripts/clean.sh` and the start of a new simulation also remove `out.txt.zst` ([#643](https://github.com/artis-mcrt/artis/pull/643)).
- A failed assertion or a fatal error stops all ranks at once. Before, the other ranks could wait until the end of
  the wall time ([#708](https://github.com/artis-mcrt/artis/pull/708)).
- The start of a run with a large model is faster. A model with 125000 cells needed 35 s in place of 151 s to read
  its input files ([#493](https://github.com/artis-mcrt/artis/pull/493)).
- The build needs no GSL library, and the make option `EIGEN=OFF` is gone
  ([#489](https://github.com/artis-mcrt/artis/pull/489)). Two small headers replace the copy of Boost.Math
  ([#505](https://github.com/artis-mcrt/artis/pull/505)).
- The Makefile finds libzstd, and `ZSTD=OFF` or `ZSTD=ON` controls it
  ([#632](https://github.com/artis-mcrt/artis/pull/632)). See [Make options](makeoptions.md).
- `PGO=GENERATE` and `PGO=USE` give a profile-guided optimisation
  ([#428](https://github.com/artis-mcrt/artis/pull/428)).
- `make` stops with a clear message when `artisoptions.h` is absent
  ([#526](https://github.com/artis-mcrt/artis/pull/526)), when `mpicxx` is absent, or when gcc is older than
  version 14 ([#698](https://github.com/artis-mcrt/artis/pull/698)).
- `make compile_commands.json` writes the compilation database
  ([#550](https://github.com/artis-mcrt/artis/pull/550)).
- `exspec` runs with more than one MPI rank ([#608](https://github.com/artis-mcrt/artis/pull/608)).
- `tests/runtest.sh` runs one end-to-end test ([#718](https://github.com/artis-mcrt/artis/pull/718)).
- The documentation site is new ([#720](https://github.com/artis-mcrt/artis/pull/720)).
- The job scripts:
    - `scripts/submit.sh` finds the cluster from the host name and submits the job script. The job scripts of JUWELS
      and MeluXina are gone ([#702](https://github.com/artis-mcrt/artis/pull/702)).
    - A job that completed only a part of the timesteps does not submit itself again without end
      ([#527](https://github.com/artis-mcrt/artis/pull/527)).
    - The job scripts queue the `exspec` job only when the simulation is complete
      ([#526](https://github.com/artis-mcrt/artis/pull/526)), and they build `exspec` together with `sn3d`
      ([#648](https://github.com/artis-mcrt/artis/pull/648)).
    - The job scripts log the core hours of each job ([#547](https://github.com/artis-mcrt/artis/pull/547)).
    - A job script stops when `SLURM_SUBMIT_DIR` has no value, e.g. when you run it by hand. Before, it changed to
      the home folder and worked there ([#604](https://github.com/artis-mcrt/artis/pull/604)).
    - `exspec-after.sh` starts one zstd process for each file, which makes the compression faster ([#587](https://github.com/artis-mcrt/artis/pull/587)).
    - The Leonardo cluster has job scripts (commit
      [3cf00789](https://github.com/artis-mcrt/artis/commit/3cf00789)).
    - `scripts/run-exspec-if-needed.sh` runs `exspec` only when `emission.out` is absent (commit
      [379bf44d](https://github.com/artis-mcrt/artis/commit/379bf44d)).
- `artistools filterpackets` and `artistools corehours` replace the scripts `filterpacketsescapedrpkt.py` and
  `sumcorehourslogs.py` ([#707](https://github.com/artis-mcrt/artis/pull/707)). `scripts/movefiles.sh` is gone,
  because `sn3d` writes the job folders itself ([#609](https://github.com/artis-mcrt/artis/pull/609)).
- `scripts/combine_allrank_out_files.py` makes the combined files of all ranks from the files of each rank
  ([#642](https://github.com/artis-mcrt/artis/pull/642)).

## v2026.5.15

Released on 2026-05-15. See the [release on GitHub](https://github.com/artis-mcrt/artis/releases/tag/v2026.5.15).

### The numerical results

- Four rate integrals have a new integration variable or a new form, so that they give no non-finite result.
  *Tests: classic, nebular, nltephotospheric.*
    - the radiative recombination coefficient ([#303](https://github.com/artis-mcrt/artis/pull/303));
    - the bound-free cooling coefficient ([#307](https://github.com/artis-mcrt/artis/pull/307));
    - the departure ratio ([#314](https://github.com/artis-mcrt/artis/pull/314));
    - the collisional recombination coefficient ([#315](https://github.com/artis-mcrt/artis/pull/315)).
- The packets of each cell run together, so that a gamma packet does not reset the cell cache
  ([#305](https://github.com/artis-mcrt/artis/pull/305)). *Tests: all.*
- The bound-free heating has no lookup table. The code integrates with the radiation temperature of the cell
  ([#316](https://github.com/artis-mcrt/artis/pull/316)). *Tests: classic.*
- The multibin radiation field model uses analytic series for the Planck integrals
  ([#320](https://github.com/artis-mcrt/artis/pull/320)). *Tests: nebular, nltephotospheric.*
- The interpolation of the XCOM photoelectric cross-section is correct
  ([#306](https://github.com/artis-mcrt/artis/pull/306)). *Tests: `kilonova_2d_xcomgammaphotoion`.*
- A build with `POL_ON = false` does no polarisation calculation
  ([#342](https://github.com/artis-mcrt/artis/pull/342)). *Tests: kilonova, nebular, nltephotospheric.*
- The kilonova preset has the spontaneous fission ([#343](https://github.com/artis-mcrt/artis/pull/343)). The code
  then enables the fission automatically, without an option ([#345](https://github.com/artis-mcrt/artis/pull/345)).
  *Tests: kilonova.*
- V48 has the correct branching of the beta-plus decay and the electron capture, and the correct positron energy
  ([#354](https://github.com/artis-mcrt/artis/pull/354)). *Tests: classic, kilonova.*
- The expansion opacity has no Doppler factor, as in the line-by-line mode
  ([#353](https://github.com/artis-mcrt/artis/pull/353)). *Tests: `kilonova_2d_expansionopac`.*
- Three changes of the arithmetic, with no change of the method
  ([#284](https://github.com/artis-mcrt/artis/pull/284), [#350](https://github.com/artis-mcrt/artis/pull/350),
  [#359](https://github.com/artis-mcrt/artis/pull/359)). *Tests: all.*

### Input and output files

- `sn3d` does not read or write `ratecoeff.dat`. It calculates the rate integrals at the start of each run
  ([#317](https://github.com/artis-mcrt/artis/pull/317)).
- `input.txt`: the lines of the grey opacity for the gamma rays, of the opacity choice, and of `rho_crit_para` have
  no effect. Compile-time options replace them ([#337](https://github.com/artis-mcrt/artis/pull/337)).
- `deposition.out`: the total deposition includes the spontaneous fission
  ([#292](https://github.com/artis-mcrt/artis/pull/292)).
- The light curve files hold no unused extra data ([#310](https://github.com/artis-mcrt/artis/pull/310)).
- The estimator files give the photoionisation rate coefficient `gamma_R` that the code uses, and
  `gamma_R_integral` with the detailed bound-free estimators
  ([#326](https://github.com/artis-mcrt/artis/pull/326)).
- The `vspecpol` and `vpkt_grid` files have no space at the end of a line
  ([#373](https://github.com/artis-mcrt/artis/pull/373)).

### Compile-time options

- New: `KEEP_ESCAPED_GAMMAS` ([#340](https://github.com/artis-mcrt/artis/pull/340)), `RPKT_GREY_TYPE` and
  `GAMMA_USE_KAPPA_GREY` ([#337](https://github.com/artis-mcrt/artis/pull/337)), and
  `VPKT_USE_EXPANSION_OPACITIES` ([#367](https://github.com/artis-mcrt/artis/pull/367)).
- `GRID_TYPE_OVERRIDE` replaces `GRID_TYPE`. With no value, the propagation grid has the dimensions of `model.txt`
  ([#334](https://github.com/artis-mcrt/artis/pull/334)).
- `GAMMA_THERMALISATION_SCHEME` and `PARTICLE_THERMALISATION_SCHEME` have separate lists of values
  ([#338](https://github.com/artis-mcrt/artis/pull/338)).
- `RPKT_USE_EXPANSION_OPACITIES` replaces `EXPANSIONOPACITIES_ON`
  ([#367](https://github.com/artis-mcrt/artis/pull/367)).
- `USE_ION_BFHEATING_ESTIMATORS` replaces `USE_LUT_BFHEATING`, and `RATECOEFF_INTEGRAL_ACCURACY` is gone
  ([#316](https://github.com/artis-mcrt/artis/pull/316)).
- `RADFIELDBINS_NU_MIN`, `RADFIELDBINS_NU_MAX`, and `RADFIELDBINS_T_E_SUPERBIN_NU_MAX` replace
  `nu_lower_first_initial`, `nu_upper_last_initial`, and `nu_upper_superbin`. `T_R_min` and `T_R_max` are gone
  ([#324](https://github.com/artis-mcrt/artis/pull/324)).
- Removed: `SEPARATE_STIMRECOMB` ([#301](https://github.com/artis-mcrt/artis/pull/301)), `RECORD_LINESTAT`
  ([#341](https://github.com/artis-mcrt/artis/pull/341)), and `DECAY_SPONTFISSION_ON`
  ([#345](https://github.com/artis-mcrt/artis/pull/345)).

### Build, programs, and scripts

- The integration and the root finder always use Boost.Math, or the Simpson rule in a GPU build. The GSL versions are
  gone ([#348](https://github.com/artis-mcrt/artis/pull/348)).
- The build needs gcc 14 or a later version ([#363](https://github.com/artis-mcrt/artis/pull/363)).
- The job scripts run `./artis/sn3d` directly, without a symlink in the model folder
  ([#347](https://github.com/artis-mcrt/artis/pull/347)).

## v2026.2.17

Released on 2026-02-17. See the [release on GitHub](https://github.com/artis-mcrt/artis/releases/tag/v2026.2.17).

### The numerical results

- Boost.Math does the integration and the root finding, and Eigen solves the matrix equations of the Spencer-Fano
  solver and of the NLTE populations. The repository holds both libraries, and they are the default
  ([#281](https://github.com/artis-mcrt/artis/pull/281), [#282](https://github.com/artis-mcrt/artis/pull/282),
  [#286](https://github.com/artis-mcrt/artis/pull/286), [#287](https://github.com/artis-mcrt/artis/pull/287)).
  `make BOOST=OFF EIGEN=OFF` selects the GSL code. *Tests: all.*
- The `nltenebular` preset removes unimportant ion stages from the NLTE solver and tests for population inversions
  ([#290](https://github.com/artis-mcrt/artis/pull/290)). *Tests: nebular.*

### Input files

- The code looks for an input file in `./`, in `./data/`, and in `./artis/data/`
  ([#276](https://github.com/artis-mcrt/artis/pull/276)).

### Compile-time options

- Removed: `IONGAMMA_POPFRAC_LEVELS_INCLUDED` ([#274](https://github.com/artis-mcrt/artis/pull/274)) and
  `NLTEPOP_FAILURE_USE_FIND_UPPERMOST_ION_FOR_LTE_RESET`, which is always on
  ([#291](https://github.com/artis-mcrt/artis/pull/291)).
- `RATECOEFF_INTEGRAL_ACCURACY` also has the function of `CONTINUUM_NU_INTEGRAL_ACCURACY`, which is gone (commit
  [48149029](https://github.com/artis-mcrt/artis/commit/48149029)).

## v2026.01.08

Released on 2026-01-08. See the [release on GitHub](https://github.com/artis-mcrt/artis/releases/tag/v2026.01.08).

### The numerical results

- `data/` holds gamma-ray tables for the alpha and the beta decays
  ([#273](https://github.com/artis-mcrt/artis/pull/273)). *Tests: kilonova.*

## v2026.01.05

Released on 2026-01-05. See the [release on GitHub](https://github.com/artis-mcrt/artis/releases/tag/v2026.01.05).

### The numerical results

- `NLEVELS_REQUIRETRANSITIONS` adds the missing transitions correctly
  ([#270](https://github.com/artis-mcrt/artis/pull/270)). *Tests: kilonova, nebular, nltephotospheric.*
- A packet with a position error at a boundary of a 3D grid moves to the correct cell again, which removes the
  boundary errors ([#271](https://github.com/artis-mcrt/artis/pull/271)).
- A 1D spherical grid handles the empty cells correctly.

## v2025.11.11

Released on 2025-11-11. See the [release on GitHub](https://github.com/artis-mcrt/artis/releases/tag/v2025.11.11).

### The numerical results

- The abundances of the stable nuclei that have a column in `model.txt` are correct
  ([#238](https://github.com/artis-mcrt/artis/pull/238)). *Tests: kilonova.*
- A standard nuclide that the input composition does not give, e.g. Ni57 or Co57, gets an abundance of zero. Before,
  the value was not initialised (commit [03d1131e](https://github.com/artis-mcrt/artis/commit/03d1131e)).
- The expansion opacity subtracts the stimulated emission, as the line-by-line mode does
  ([#232](https://github.com/artis-mcrt/artis/pull/232)). *Tests: `kilonova_2d_expansionopac`.*
- A pellet selects the frequency of its gamma ray when the code makes the packet
  ([#246](https://github.com/artis-mcrt/artis/pull/246)). *Tests: all.*
- The nuclide list holds all daughter nuclei ([#257](https://github.com/artis-mcrt/artis/pull/257)), and the decay
  solver has a new structure ([#258](https://github.com/artis-mcrt/artis/pull/258)). *Tests: all.*
- The decay network has the spontaneous fission of Cf254, with the option `DECAY_SPONTFISSION_ON`
  ([#261](https://github.com/artis-mcrt/artis/pull/261)).

### Output files

- The packet files hold the position of the true emission (`trueem_pos`) in place of the true emission velocity
  ([#227](https://github.com/artis-mcrt/artis/pull/227)).
- The estimator files give the isotopic abundances in the order of the mass number
  ([#239](https://github.com/artis-mcrt/artis/pull/239)).

### Compile-time options

- `RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY` is a `std::optional<float>`
  ([#226](https://github.com/artis-mcrt/artis/pull/226)).
