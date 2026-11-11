# Methods and references

This page gives each method of ARTIS: a short description, the options that select it, and its reference. Use it when
you describe the method of your results. [How a simulation works](simulation.md) gives the sequence of the steps, and
[Compile-time options](options.md) describes each option in full.

The log file of {term}`rank` 0 holds the git version and the options of a run (see [Running sn3d](running.md)).
`artistools describe <model folder>` gives the main options of a run together with their effect.

## The energy packets

Code: `packet.cc` and `decay.cc` (see [Code map](codemap.md)).

The radiation field is a set of indivisible energy packets. A packet changes its type when its energy changes its
form, and the code follows it until it escapes or the simulation ends. [How a simulation works](simulation.md) gives
the types of the packets.

The sum of the energy of the {term}`pellets <pellet>` of a rank is the energy that the model releases in the simulation:

$$
\varepsilon = \frac{1}{N_\mathrm{rank}} \sum_\mathrm{cells} V(t_\mathrm{min}) \, \rho(t_\mathrm{min}) \,
\left( e_\mathrm{decay} + q \right)
$$

- $\varepsilon$ is the energy of one pellet in the comoving frame.
- $N_\mathrm{rank}$ is the number of packets of the rank, so each rank holds the energy of the whole model.
- $V$ and $\rho$ are the volume and the density of a {term}`propagation cell`.
- $e_\mathrm{decay}$ is the decay energy per unit mass of the cell in the time range of the simulation.
- $q$ is the `q` column of `model.txt`, multiplied by $t_\mathrm{model} / t_\mathrm{min}$ for the adiabatic loss. It is
  zero without [INITIAL_PACKETS_ON](options.md#initial_packets_on).

The code selects a cell in proportion to its energy, and then a decay chain. With
[UNIFORM_PELLET_ENERGIES](options.md#uniform_pellet_energies), it adds one exponential sample of the mean lifetime for
each nuclide of the chain to the snapshot time of `model.txt`. It samples again until the result is inside the time
range of the simulation.

- Options:
    - [NUM_PACKETS](options.md#num_packets) sets the number of packets of all ranks together.
    - [UNIFORM_PELLET_ENERGIES](options.md#uniform_pellet_energies) gives each pellet the same energy. With the value
      false, the decay times are uniform instead, and the energy of each pellet follows the decay power at its time.
    - [INITIAL_PACKETS_ON](options.md#initial_packets_on) adds the initial energy of the cells to the energy that the
      packets share. The number of packets stays the same.
- References: Lucy (2002), A&A, 384, 725-735, doi:10.1051/0004-6361:20011756 describes the packets. Sim (2007),
  MNRAS, 375, 154-162, doi:10.1111/j.1365-2966.2006.11271.x and Kromer & Sim (2009), MNRAS, 398, 1809-1826,
  doi:10.1111/j.1365-2966.2009.15256.x describe an early version of the code.

## The opacity of the r-packets

Code: `rpkt.cc` (see [Code map](codemap.md)).

An {term}`r-packet` moves until it reaches a line, a continuum absorption, or an electron scattering. The default mode gives
each line its own {term}`Sobolev <Sobolev approximation>` optical depth. The expansion of the ejecta shifts the comoving frequency of the packet
through the lines one after the other. The continuum opacity is the sum of the electron scattering, the
bound-free, and the free-free opacities.

The Sobolev optical depth of a line from the lower level $l$ to the upper level $u$ is:

$$
\tau_{lu} = \max\left[ 0, \; \left( B_{lu} n_l - B_{ul} n_u \right) \frac{hc}{4\pi} \, t \right]
$$

$B_{lu}$ and $B_{ul}$ are the Einstein coefficients, $n_l$ and $n_u$ are the level populations, and $t$ is the time at
which the packet reaches the line. A packet samples an optical depth $\tau_\mathrm{rnd} = -\ln \xi$ with $\xi$ a
uniform random number in $(0, 1)$. The event occurs where the sum of the line optical depths and of the continuum
optical depth $\chi_\mathrm{cont} \, d$ reaches $\tau_\mathrm{rnd}$. The distance $d$ is the path of the packet.

- Options:
    - [USE_RELATIVISTIC_DOPPLER_SHIFT](options.md#use_relativistic_doppler_shift) selects the full Doppler factor.
    - [PHIXS_CLASSIC_NO_INTERPOLATION](options.md#phixs_classic_no_interpolation) takes the photoionisation cross
      section from the table point at or below the frequency. Without it, the code interpolates between two table
      points.
    - [DIPOLE](options.md#dipole) selects the direction of an electron scattering.
    - [NU_MIN_R](options.md#nu_min_r) and [NU_MAX_R](options.md#nu_max_r) set the frequency limits of the spectra.
      They also limit the frequency that a {term}`k-packet` samples from a Planck function, so they change the transport.
- The alternatives are the [binned opacities](#the-binned-opacities-and-the-thermalisation-probability) and the
  [grey opacity of a thick cell](#the-grey-opacity-of-a-thick-cell).

## The macro atom and the k-packets

Code: `macroatom.cc` and `kpkt.cc` (see [Code map](codemap.md)).

A {term}`macro atom` makes internal radiative and collisional transitions between the levels and the ions. It stops when it
emits an r-packet or becomes a k-packet. Each type of absorption has its own result:

- In the default mode, an absorption in a line activates a macro atom.
- A bound-free absorption activates a macro atom with the probability $\min(1, \nu_\mathrm{edge} / \nu)$, where
  $\nu_\mathrm{edge}$ is the threshold frequency of the continuum. It makes a k-packet otherwise.
- A free-free absorption makes a k-packet.

A k-packet holds thermal energy. The code adds the
cooling rate of each process of the cell, and it samples one of them. That process converts the k-packet into
radiation or into a macro atom. The cooling processes are the free-free emission, the free-bound emission, the
collisional excitation, and the collisional ionisation.

- Options: [RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY](options.md#rpkt_boundbound_thermalisation_probability)
  replaces the macro atom for an absorption in a line when it has a value (see
  [the binned opacities](#the-binned-opacities-and-the-thermalisation-probability)). A bound-free absorption then
  still activates a macro atom. A k-packet does not: with a value of this option, each k-packet emits an r-packet
  directly from the thermal emission of the cell. The thermal balance has no term for this line thermalisation, so
  a nonzero value needs LTE timesteps only (`num_lte_timesteps` of `input.txt`).
- References: Lucy (2002), A&A, 384, 725-735, doi:10.1051/0004-6361:20011756 and Lucy (2003), A&A, 403, 261-275,
  doi:10.1051/0004-6361:20030357.

## The gas state

The grid update solves the state of each {term}`model grid cell` that contains matter: the electron temperature T_e, the
electron density, and the ion and level populations. Kromer & Sim (2009), MNRAS, 398, 1809-1826,
doi:10.1111/j.1365-2966.2009.15256.x describes the classic mode. Shingles et al. (2020), MNRAS, 492, 2029-2043,
doi:10.1093/mnras/stz3412 describes the mode with the non-LTE ({term}`NLTE`) populations. [Option presets](presets.md) says
which {term}`preset` selects which mode.

### The populations without NLTE levels

Code: `ltepop.cc` (see [Code map](codemap.md)).

The first `num_lte_timesteps` {term}`timesteps <timestep>` of `input.txt` solve the ionisation balance in local thermodynamic equilibrium
({term}`LTE`). After these timesteps, an element without NLTE levels gets a photoionisation equilibrium for its ions and
Boltzmann populations for its levels. An element that [FORCE_SAHA_ION_BALANCE](options.md#force_saha_ion_balance)
selects keeps the Saha ionisation balance. The kilonova preset selects each element except Sr, Y, and Zr.

- Options:
    - [LTEPOP_EXCITATION_USE_TJ](options.md#ltepop_excitation_use_tj) selects the temperature of the Boltzmann
      factor.
    - [USE_LUT_PHOTOION](options.md#use_lut_photoion) selects the tabulated photoionisation coefficients.
    - [FORCE_SAHA_ION_BALANCE](options.md#force_saha_ion_balance) selects the Saha ionisation balance of an element.
    - [SINGLE_LEVEL_TOP_ION](options.md#single_level_top_ion) keeps one level for the highest ion stage of each
      element.
- Reference: the Boltzmann and the Saha relations are standard, e.g. Mihalas (1978), Stellar Atmospheres, 2nd ed.,
  W. H. Freeman, San Francisco.

### The NLTE populations

Code: `nltepop.cc` (see [Code map](codemap.md)).

The solver puts the radiative rates, the collisional rates, and the non-thermal ionisation and excitation rates of
all the levels of all the ions of an element into one rate matrix. The solution of the statistical equilibrium gives the level
populations $n_i$:

$$
\sum_{j \neq i} n_j R_{ji} - n_i \sum_{j \neq i} R_{ij} = 0
$$

$R_{ij}$ is the total rate for one particle from level $i$ to level $j$. With the LU decomposition, the condition
$\sum_i n_i = n_\mathrm{element}$ replaces the first equation of the set. The GTH state elimination
([NLTE_USE_GTH_SOLVER](options.md#nlte_use_gth_solver)) solves the full set and then scales the populations to this sum.

Each ion has a ground level and a number of excited levels in full NLTE. An ion with more levels above
them also has one {term}`superlevel` for those levels. An autoionising level is not in the superlevel, and it keeps its own
row in the matrix. The cell repeats the solutions of the populations, of T_e, and of the
non-thermal electrons until T_e and the electron density change by less than a tolerance. With
[NLTE_TE_NNE_USE_ANDERSON_ACCEL](options.md#nlte_te_nne_use_anderson_accel) or
[ENABLE_CHARGE_TRANSFER_REACTIONS](options.md#enable_charge_transfer_reactions), the iteration also tests the change of
the populations of the significant ions. With the acceleration, it also needs an estimate of the remaining error below
the tolerance. The iteration also stops
after [NLTE_TE_NNE_MAXITER](options.md#nlte_te_nne_maxiter) + 1 passes. The cell then keeps the last solution, and the
log of the rank has a message with the text `failed to converge`.

From the timestep of [NLTE_TIME_DEPENDENT_FIRST_TIMESTEP](options.md#nlte_time_dependent_first_timestep), the
equations of the ion populations have a time term. The ion populations are then not in a statistical equilibrium,
but the level populations inside each ion are. The equation of the ground level of each ion, except the lowest ion,
becomes:

$$
\sum_{j \neq g} n_j R_{jg} - n_g \sum_{j \neq g} R_{gj} - \frac{n_\mathrm{ion}}{\Delta t}
= -\frac{n_\mathrm{ion}^\mathrm{old}}{\Delta t}
$$

$n_\mathrm{ion}$ is the sum of the populations of the ion, and $n_\mathrm{ion}^\mathrm{old}$ is its value at the
grid update before. $\Delta t$ is the time between the two grid updates. The sum over the levels of an ion removes the transitions inside the ion, so the sum is the
backward Euler step of the ion population. The
time term needs a complete NLTE solution of the cell from the timestep before. A cell without one keeps the
equations with no time term for that timestep, e.g. after an LTE timestep, after the {term}`grey opacity`, or after a
fallback to LTE.

If the solution of an element fails, the solver removes one ion and solves again, or the element falls back to
LTE. The log of the rank has a warning for each fallback to LTE.

- Options:
    - [ION_NLEVELS_EXCITED_NLTE](options.md#ion_nlevels_excited_nlte) sets the number of the excited levels in full
      NLTE.
    - [NLEVELS_ADD_MISSING_COLLISIONAL_TRANSITIONS](options.md#nlevels_add_missing_collisional_transitions) adds
      collisional transitions where the data has none.
    - [NLTE_USE_GTH_SOLVER](options.md#nlte_use_gth_solver) selects the Grassmann-Taksar-Heyman state elimination
      in place of the lower-upper (LU) decomposition. An element that
      [FORCE_SAHA_ION_BALANCE](options.md#force_saha_ion_balance) selects, and a solution with the time term, always
      use the LU decomposition.
    - [NLTE_TE_NNE_MAXITER](options.md#nlte_te_nne_maxiter), [NLTE_TE_NNE_RELTOL](options.md#nlte_te_nne_reltol),
      and [NLTE_TE_NNE_USE_ANDERSON_ACCEL](options.md#nlte_te_nne_use_anderson_accel) control the iteration of
      the cell.
    - [STRICT_POPULATION_CHECKING](options.md#strict_population_checking) and
      [NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION](options.md#nlte_limit_ion_stages_max_levelpop_over_elementpop_remove_ion)
      control the test of a solution and the removal of an ion.
    - [NLTE_TIME_DEPENDENT_FIRST_TIMESTEP](options.md#nlte_time_dependent_first_timestep) adds a time term to the
      ionisation balance and to the thermal balance.
    - [ENABLE_CHARGE_TRANSFER_REACTIONS](options.md#enable_charge_transfer_reactions) adds the charge transfer
      reactions.
- References:
    - the solver: Shingles et al. (2020), MNRAS, 492, 2029-2043, section 2.3, doi:10.1093/mnras/stz3412;
    - the state elimination: Grassmann, Taksar & Heyman (1985), Operations Research, 33, 1107-1116,
      doi:10.1287/opre.33.5.1107;
    - the acceleration of the iteration: Anderson (1965), Journal of the ACM, 12, 547-560,
      doi:10.1145/321296.321305;
    - the time term: Pognan, Jerkstrand & Grumer (2022), MNRAS, 510, 3806-3837, equations 8 and 17,
      doi:10.1093/mnras/stab3674;
    - the charge transfer reactions: Kingdon & Ferland (1996), ApJS, 106, 205-211, doi:10.1086/192335 and Butler &
      Dalgarno (1980), ApJ, 241, 838-843, doi:10.1086/158395.

### The thermal balance

Code: `thermalbalance.cc` (see [Code map](codemap.md)).

A cell in an LTE timestep, and a cell with the grey opacity, have no thermal balance. The code sets their T_e to the
temperature of the mean intensity. The first timestep of a new simulation is different: each cell keeps its initial
temperature.

The code calculates the initial temperature from the energy of the decays between the snapshot time of `model.txt`
and the midpoint of the first timestep. With [INITIAL_PACKETS_ON](options.md#initial_packets_on), it adds the energy
of the `q` column of `model.txt`. It takes all of that energy as radiation that stays in the cell. That radiation
loses energy to the expansion in proportion to 1/t. A decay at the time `t_decay` thus keeps the fraction
`t_decay / t_mid` of its energy at the midpoint `t_mid` of the first timestep. The `q` energy keeps the fraction
`t_model / t_mid`, where `t_model` is the snapshot time.

For each other cell, the code solves this equation for the electron temperature. One grid update changes T_e by a
maximum factor of two:

$$
H_\mathrm{bf} + H_\mathrm{ff} + H_\mathrm{coll} + H_\mathrm{dep} =
C_\mathrm{ff} + C_\mathrm{fb} + C_\mathrm{coll} + C_\mathrm{ad} + C_\mathrm{time}
$$

The heating rates $H$ are the bound-free heating, the free-free heating, the collisional de-excitation, and the
non-thermal heating. The non-thermal heating is:

$$
H_\mathrm{dep} = f_\mathrm{heat} \, D_\mathrm{lepton} + D_\mathrm{alpha} + D_\mathrm{fission}
$$

$D$ is a deposition rate density, and $f_\mathrm{heat}$ is the heating fraction of the non-thermal solution. It
includes the excitation energy of the transitions outside the stored list
([MAX_NT_EXCITATIONS_STORED](options.md#max_nt_excitations_stored)), which the packets convert to k-packets. With
[COLLISIONAL_HEATING_FROM_LEVELPOPS](options.md#collisional_heating_from_levelpops), the stored list holds no
transition of an element without NLTE levels, so the excitation energy of such an element is also heating. With
[USE_ION_BFHEATING_ESTIMATORS](options.md#use_ion_bfheating_estimators), a Monte Carlo {term}`estimator` corrects the
bound-free heating (see [the radiation field model](#the-radiation-field-model)). Without
[COLLISIONAL_HEATING_FROM_LEVELPOPS](options.md#collisional_heating_from_levelpops), the collisional heating comes
from a Monte Carlo estimator, which also holds the collisional recombination heating.

The cooling rates $C$ are the free-free, the free-bound, and the collisional cooling, and the adiabatic cooling of
the expansion. With the gas pressure $p = (n_\mathrm{ion} + n_e) \, k_\mathrm{B} T_e$ and a volume that grows as
$t^3$, the adiabatic cooling is $C_\mathrm{ad} = 3 p / t$.

From the timestep of [NLTE_TIME_DEPENDENT_FIRST_TIMESTEP](options.md#nlte_time_dependent_first_timestep),
the cooling also holds the change of the thermal energy of the gas in the timestep:

$$
C_\mathrm{time} = \frac{3}{2} k_\mathrm{B} \,
\frac{n_\mathrm{tot} \left( T_e - T_e^\mathrm{old} \right) + T_e^\mathrm{old} \left( n_e - n_e^\mathrm{old} \right)}{\Delta t}
$$

$n_\mathrm{tot}$ is $n_\mathrm{ion} + n_e$. This term has the same condition as the time
term of the ion populations: the cell needs a complete NLTE solution from the timestep before.

The code keeps T_e between [MINTEMP](options.md#mintemp) and [MAXTEMP](options.md#maxtemp). A cell at one of these
limits has no thermal balance, and the log of its rank has a warning (see
[Troubleshooting](troubleshooting.md#the-log)).

- Options: [COLLISIONAL_HEATING_FROM_LEVELPOPS](options.md#collisional_heating_from_levelpops),
  [BFCOOLING_USELEVELPOPNOTIONPOP](options.md#bfcooling_uselevelpopnotionpop),
  [USE_ION_BFHEATING_ESTIMATORS](options.md#use_ion_bfheating_estimators), and
  [WRITE_ION_HEATING_COOLING_RATES](options.md#write_ion_heating_cooling_rates).

## The radiation field model

Code: `radfield.cc` (see [Code map](codemap.md)).

The packets add to the estimators of the radiation field of each cell that they enter. The mean intensity of a cell
is:

$$
J = \frac{1}{4\pi \, V \, \Delta t \, N_\mathrm{ranks}} \sum_\mathrm{paths} \varepsilon_\mathrm{cmf} \, ds
$$

$V$ is the volume of the model grid cell at the midpoint of the timestep of the packets, and $\Delta t$ is the width
of that timestep. $\varepsilon_\mathrm{cmf}$ is the comoving energy of a packet, and $ds$ is the length of a path in
the cell. The sum holds the paths of all ranks. Each rank holds the energy of the whole model, so the code divides by
$N_\mathrm{ranks}$. The estimator of $\nu J$ has the additional factor $\nu_\mathrm{cmf}$ in the sum.

The grid update fits a dilute
blackbody to the full spectrum of each cell, with a dilution factor W and a radiation temperature T_R. With the
multibin model, it also fits a dilute blackbody to each frequency bin. The last bin is an exception: its temperature
is the electron temperature of the cell, and only its dilution factor comes from the fit.

The rates of the photoionisation and of the bound-free heating use these fits. The photoionisation rate of a
continuum with a detailed bound-free estimator is an exception (see the next section). The bound-free heating of
that continuum still uses the fits.

[USE_ION_BFHEATING_ESTIMATORS](options.md#use_ion_bfheating_estimators) is a second exception, and the classic preset
sets it. The packets then add to a Monte Carlo estimator of the bound-free heating of each ground continuum. The
code divides that estimator by the bound-free heating of the same continuum from the fits. It multiplies the
bound-free heating of each level from the fits by this ratio of the ground continuum with the nearest threshold. If
the heating from the fits is zero, the ratio becomes 1 and the log has a warning.

A cell in an LTE timestep, and a cell with the grey opacity, get no fit. The code then sets T_R to the temperature
of the mean intensity, and W to 1. In the first timestep of a new simulation, T_R is the initial temperature of the
cell (see [the thermal balance](#the-thermal-balance)).

- Options:
    - [MULTIBIN_RADFIELD_FIRST_TIMESTEP](options.md#multibin_radfield_first_timestep) sets the first timestep that
      reads the fits of the bins.
    - [RADFIELDBINCOUNT](options.md#radfieldbincount), [RADFIELDBINS_NU_MIN](options.md#radfieldbins_nu_min),
      [RADFIELDBINS_NU_MAX](options.md#radfieldbins_nu_max), and
      [RADFIELDBINS_T_E_SUPERBIN_NU_MAX](options.md#radfieldbins_t_e_superbin_nu_max) set the bins.
    - [DETAILED_LINE_ESTIMATORS_ON](options.md#detailed_line_estimators_on) stores the estimators of selected lines.
- Reference: Shingles et al. (2020), MNRAS, 492, 2029-2043, section 2.2, doi:10.1093/mnras/stz3412.

## The photoionisation rates and the bound-free estimators

Code: `ratecoeff.cc` and `radfield.cc` (see [Code map](codemap.md)).

The photoionisation rate coefficient of a level comes from one of three sources:

- a table for a Planck radiation field, which the code multiplies by the dilution factor W and by the ratio of a Monte
  Carlo estimator to the blackbody rate ([USE_LUT_PHOTOION](options.md#use_lut_photoion));
- an integral of the cross section over the radiation field model of the cell;
- a detailed bound-free estimator, which the packets build for each selected continuum
  ([DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP](options.md#detailed_bf_estimators_first_timestep) and
  [DETAILED_BF_ESTIMATORS_LEVEL_SELECTED](options.md#detailed_bf_estimators_level_selected)).

The detailed estimator of a continuum is:

$$
\gamma = \frac{1}{h \, V \, \Delta t \, N_\mathrm{ranks}}
\sum_\mathrm{paths} \sigma_\mathrm{contr}(\nu_\mathrm{cmf}) \, \frac{\varepsilon_\mathrm{cmf}}{\nu_\mathrm{cmf}} \, ds
$$

$\sigma_\mathrm{contr}$ is the cross section of the continuum, multiplied by the probability of its target level and
by a factor for the stimulated recombination. The volume, the width, and $\varepsilon_\mathrm{cmf}$ are the same as
for the mean intensity.

If the {term}`model folder` holds `recombrates.txt`, the code calibrates the total recombination rate of an ion to that file
at the temperature [RECOMBCALIBRATION_T_ELEC](options.md#recombcalibration_t_elec).
[Input files](inputfiles.md#recombratestxt) gives the rules of the calibration, and the log gives the result for each
ion.

- Reference: Shingles et al. (2020), MNRAS, 492, 2029-2043, section 2.2, doi:10.1093/mnras/stz3412.

## The non-thermal electrons

Code: `nonthermal.cc` (see [Code map](codemap.md)).

The deposited energy of the gamma rays, of the positrons, and of the electrons makes fast leptons. These leptons
lose their energy to the heating, the ionisation, and the excitation of the gas. The
Spencer-Fano solver finds the energy spectrum $y(E)$ of these electrons on a uniform energy grid. The code solves the
integral form of the equation of Kozma & Fransson (1992), ApJ, 390, 602-621, doi:10.1086/171311, equation 7:

$$
\begin{aligned}
& L(E) \, y(E)
+ \sum_{\mathrm{exc}} n_l \int_E^{E + \epsilon_\mathrm{exc}} y(E') \, \sigma_\mathrm{exc}(E') \, dE' \\
& + \sum_{\mathrm{shells}} n_\mathrm{ion} \int_E^{E_\mathrm{max}} y(E')
\int_{\max(E' - E, \, I)}^{(E' + I)/2} \sigma_\mathrm{ion}(E', \epsilon) \, d\epsilon \, dE' \\
& - \sum_{\mathrm{shells}} n_\mathrm{ion} \int_{2E + I}^{E_\mathrm{max}} y(E')
\int_{E + I}^{(E' + I)/2} \sigma_\mathrm{ion}(E', \epsilon) \, d\epsilon \, dE'
= \int_E^{E_\mathrm{max}} S(E') \, dE'
\end{aligned}
$$

- $L(E)$ is the energy loss rate of a fast electron to the free thermal electrons.
- $\epsilon_\mathrm{exc}$ and $\sigma_\mathrm{exc}$ are the energy and the cross section of an excitation from a
  level with the population $n_l$.
- $I$ is the ionisation potential of a shell, and $\epsilon$ is the energy that an ionisation transfers.
- $\sigma_\mathrm{ion}(E', \epsilon)$ is the differential ionisation cross section.
- $S(E)$ is the source of electrons. It is constant over a narrow range at the top of the energy grid.

The first ionisation integral holds the primary electrons that fall below $E$. The second integral holds the secondary
electrons that start above $E$. The code writes the equation as an upper triangular matrix and solves it by
back-substitution. The spectrum gives the
three fractions of the deposited energy, and the non-thermal ionisation and excitation rates for the NLTE solver.
With [SF_AUGER_CONTRIBUTION_ON](options.md#sf_auger_contribution_on), the equation also has the source term of the
Auger electrons. The work function approximation is the alternative. It gives ionisation rates and no excitation
rates, and the thermal balance then counts all of the deposited lepton energy as heat.

- Options:
    - [NT_SCHEME](options.md#nt_scheme) selects the {term}`Spencer-Fano solution`, the work function approximation, or no
      non-thermal ionisation. With no non-thermal ionisation, the run also has no non-thermal excitation, and all of
      the deposited lepton energy heats the gas.
    - [SF_RE_SOLVE_WITHIN_TIMESTEP_MIN_NNEPERION_FRACDIFF](options.md#sf_re_solve_within_timestep_min_nneperion_fracdiff)
      controls a new solution inside a timestep.
    - [NTEXCITATION_MAXNLEVELS_LOWER](options.md#ntexcitation_maxnlevels_lower),
      [NTEXCITATION_MAXNLEVELS_UPPER](options.md#ntexcitation_maxnlevels_upper), and
      [MAX_NT_EXCITATIONS_STORED](options.md#max_nt_excitations_stored) limit the excitation transitions.
    - [NT_USE_VALENCE_IONPOTENTIAL](options.md#nt_use_valence_ionpotential),
      [NT_MAX_AUGER_ELECTRONS](options.md#nt_max_auger_electrons), and
      [SF_AUGER_CONTRIBUTION_ON](options.md#sf_auger_contribution_on) control the ionisation potentials and the
      Auger electrons.
- References:
    - the equation: Spencer & Fano (1954), Phys. Rev., 93, 1172-1181, doi:10.1103/PhysRev.93.1172;
    - the application to supernovae: Kozma & Fransson (1992), ApJ, 390, 602-621, doi:10.1086/171311;
    - the solver of ARTIS and the Auger source term: Shingles et al. (2020), MNRAS, 492, 2029-2043, section 2.5,
      doi:10.1093/mnras/stz3412;
    - the work function approximation: Axelrod (1980), PhD thesis, University of California, Santa Cruz.

## The decays and the thermalisation of the particles

Code: `decay.cc` and `update_packets.cc` (see [Code map](codemap.md)).

The decay network holds the alpha decays, the beta-minus decays, the beta-plus decays, the electron captures, and
the spontaneous fission decays. The code gets the abundance of each nuclide at
each time from the solutions of the Bateman equation. The `data` folder of the repository holds the decay data and
the gamma-ray line tables of some nuclides. [Input files](inputfiles.md#the-gamma-ray-line-files) says which lines a
nuclide without a table gets. The decay particles deposit their energy in one of three ways:

- They deposit all of their energy at once.
- The code transports them, and they lose their energy over time.
- With an analytic {term}`thermalisation` efficiency, a packet of particles deposits its energy at once with the probability
  of that efficiency. If not, the packet escapes and deposits nothing.

The fragments of a spontaneous fission are an exception. They always deposit all of their energy at once.

The deposited energy of the alpha particles and of the fission fragments is heat only. The packet becomes a
k-packet, and its energy does not go into the Spencer-Fano solution. The deposited energy of the electrons and of
the positrons goes to the non-thermal electrons (see [the non-thermal electrons](#the-non-thermal-electrons)).

- Options: [PARTICLE_THERMALISATION_SCHEME](options.md#particle_thermalisation_scheme) selects the method:
    - `INSTANTFULLDEPOSITION` deposits the energy at once.
    - `TIMEDEPENDENT` transports the particles with the Monte Carlo method.
    - `TIMEDEPENDENT_WITH_ADIABATIC_LOSS` adds the adiabatic loss of the kinetic energy to the collisional loss. Only
      the collisional share of the lost energy heats the gas.
    - `TIMEDEPENDENT_WITH_ADIABATIC_LOSS_ULTRARELATIVISTICLIMIT` uses the rate E/t for the adiabatic loss. The
      kilonova preset selects this value.
    - `TIMEDEPENDENTWITHGAMMAPRODUCTS` also transports the electrons and the positrons that the interactions of the
      gamma rays make.
    - `BARNES` and `WOLLAEGER` use an analytic thermalisation efficiency.
- References:
    - the decays and the time-dependent thermalisation: Shingles et al. (2023), ApJL, 954, L41,
      doi:10.3847/2041-8213/acf29a;
    - the analytic efficiencies: Barnes, Kasen, Wu & Martínez-Pinedo (2016), ApJ, 829, 110,
      doi:10.3847/0004-637X/829/2/110 and Wollaeger, Korobkin, Fontes, Rosswog, Even & Fryer (2018), MNRAS, 478,
      3298-3334, doi:10.1093/mnras/sty1018.

## The gamma rays

Code: `gammapkt.cc` (see [Code map](codemap.md)).

A {term}`gamma packet` gives its energy to the gas by Compton scattering, by photoelectric absorption, or by pair production.
Every preset transports the gamma packets with frequency-dependent opacities. The alternatives are a grey opacity
and an analytic thermalisation efficiency. With that efficiency, a gamma packet deposits all of its energy at once
with the probability of the efficiency, and it escapes if not.

- Options:
    - [GAMMA_THERMALISATION_SCHEME](options.md#gamma_thermalisation_scheme) selects the transport or an analytic
      efficiency.
    - [GAMMA_USE_XCOM_PHOTOION](options.md#gamma_use_xcom_photoion) selects the XCOM table for the photoelectric
      absorption.
    - [GAMMA_USE_KAPPA_GREY](options.md#gamma_use_kappa_grey) selects a grey opacity.
    - [KEEP_ESCAPED_GAMMAS](options.md#keep_escaped_gammas) writes the gamma-ray light curve and spectrum.
- References:
    - the photoelectric power laws: Ambwani & Sutherland (1988), ApJ, 325, 820-827, eq. 2, doi:10.1086/166052, after
      Veigele (1973), Atomic Data and Nuclear Data Tables, 5, 51-111, doi:10.1016/S0092-640X(73)80015-4;
    - the analytic efficiencies: Barnes, Kasen, Wu & Martínez-Pinedo (2016), ApJ, 829, 110,
      doi:10.3847/0004-637X/829/2/110, Wollaeger, Korobkin, Fontes, Rosswog, Even & Fryer (2018), MNRAS, 478,
      3298-3334, doi:10.1093/mnras/sty1018, and Guttman, Shenhar, Sarkar & Waxman (2024), MNRAS, 533, 994-1011,
      doi:10.1093/mnras/stae1795.

## The grey opacity of a thick cell

Code: `grid.cc` and `rpkt.cc` (see [Code map](codemap.md)).

`input.txt` holds an optical depth limit and a number of grey timesteps (see [Plan a simulation](planning.md)). In
these first timesteps, a cell with a grey optical depth to the outer boundary at or above the limit is thick. An
r-packet in a thick cell gets a grey opacity, and a k-packet there emits an r-packet from a Planck function.

- Options: [RPKT_GREY_TYPE](options.md#rpkt_grey_type) selects the grey opacity.
- References: Tanaka, Kato, Gaigalas & Kawaguchi (2020), MNRAS, 496, 1369-1392, doi:10.1093/mnras/staa1576 and Just,
  Kullmann, Goriely, Bauswein, Janka & Collins (2022), MNRAS, 510, 2820-2840, doi:10.1093/mnras/stab3327 give the
  fits of two of the three types.

## The polarisation and the virtual packets

Code: `vpkt.cc` (see [Code map](codemap.md)).

With [POL_ON](options.md#pol_on), the code follows the Stokes parameters of the packets. With
[VPKT_ON](options.md#vpkt_on), it also sends {term}`virtual packets <virtual packet>` towards the observer directions of `vpkt.txt` at an
emission or a scattering event of a real packet. An event in a cell with the grey opacity
sends none, and `vpkt.txt` gives the range of the arrival times and the wavelength ranges of the events that send
them. `vpkt.txt` can also give a limit of the grey optical depth. A cell above that limit sends no virtual packet,
also after the timesteps of the grey opacity.

Each virtual packet adds its rest-frame energy to the spectrum of its direction with the weight `pn * exp(-tau)`.
`tau` is the optical depth from the event to the observer. `pn` is the probability for each unit of solid angle that
the event sends the packet in the observer direction:

- An emission is isotropic in the comoving frame, so `pn` is `1 / (4 pi)`.
- An electron scattering event has the dipole phase function, so `pn` is
  `3 / (16 pi) * (1 + mu^2 + (mu^2 - 1) * q)`. `mu` is the cosine of the scattering angle, and `q` is the Stokes
  parameter Q/I of the packet in the scattering plane.

The code divides `pn` by the square of the Doppler factor `nu_cmf / nu_rf` of the observer direction. That factor
transforms the solid angle from the comoving frame to the rest frame. The sum of the virtual packets gives the
spectrum of each direction with less noise than the real packets give.

- Options: [POL_ON](options.md#pol_on), [VPKT_ON](options.md#vpkt_on), and
  [VPKT_WRITE_CONTRIBS](options.md#vpkt_write_contribs).
- Reference: Bulla, Sim & Kromer (2015), MNRAS, 450, 967-981, doi:10.1093/mnras/stv657.

## The binned opacities and the thermalisation probability

Code: `rpkt.cc` and `vpkt.cc` (see [Code map](codemap.md)).

The binned opacities replace the line-by-line opacities with one opacity for each wavelength bin. Each line adds a
weight to the sum of its bin, and an option selects the weight from the Sobolev optical depth of the line. A
thermalisation probability for each bound-bound absorption can replace the macro atom. With that probability, the
packet is emitted again at once as an r-packet with a new frequency from the thermal emission of the cell. It does
not become a k-packet. With the remaining probability, the packet scatters and keeps its comoving frequency.

- Options:
    - [RPKT_USE_BINNED_OPACITIES](options.md#rpkt_use_binned_opacities) and
      [VPKT_USE_BINNED_OPACITIES](options.md#vpkt_use_binned_opacities) select the binned opacities for the real and
      for the virtual packets.
    - [BINNED_OPACITY_LINE_WEIGHT](options.md#binned_opacity_line_weight) and
      [LINE_BINNED_CAPPED_TAU_MAX](options.md#line_binned_capped_tau_max) select the weight of a line.
    - [RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY](options.md#rpkt_boundbound_thermalisation_probability) sets the
      thermalisation probability.
- References:
    - the {term}`expansion opacity`: Eastman & Pinto (1993), ApJ, 412, 731-751, doi:10.1086/172957;
    - the line-binned opacity: Fontes, Fryer, Hungerford, Wollaeger & Korobkin (2020), MNRAS, 493, 4143-4171,
      doi:10.1093/mnras/staa485;
    - these modes in ARTIS: Shingles et al. (in prep).
