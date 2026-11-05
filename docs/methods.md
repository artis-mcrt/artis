# Methods and references

This page gives each method of ARTIS: a short description, the options that select it, and its reference. Use it when
you describe the method of your results. [How a simulation works](simulation.md) gives the sequence of the steps, and
[Compile-time options](options.md) describes each option in full.

The log file of rank 0 holds the git version and the options of a run (see [Running sn3d](running.md)).
`artistools describe <model folder>` gives the main options of a run together with their effect.

## The energy packets

The radiation field is a set of indivisible energy packets. A packet changes its type when its energy changes its
form, and the code follows it until it escapes or the simulation ends. [How a simulation works](simulation.md) gives
the types of the packets.

- Options: [NUM_PACKETS](options.md#num_packets), [UNIFORM_PELLET_ENERGIES](options.md#uniform_pellet_energies),
  and [INITIAL_PACKETS_ON](options.md#initial_packets_on).
- References: Lucy (2002), A&A, 384, 725-735, doi:10.1051/0004-6361:20011756 describes the packets. An early version
  of the code is described in Sim (2007), MNRAS, 375, 154-162, doi:10.1111/j.1365-2966.2006.11271.x and Kromer & Sim
  (2009), MNRAS, 398, 1809-1826, doi:10.1111/j.1365-2966.2009.15256.x.

## The opacity of the r-packets

An r-packet moves until it reaches a line, a continuum absorption, or an electron scattering. The default mode gives
each line its own Sobolev optical depth. The expansion of the ejecta shifts the comoving frequency of the packet
through the lines one after the other. The continuum opacity is the sum of the bound-free and the free-free opacities.

- Options:
    - [USE_RELATIVISTIC_DOPPLER_SHIFT](options.md#use_relativistic_doppler_shift) selects the full Doppler factor.
    - [PHIXS_CLASSIC_NO_INTERPOLATION](options.md#phixs_classic_no_interpolation) selects the table lookup of the
      photoionisation cross section.
    - [DIPOLE](options.md#dipole) selects the direction of an electron scattering.
    - [NU_MIN_R](options.md#nu_min_r) and [NU_MAX_R](options.md#nu_max_r) set the frequency limits of the spectra.
- The alternatives are the [binned opacities](#the-binned-opacities-and-the-thermalisation-probability) and the
  [grey opacity of a thick cell](#the-grey-opacity-of-a-thick-cell).

## The macro atom and the k-packets

After an absorption, a macro atom makes internal radiative and collisional transitions between the levels and the
ions. It stops when it emits an r-packet or becomes a k-packet. A k-packet holds thermal energy. The code adds the
cooling rate of each process of the cell, and it samples one of them. That process converts the k-packet into
radiation or into a macro atom. The cooling processes are the free-free emission, the free-bound emission, the
collisional excitation, and the collisional ionisation.

- Options: [RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY](options.md#rpkt_boundbound_thermalisation_probability)
  replaces the macro atom when it has a value (see
  [the binned opacities](#the-binned-opacities-and-the-thermalisation-probability)).
- References: Lucy (2002), A&A, 384, 725-735, doi:10.1051/0004-6361:20011756 and Lucy (2003), A&A, 403, 261-275,
  doi:10.1051/0004-6361:20030357.

## The gas state

The grid update solves the state of each model grid cell that contains matter: the electron temperature T_e, the
electron density, and the ion and level populations. Kromer & Sim (2009), MNRAS, 398, 1809-1826,
doi:10.1111/j.1365-2966.2009.15256.x describes the classic mode. Shingles et al. (2020), MNRAS, 492, 2029-2043,
doi:10.1093/mnras/stz3412 describes the mode with the non-LTE (NLTE) populations. [Option presets](presets.md) says
which preset selects which mode.

### The populations without NLTE levels

The first `num_lte_timesteps` timesteps of `input.txt` solve the ionisation balance in local thermodynamic equilibrium
(LTE). After these timesteps, an element without NLTE levels gets a photoionisation equilibrium for its ions and
Boltzmann populations for its levels.

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

The solver puts the radiative rates, the collisional rates, and the non-thermal ionisation rates of all the levels of
all the ions of an element into one rate matrix. The solution of the statistical equilibrium gives the level
populations. Each ion has a ground level, a number of excited levels in full NLTE, and one superlevel for the levels
above them. The cell repeats the solutions of the populations, of T_e, and of the non-thermal electrons until T_e and
the electron density change by less than a tolerance.

If the solution of an element fails, the solver removes one ion and solves again, or the element falls back to
LTE. The log of the rank has a warning for each fallback to LTE.

- Options:
    - [ION_NLEVELS_EXCITED_NLTE](options.md#ion_nlevels_excited_nlte) sets the number of the excited levels in full
      NLTE.
    - [NLEVELS_ADD_MISSING_COLLISIONAL_TRANSITIONS](options.md#nlevels_add_missing_collisional_transitions) adds
      collisional transitions where the data has none.
    - [NLTE_USE_GTH_SOLVER](options.md#nlte_use_gth_solver) selects the Grassmann-Taksar-Heyman state elimination
      in place of the LU decomposition.
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
    - the time term: Pognan, Jerkstrand & Grumer (2022), MNRAS, 510, 3806-3837, eqs. 8 and 17,
      doi:10.1093/mnras/stab3674;
    - the charge transfer reactions: Kingdon & Ferland (1996), ApJS, 106, 205-211, doi:10.1086/192335 and Butler &
      Dalgarno (1980), ApJ, 241, 838-843, doi:10.1086/158395.

### The thermal balance

The code solves heating = cooling for the electron temperature of each cell. The heating rates are the bound-free
heating, the free-free heating, the collisional de-excitation, and the non-thermal heating. The code keeps T_e
between [MINTEMP](options.md#mintemp) and [MAXTEMP](options.md#maxtemp). A cell at one of these limits has no thermal
balance, and the log of its rank has a warning (see [Troubleshooting](troubleshooting.md#the-log)).

- Options: [COLLISIONAL_HEATING_FROM_LEVELPOPS](options.md#collisional_heating_from_levelpops),
  [BFCOOLING_USELEVELPOPNOTIONPOP](options.md#bfcooling_uselevelpopnotionpop),
  [USE_ION_BFHEATING_ESTIMATORS](options.md#use_ion_bfheating_estimators), and
  [WRITE_ION_HEATING_COOLING_RATES](options.md#write_ion_heating_cooling_rates).

## The radiation field model

The packets add to the estimators of the radiation field of each cell that they enter. The grid update fits a dilute
blackbody to the full spectrum of each cell, with a dilution factor W and a radiation temperature T_R. With the
multibin model, it also fits a dilute blackbody to each frequency bin. The rates of the photoionisation and of the
bound-free heating use these fits.

- Options:
    - [MULTIBIN_RADFIELD_FIRST_TIMESTEP](options.md#multibin_radfield_first_timestep) sets the first timestep that
      reads the fits of the bins.
    - [RADFIELDBINCOUNT](options.md#radfieldbincount), [RADFIELDBINS_NU_MIN](options.md#radfieldbins_nu_min),
      [RADFIELDBINS_NU_MAX](options.md#radfieldbins_nu_max), and
      [RADFIELDBINS_T_E_SUPERBIN_NU_MAX](options.md#radfieldbins_t_e_superbin_nu_max) set the bins.
    - [DETAILED_LINE_ESTIMATORS_ON](options.md#detailed_line_estimators_on) stores the estimators of selected lines.
- Reference: Shingles et al. (2020), MNRAS, 492, 2029-2043, section 2.2, doi:10.1093/mnras/stz3412.

## The photoionisation rates and the bound-free estimators

The photoionisation rate coefficient of a level comes from one of three sources:

- a table for a Planck radiation field, which the code multiplies by the dilution factor W and by the ratio of a Monte
  Carlo estimator to the blackbody rate ([USE_LUT_PHOTOION](options.md#use_lut_photoion));
- an integral of the cross section over the radiation field model of the cell;
- a detailed bound-free estimator, which the packets build for each selected continuum
  ([DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP](options.md#detailed_bf_estimators_first_timestep) and
  [DETAILED_BF_ESTIMATORS_LEVEL_SELECTED](options.md#detailed_bf_estimators_level_selected)).

If the model folder holds `recombrates.txt`, the code calibrates the total recombination rate of each ion to that file
at the temperature [RECOMBCALIBRATION_T_ELEC](options.md#recombcalibration_t_elec).

- Reference: Shingles et al. (2020), MNRAS, 492, 2029-2043, section 2.2, doi:10.1093/mnras/stz3412.

## The non-thermal electrons

The fast leptons of the decays lose their energy to the heating, the ionisation, and the excitation of the gas. The
Spencer-Fano solver finds the energy spectrum of these electrons on a uniform energy grid. The spectrum gives the
three fractions of the deposited energy, and the non-thermal ionisation and excitation rates for the NLTE solver.
The source term of the Auger electrons is an addition to the equation. The work function approximation is the
alternative. It gives ionisation rates and no excitation rates.

- Options:
    - [NT_SCHEME](options.md#nt_scheme) selects the Spencer-Fano solution, the work function approximation, or no
      non-thermal ionisation.
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

The decay network holds the alpha, the beta, and the fission decays. The code gets the abundance of each nuclide at
each time from the solutions of the Bateman equation. The `data` folder of the repository holds the decay data and
the gamma-ray lines of each nuclide (see [Atomic and model data](data.md)). The decay particles deposit their energy
at once, or the code transports them and they lose their energy over time.

- Options: [PARTICLE_THERMALISATION_SCHEME](options.md#particle_thermalisation_scheme) selects the instant
  deposition, a transport with the Monte Carlo method, or an analytic thermalisation efficiency.
- References:
    - the decays and the time-dependent thermalisation: Shingles et al. (2023), ApJL, 954, L41,
      doi:10.3847/2041-8213/acf29a;
    - the analytic efficiencies: Barnes, Kasen, Wu & Martínez-Pinedo (2016), ApJ, 829, 110,
      doi:10.3847/0004-637X/829/2/110 and Wollaeger, Korobkin, Fontes, Rosswog, Even & Fryer (2018), MNRAS, 478,
      3298-3334, doi:10.1093/mnras/sty1018.

## The gamma rays

A gamma packet gives its energy to the gas by Compton scattering, by photoelectric absorption, or by pair production.
Every preset transports the gamma packets with frequency-dependent opacities. The alternatives are a grey opacity
and an analytic thermalisation efficiency, which deposits a fraction of the gamma energy at once.

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

`input.txt` holds an optical depth limit and a number of grey timesteps (see [Plan a simulation](planning.md)). In
these first timesteps, a cell with a grey optical depth to the outer boundary at or above the limit is thick. An
r-packet in a thick cell gets a grey opacity, and a k-packet there emits an r-packet from a Planck function.

- Options: [RPKT_GREY_TYPE](options.md#rpkt_grey_type) selects the grey opacity.
- References: Tanaka, Kato, Gaigalas & Kawaguchi (2020), MNRAS, 496, 1369-1392, doi:10.1093/mnras/staa1576 and Just,
  Kullmann, Goriely, Bauswein, Janka & Collins (2022), MNRAS, 510, 2820-2840, doi:10.1093/mnras/stab3327 give the
  fits of two of the three types.

## The polarisation and the virtual packets

The code follows the Stokes parameters of the packets. At each emission and scattering event of a real packet, it
also sends virtual packets towards the observer directions of `vpkt.txt`. The energy of each virtual packet gets the
weight of its optical depth to the observer. The sum gives the spectrum of each direction with less noise than the
real packets give.

- Options: [POL_ON](options.md#pol_on), [VPKT_ON](options.md#vpkt_on), and
  [VPKT_WRITE_CONTRIBS](options.md#vpkt_write_contribs).
- Reference: Bulla, Sim & Kromer (2015), MNRAS, 450, 967-981, doi:10.1093/mnras/stv657.

## The binned opacities and the thermalisation probability

The binned opacities replace the line-by-line opacities with one opacity for each wavelength bin. Each line adds a
weight to the sum of its bin, and an option selects the weight from the Sobolev optical depth of the line. A
thermalisation probability for each bound-bound absorption can replace the macro atom. The absorbed energy then
becomes a k-packet with that probability, and the packet scatters with the remaining probability.

- Options:
    - [RPKT_USE_BINNED_OPACITIES](options.md#rpkt_use_binned_opacities) and
      [VPKT_USE_BINNED_OPACITIES](options.md#vpkt_use_binned_opacities) select the binned opacities for the real and
      for the virtual packets.
    - [BINNED_OPACITY_LINE_WEIGHT](options.md#binned_opacity_line_weight) and
      [LINE_BINNED_CAPPED_TAU_MAX](options.md#line_binned_capped_tau_max) select the weight of a line.
    - [RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY](options.md#rpkt_boundbound_thermalisation_probability) sets the
      thermalisation probability.
- References:
    - the expansion opacity: Eastman & Pinto (1993), ApJ, 412, 731-751, doi:10.1086/172957;
    - the line-binned opacity: Fontes, Fryer, Hungerford, Wollaeger & Korobkin (2020), MNRAS, 493, 4143-4171,
      doi:10.1093/mnras/staa485;
    - these modes in ARTIS: Shingles et al. (in prep).
