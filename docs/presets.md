# Option presets

`artisoptions.h` holds the compile-time options of a run (see [Compile-time options](options.md)). It is a symlink to
one of the `artisoptions_*.h` {term}`presets <preset>`, and each preset is a complete set of values for one type of simulation. This
page describes the presets of the four types below. Select one with a symlink:

```sh
ln -s artisoptions_nltenebular.h artisoptions.h
```

The build needs this symlink, and `artisoptions.h` is gitignored. Do not copy the preset. The symlink keeps the options
correct when the preset changes, e.g. after a change of branch.

Start from the preset that is nearest to your simulation, and then change the options that your model needs.

## Which preset

| Preset | Transient | Phase | Ionisation and excitation | Non-thermal electrons | Particle thermalisation |
| --- | --- | --- | --- | --- | --- |
| `classic` | Type Ia supernova | photospheric | LTE excitation, tabulated photoionisation | off | instant |
| `kilonova_lte` | kilonova | all | LTE | off | time dependent |
| `nltenebular` | supernova | nebular | NLTE, multibin radiation field | Spencer-Fano | instant |
| `nltephotospheric` | supernova | photospheric | NLTE, multibin radiation field | Spencer-Fano | instant |

Most presets hold a value of `NUM_PACKETS` for a production run with approximately 1000 MPI {term}`ranks <rank>`. The
`nltephotospheric` preset has a much smaller value. Change the value for your rank count, because each rank keeps its
share of the packets in memory.

## The options that differ

This table holds each option whose value differs between the presets, and `scripts/make_options_docs.py` writes it.
[Compile-time options](options.md) describes each option.

<!-- begin the table that scripts/make_options_docs.py writes -->

| Option | classic | kilonova_lte | nltenebular | nltephotospheric |
| --- | --- | --- | --- | --- |
| [NUM_PACKETS](options.md#num_packets) | `5e8` | `2e8` | `1e9` | `9.6e6` |
| [INITIAL_PACKETS_ON](options.md#initial_packets_on) | `true` | `true` | `false` | `true` |
| [CUBOID_NCOORDGRID_X](options.md#cuboid_ncoordgrid_x) | `100` | `50` | `50` | `50` |
| [CUBOID_NCOORDGRID_Y](options.md#cuboid_ncoordgrid_y) | `100` | `50` | `50` | `50` |
| [CUBOID_NCOORDGRID_Z](options.md#cuboid_ncoordgrid_z) | `100` | `50` | `50` | `50` |
| [USE_CALCULATED_MEANATOMICWEIGHT](options.md#use_calculated_meanatomicweight) | `false` | `true` | `false` | `false` |
| [SINGLE_LEVEL_TOP_ION](options.md#single_level_top_ion) | `true` | `false` | `false` | `false` |
| [PHIXS_CLASSIC_NO_INTERPOLATION](options.md#phixs_classic_no_interpolation) | `true` | `false` | `false` | `false` |
| [RATECOEFF_TABLESIZE](options.md#ratecoeff_tablesize) | `100` | `200` | `100` | `100` |
| [MINTEMP](options.md#mintemp) | `3500.` | `500.` | `1000.` | `3500.` |
| [MAXTEMP](options.md#maxtemp) | `140000.` | `150000.` | `30000.` | `140000.` |
| [MINPOP](options.md#minpop) | `1e-30` | `1e-40` | `1e-40` | `1e-40` |
| [ION_NLEVELS_EXCITED_NLTE](options.md#ion_nlevels_excited_nlte) | `{ return 0; }` | `{ return 0; }` | see the option | see the option |
| [NLEVELS_ADD_MISSING_COLLISIONAL_TRANSITIONS](options.md#nlevels_add_missing_collisional_transitions) | `{ return 0; }` | see the option | see the option | see the option |
| [LTEPOP_EXCITATION_USE_TJ](options.md#ltepop_excitation_use_tj) | `true` | `false` | `false` | `false` |
| [FORCE_SAHA_ION_BALANCE](options.md#force_saha_ion_balance) | `{ return false; }` | see the option | `{ return false; }` | `{ return false; }` |
| [NLTE_TE_NNE_USE_ANDERSON_ACCEL](options.md#nlte_te_nne_use_anderson_accel) | `false` | `true` | `false` | `false` |
| [NLTE_TE_NNE_RELTOL](options.md#nlte_te_nne_reltol) | `0.04` | `0.01` | `0.04` | `0.04` |
| [NLTE_USE_GTH_SOLVER](options.md#nlte_use_gth_solver) | `false` | `true` | `true` | `false` |
| [STRICT_POPULATION_CHECKING](options.md#strict_population_checking) | `false` | `true` | `true` | `true` |
| [STRICT_POPULATION_CHECKING_INVERSION_FACTOR_SOLVER_FAIL](options.md#strict_population_checking_inversion_factor_solver_fail) | `1000.` | `50.` | `50.` | `50.` |
| [STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING](options.md#strict_population_checking_inversion_factor_printout_warning) | `10.` | `2.` | `2.` | `2.` |
| [NLTE_LIMIT_ION_STAGES_MAX_LEVELPOP_OVER_ELEMENTPOP_REMOVE_ION](options.md#nlte_limit_ion_stages_max_levelpop_over_elementpop_remove_ion) | `no value` | `1e-9` | `1e-9` | `1e-9` |
| [COLLISIONAL_HEATING_FROM_LEVELPOPS](options.md#collisional_heating_from_levelpops) | `false` | `true` | `true` | `true` |
| [BFCOOLING_USELEVELPOPNOTIONPOP](options.md#bfcooling_uselevelpopnotionpop) | `false` | `false` | `true` | `true` |
| [USE_ION_BFHEATING_ESTIMATORS](options.md#use_ion_bfheating_estimators) | `true` | `false` | `false` | `false` |
| [USE_LUT_PHOTOION](options.md#use_lut_photoion) | `true` | `false` | `false` | `false` |
| [MULTIBIN_RADFIELD_FIRST_TIMESTEP](options.md#multibin_radfield_first_timestep) | `no value` | `no value` | `12` | `12` |
| [DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP](options.md#detailed_bf_estimators_first_timestep) | `no value` | `no value` | `13` | `13` |
| [DETAILED_BF_ESTIMATORS_LEVEL_SELECTED](options.md#detailed_bf_estimators_level_selected) | `{ return false; }` | `{ return false; }` | `{ return true; }` | see the option |
| [NT_SCHEME](options.md#nt_scheme) | `NonThermalScheme::NT_OFF` | `NonThermalScheme::NT_OFF` | `NonThermalScheme::NT_SPENCERFANO` | `NonThermalScheme::NT_SPENCERFANO` |
| [PARTICLE_THERMALISATION_SCHEME](options.md#particle_thermalisation_scheme) | see the option | see the option | see the option | see the option |
| [NU_MIN_R](options.md#nu_min_r) | `1e14` | `1e13` | `1e13` | `1e13` |
| [NU_MAX_R](options.md#nu_max_r) | `5e15` | `5e16` | `5e15` | `5e15` |
| [DIPOLE](options.md#dipole) | `true` | `false` | `false` | `false` |
| [POL_ON](options.md#pol_on) | `true` | `false` | `false` | `false` |
| [USE_RELATIVISTIC_DOPPLER_SHIFT](options.md#use_relativistic_doppler_shift) | `false` | `true` | `false` | `false` |
| [RPKT_GREY_TYPE](options.md#rpkt_grey_type) | `RpktGreyType::FEGROUP_APPROX` | `RpktGreyType::TANAKA2020_ELECTRONFRAC` | `RpktGreyType::FEGROUP_APPROX` | `RpktGreyType::FEGROUP_APPROX` |

<!-- end the table that scripts/make_options_docs.py writes -->

## classic

The historic mode of ARTIS for Type Ia supernovae, as
[Kromer & Sim (2009), MNRAS, 398, 1809-1826, doi:10.1111/j.1365-2966.2009.15256.x](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1809K/abstract)
describes it. The preset uses the tabulated photoionisation coefficients (`USE_LUT_PHOTOION`), {term}`LTE` excitation with the
radiation temperature (`LTEPOP_EXCITATION_USE_TJ`), one level for the top ion, and no non-thermal solver. Polarisation
is on, and the packets get a dipole scattering distribution. With `GRID_TYPE_OVERRIDE` set to `GridType::CARTESIAN3D`, the 3D grid
of a 1D or a 2D model has 100 cells per coordinate.

## kilonova_lte

Kilonovae in LTE, as
[Shingles et al. (2023), ApJL, 954, L41, doi:10.3847/2041-8213/acf29a](https://ui.adsabs.harvard.edu/abs/2023ApJ...954L..41S/abstract)
describes the decay and the {term}`thermalisation`. The preset has the time-dependent thermalisation of the decay particles,
a {term}`grey opacity` from a fit to the electron fraction, the relativistic Doppler shift, and a low `MINTEMP` for the
late ejecta. The frequency range of the {term}`r-packets <r-packet>` reaches further into the ultraviolet than the supernova presets.

## nltenebular

The nebular phase of a supernova, as
[Shingles et al. (2020), MNRAS, 492, 2029-2043, doi:10.1093/mnras/stz3412](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.2029S/abstract)
describes the method. The preset has these methods:

- the {term}`NLTE` populations, with 197 excited levels of Fe II and 80 of each other ion;
- the Spencer-Fano solver;
- the multibin radiation field model from {term}`timestep` 12;
- the detailed bound-free {term}`estimators <estimator>` from timestep 13.

`MAXTEMP` is low, and the cells get no initial packets.

[Atomic and model data](data.md) names the dataset for a production run with this preset.

## nltephotospheric

The photospheric phase of a supernova with the NLTE populations, the Spencer-Fano solver, and the multibin radiation
field model. The preset keeps 100 excited levels for each ion of an element below calcium, and 200 for calcium and the heavier elements.
