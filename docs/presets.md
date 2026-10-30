# Option presets

`artisoptions.h` holds the compile-time options of a run (see [Compile-time options](options.md)). The repository has
five presets of these options. Each preset is a complete set of values for one type of simulation. Select one with a
symlink:

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
| `kilonova_nlte` | kilonova | all | NLTE for Sr, Y, and Zr | Spencer-Fano | time dependent |
| `nltenebular` | supernova | nebular | NLTE, multibin radiation field | Spencer-Fano | instant |
| `nltephotospheric` | supernova | photospheric | NLTE, multibin radiation field | Spencer-Fano | instant |

Each preset holds a value of `NUM_PACKETS` for a production run with approximately 1000 MPI ranks. The
`nltephotospheric` preset has a much smaller value. Change the value for your rank count, because each rank keeps its
share of the packets in memory.

## classic

The historic mode of ARTIS for Type Ia supernovae, as
[Kromer & Sim (2009), MNRAS, 398, 1809-1826, doi:10.1111/j.1365-2966.2009.15256.x](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1809K/abstract)
describes it. The preset uses the tabulated photoionisation coefficients (`USE_LUT_PHOTOION`), LTE excitation with the
radiation temperature (`LTEPOP_EXCITATION_USE_TJ`), one level for the top ion, and no non-thermal solver. Polarisation
is on, and the packets get a dipole scattering distribution. The 3D grid of a 1D or a 2D model has 100 cells per
coordinate.

The tests `classicmode_1d_3dgrid` and `classicmode_3d` use this preset.

## kilonova_lte

Kilonovae in LTE, as
[Shingles et al. (2023), ApJL, 954, L41, doi:10.3847/2041-8213/acf29a](https://ui.adsabs.harvard.edu/abs/2023ApJ...954L..41S/abstract)
describes the decay and the thermalisation. The preset has the time-dependent thermalisation of the decay particles,
a grey opacity from a fit to the electron fraction, the relativistic Doppler shift, and a low `MINTEMP` for the
late ejecta. The frequency range of the r-packets reaches further into the ultraviolet than the supernova presets.

The tests `kilonova_1d`, `kilonova_2d`, `kilonova_2d_expansionopac`, `kilonova_2d_barnesthermalisation`, and
`kilonova_2d_xcomgammaphotoion` use this preset. `kilonova_1d` is the quickest test of the repository, so the
[Quickstart](quickstart.md) uses it.

## kilonova_nlte

The same kilonova setup as `kilonova_lte`, and additionally the NLTE populations of strontium, yttrium, and zirconium
(250 excited levels each) and the Spencer-Fano solver for the non-thermal electrons.

The test `kilonova_1d_timedepnlte` uses this preset.

## nltenebular

The nebular phase of a supernova, as
[Shingles et al. (2020), MNRAS, 492, 2029-2043, doi:10.1093/mnras/stz3412](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.2029S/abstract)
describes the method. The preset has the NLTE populations (197 excited levels of Fe II, and 80 of the other ions), the
Spencer-Fano solver, the multibin radiation field model from timestep 12, and the detailed bound-free estimators from
timestep 13. `MAXTEMP` is low, and the cells get no initial packets.

The tests `nebular_1d_3dgrid` and `nebular_1d_3dgrid_limitbfest` use this preset. `nebular_1d_3dgrid` covers the
physics that most changes touch, so it is the default test for a change of the code.

## nltephotospheric

The photospheric phase of a supernova with the NLTE populations, the Spencer-Fano solver, and the multibin radiation
field model. The preset keeps 100 excited levels for each ion below calcium and 200 for the heavier ions. `MAXTEMP` is
high for the hot early ejecta.

The test `nltephotospheric_dynamic_ion_range_1d_1dgrid` uses this preset.
