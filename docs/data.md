# Atomic and model data

A simulation reads three sets of data:

1. the physics data of the repository;
2. an atomic database;
3. a model of the ejecta.

[Input files](inputfiles.md) gives the format of each file. The atomic database and the model go in the {term}`model folder`.

## The physics data of the repository

`data/` holds the data that comes with the code:

- the decay data of the alpha, beta-minus, and fission decays;
- the gamma-ray spectrum of each nuclide;
- the collisional ionisation and Auger tables;
- the electron binding energies;
- the rates of the charge transfer reactions.

Keep the source code in a subfolder `artis` of the model folder. The code then finds these files in `artis/data`, and
you copy nothing.

## The atomic database

A run reads the levels, the lines, and the photoionisation cross-sections of each ion from the model folder:

- `adata.txt`: the energy levels of each ion;
- `transitiondata.txt`: the lines of each ion;
- `phixsdata_v2.txt` or `phixsdata.txt`: the photoionisation cross-sections;
- `compositiondata.txt`: the filter on the elements, the ion stages, and the level counts.

`compositiondata.txt` selects what the code reads from a database, so one database serves models of different
composition. A larger level count gives a slower run and a larger memory use.

The end-to-end tests use reduced datasets with few elements (see [Testing](testing.md)). They are not fit for a
production run.

### A database for a production run

For the nebular phase of a Type Ia supernova, use the dataset of
[Shingles et al. (2022), MNRAS, 512, 6150-6163, doi:10.1093/mnras/stac902](https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.6150S/abstract).
Download it into the model folder and unpack it there:

```sh
curl -fLO https://github.com/artis-mcrt/artis/releases/download/v2026.5.15/atomicdata_shinglesetal2022.tar.zst
tar --zstd -xf atomicdata_shinglesetal2022.tar.zst
```

GNU tar calls the `zstd` command for a `.zst` archive. The package of the library, e.g. `libzstd-dev`, supplies no such
command, so install the package `zstd` as well.

The dataset holds 41 ions of the 11 elements of a Type Ia model:

- C I to C IV, O I to O IV, Si I to Si IV, S I to S IV, Ar I to Ar IV, and Ca I to Ca IV;
- Ne I to Ne III, and Mg I to Mg II;
- Fe I to Fe V, Co II to Co IV, and Ni II to Ni V.

Use it with the `nltenebular` {term}`preset` of [Option presets](presets.md). The dataset also holds `recombrates.txt`, so the
code calibrates the recombination rate of each ion at `RECOMBCALIBRATION_T_ELEC`.

Column 5 of its `compositiondata.txt` is -1, so a run reads every level of each ion. Put a level count there for a
quicker run. A build without libzstd reads no compressed input file, so decompress each `.zst` file of the archive.

The maintainers will publish the databases of the other model types. Until then, make a database with
[artisatomic](https://github.com/artis-mcrt/artisatomic), which converts atomic data of several published sources into
the ARTIS format. An archive that artisatomic made holds an `atomic_data_logs` folder with the log of that conversion for each ion.

## The model of the ejecta

A run needs two model files, and some runs need more:

- `model.txt`: the grid, the density of each cell, and the nuclear composition;
- `abundances.txt`: the elemental mass fractions of each cell;
- `Ye.txt`: the electron fraction of a kilonova model;
- `recombrates.txt`: the tabulated recombination rates, e.g. for the nebular test.

artistools also writes `gridcontributions.txt` for a kilonova model, but `sn3d` does not read that file.

The `tests/*_inputfiles` folders hold small models of each type: 1D spherical, 2D cylindrical, and 3D Cartesian. They
are the quickest examples of the format.

[artistools](https://github.com/artis-mcrt/artistools) makes an ARTIS model from the output of a hydrodynamic
simulation or of a nuclear network. `artistools inputmodel` lists the converters, e.g. from CMFGEN, from the tracer
particles of an end-to-end kilonova simulation, and from a single trajectory. `artistools inputmodel describe` gives
the mass, the velocity structure, and the abundances of a model.
