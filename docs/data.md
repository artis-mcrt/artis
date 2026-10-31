# Atomic and model data

A simulation reads three sets of data:

1. the physics data of the repository;
2. an atomic database;
3. a model of the ejecta.

[Input files](inputfiles.md) gives the format of each file. All of them go in the model folder.

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

### The datasets of the tests

The end-to-end tests download a reduced dataset from a release of this repository.
`tests/atomic_data_release.txt` names the release, and `tests/setupfuncs.sh` holds the download command.

| Archive | Elements | Tests |
| --- | --- | --- |
| `atomicdata_classic.tar.xz` | the classic dataset, with cross-sections in the v1 format | `classicmode_1d_3dgrid`, `classicmode_3d` |
| `atomicdata_feconi.tar.xz` | Fe, Co, Ni | `kilonova_1d`, `kilonova_2d*`, `nebular_1d_3dgrid*` |
| `atomicdata_hefeconi_fe_i_to_vii.tar.xz` | He, Fe, Co, Ni, with Fe I to Fe VII | `nltephotospheric_dynamic_ion_range_1d_1dgrid` |
| `atomicdata_sryzrlace.tar.zst` | Sr, Y, Zr, La, Ce | `kilonova_1d_timedepnlte` |

Each dataset holds few elements, so that a test runs in about ten minutes. The datasets of the tests are not fit for a
production run.

### A database for a production run

The maintainers will publish the full atomic databases. Until then, ask them for a database in an
[issue](https://github.com/artis-mcrt/artis/issues).

[artisatomic](https://github.com/artis-mcrt/artisatomic) converts atomic data of several published sources into the
ARTIS format. The `atomic_data_logs` folder of each archive holds the log of that conversion for each ion.

## The model of the ejecta

A run needs two model files, and some runs need more:

- `model.txt`: the grid, the density of each cell, and the nuclear composition;
- `abundances.txt`: the elemental mass fractions of each cell;
- `Ye.txt` and `gridcontributions.txt`: the electron fraction and the tracer particle contributions of a kilonova
  model;
- `recombrates.txt`: the tabulated recombination rates, e.g. for the nebular test.

The `tests/*_inputfiles` folders hold small models of each type: 1D spherical, 2D cylindrical, and 3D Cartesian. They
are the quickest examples of the format.

[artistools](https://github.com/artis-mcrt/artistools) makes an ARTIS model from the output of a hydrodynamic
simulation or of a nuclear network. `artistools inputmodel` lists the converters, e.g. from CMFGEN, from the tracer
particles of an end-to-end kilonova simulation, and from a single trajectory. `artistools inputmodel describe` gives
the mass, the velocity structure, and the abundances of a model.
