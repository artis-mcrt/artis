# Input files

Put these files in the {term}`model folder`. Keep the ARTIS source folder, or a symlink to it, in the model folder with the
name `artis`. The code looks for an input file in the model folder, then in `data/`, and then in `artis/data/`. It
thus finds the physics data that comes with the source, and you copy nothing (see
[Atomic and model data](data.md)).

A build with libzstd also reads a text input file in zstd compressed form, e.g. `model.txt.zst`,
`transitiondata.txt.zst`, or `packets/packets00_0000.out.zst` for `exspec`. The plain file has priority when both
exist.

| File | Necessary | Content |
| --- | --- | --- |
| [input.txt](#inputtxt) | yes, or `input-newrun.txt` | the run-time parameters |
| [model.txt](#modeltxt) | yes | the grid, the densities, and the nuclide mass fractions |
| [abundances.txt](#abundancestxt) | yes | the element mass fractions of each cell |
| [compositiondata.txt](#compositiondatatxt) | yes | the elements, the ion stages, and the level counts of the run |
| [adata.txt](#adatatxt) | yes | the energy levels of each ion |
| [transitiondata.txt](#transitiondatatxt) | yes | the lines of each ion |
| [phixsdata_v2.txt or phixsdata.txt](#phixsdata_v2txt-and-phixsdatatxt) | one of them | the photoionisation cross-sections |
| [vpkt.txt](#vpkttxt) | with `VPKT_ON` | the observer directions of the virtual packets |
| [Ye.txt](#yetxt) | with `RPKT_GREY_TYPE = TANAKA2020_ELECTRONFRAC`, if `model.txt` has no `Ye` column | the electron fraction of each cell |
| [recombrates.txt](#recombratestxt) | no | the recombination rate coefficients for the calibration |
| [autoion.txt](#autoiontxt) | no | the autoionisation rates |
| [gamma_*.txt](#the-gamma-ray-line-files) | for seven nuclides | the gamma-ray lines of each nuclide, which come with the source in `artis/data` |

In the layouts of this page, `<name>` is a value that you supply. Where this page says that a file can have comment
lines, `sn3d` ignores each blank line and each line whose first character that is not a space is `#`.

## input.txt

`input.txt` holds the run-time parameters. The format is positional: each line has a fixed meaning, and the unused
lines must stay. The file can have comment lines. `sn3d` ignores the text after the values of a line, so each line
can end with a comment.

The file has 24 lines that are not comment lines. The line numbers of this table count only these lines, and they
start at 0:

| Line | Values | Meaning |
| --- | --- | --- |
| 0 | `pre_zseed` | The seed of the random number generator. A positive value is the seed. For zero or a negative value, `sn3d` makes a random seed. Reproducible results need a fixed seed, and a build with `REPRODUCIBLE=ON` stops without a positive value. |
| 1 | `ntimesteps` | The number of timesteps of the model. |
| 2 | `timestep_start timestep_finish` | The first timestep of this job, and the stop timestep. The job does not run the stop timestep. |
| 3 | `tmin_days tmax_days` | The start time and the end time of the model in days. |
| 16 | `simulation_continued_from_saved` | The value must be 0 or 1. 0 starts a new simulation. 1 continues from the restart files of a previous job. With `timestep_start` 0, the run is always a new simulation. |
| 18 | `num_lte_timesteps` | The number of timesteps at the start of the model that use LTE populations. The value must be 1 or more. |
| 19 | `optical_depth_is_thick num_grey_timesteps` | The condition for the grey opacity treatment of a cell (see below). |
| 21 | `nprocs_exspec` | The number of packet files that `exspec` reads, one for each `sn3d` rank. `sn3d` sets this line itself. |

`sn3d` does not use the other 16 lines, but they must stay. Start from `scripts/input-template.txt`, or see
[tests/classicmode_3d_inputfiles/input-newrun.txt](https://github.com/artis-mcrt/artis/blob/main/tests/classicmode_3d_inputfiles/input-newrun.txt)
for a complete example with a comment on every line.

The {term}`timesteps <timestep>`:

- `timestep_start` must be 0 or more, and it must be less than `ntimesteps`. `timestep_finish` must not be less than
  `timestep_start`, and it must not be more than `ntimesteps`.
- Give the number of timesteps as the stop value for the whole model. Long simulations run in several jobs, and the
  wall-time mechanism (`sn3d -w`) advances the start timestep at each restart and keeps the stop value.
- A smaller stop value makes the simulation stop earlier, e.g. to examine the partial results. The simulation
  continues only after you edit `input.txt`. The last {term}`restart files <restart file>` of such a run hold the start of its last
  timestep. The continuation job therefore propagates that timestep again and writes its spectra and light curves
  again.
- `sn3d` writes the restart files at the start of each timestep after the first timestep of the job, so give a job at
  least two timesteps. A job of one timestep writes no restart files.
- `tmin_days` must be more than 0 and less than `tmax_days`. The snapshot time of `model.txt` must not be after
  `tmin_days`.

The {term}`LTE` timesteps and the grey cells:

- With a value for `MULTIBIN_RADFIELD_FIRST_TIMESTEP`, `num_lte_timesteps` must not be larger than that value.
- With the Spencer-Fano solver (`NT_SCHEME`), the first solution comes at the timestep `num_lte_timesteps` + 1.
- `sn3d` estimates the grey optical depth of each cell as its grey opacity, times its density, times the distance from
  the cell to the outer boundary of the model. A cell with a depth of
  `optical_depth_is_thick` or more gets the {term}`grey opacity` treatment. This applies to each timestep before the timestep
  `num_grey_timesteps`.

At the start of a new simulation, `sn3d` writes `input.txt` again with a comment on each line and with the value of
`nprocs_exspec`. It keeps the same content in `input-newrun.txt`. If `input.txt` is absent at the start of a run,
{term}`rank` 0 restores it from `input-newrun.txt`.

## model.txt

`model.txt` holds the grid, the density of each cell, and the mass fractions of the nuclides. `sn3d` finds the number
of dimensions from the first lines of the file.

A 1D model has spherical shells:

```text
<number of cells>
<t_model in days>
#inputcellid vel_r_max_kmps logrho X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48
<one row for each cell>
```

A 2D model has a cylindrical grid:

```text
<number of cells in r> <number of cells in z>
<t_model in days>
<vmax in cm/s>
#inputcellid pos_rcyl_mid pos_z_mid rho X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48
<one row for each cell>
```

A 3D model has a Cartesian grid of cubes:

```text
<number of cells>
<t_model in days>
<vmax in cm/s>
#inputcellid pos_x_min pos_y_min pos_z_min rho X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48
<one row for each cell>
```

The lines before the cells:

- Each number of cells must be 1 or more.
- The file can have comment lines before the number of cells and before `t_model`.
- `t_model` is the snapshot time of the model, which must be above 0. The densities and the positions of the file
  apply at this time.
- `vmax` is the velocity of the outer boundary of the grid. Put this line directly after the line of `t_model`. A 1D
  model has no such line, and its `vmax` is the outer velocity of the last cell. `vmax` must be above 0 and
  less than the speed of light.
- The corners of a 2D or a 3D grid move faster than `vmax`. The corner velocity is `vmax` times the square root of 2
  for a cylindrical grid, and `vmax` times the square root of 3 for a Cartesian grid. `sn3d` stops if that velocity
  is not less than the speed of light, unless the run has `FORCE_SPHERICAL_ESCAPE_SURFACE`.
- The number of cells of a 3D model must be a cube number, because the grid has the same number of cells on each axis.
- The header line names the columns. It is optional. The columns of the cell id, the 1D velocity, and the density
  must have the names and the column positions of the layouts above. The column `velocity_outer` is a second name for
  `vel_r_max_kmps`. The position columns are between the cell id and the density. `sn3d` accepts each name that
  starts with `pos_` for them, and it does not read the axis from the name.
- A file without a header line gets the columns of the layouts above. A row with two more values gets `X_Ni57` and
  `X_Co57` as the last two columns.

The columns of each cell:

| Column | Model | Meaning |
| --- | --- | --- |
| `#inputcellid` | all | The cell id. The ids start at 0 or at 1, and they increase by one for each cell. |
| `vel_r_max_kmps` | 1D | The velocity of the outer boundary of the shell [km/s]. Each value must be larger than the value of the cell before, and the first value must be above 0. |
| `logrho` | 1D | The base-10 logarithm of the density [g/cm^3] at `t_model`. A value of -90 or less gives an empty cell. |
| `pos_rcyl_mid`, `pos_z_mid` | 2D | The midpoint of the cell [cm] at `t_model`. |
| `pos_x_min`, `pos_y_min`, `pos_z_min` | 3D | The corner of the cell with the lowest coordinates [cm] at `t_model`. |
| `rho` | 2D, 3D | The density [g/cm^3] at `t_model`. A value of zero gives an empty cell, and a negative value stops the run. |
| `X_Fegroup` | all | The mass fraction of the iron-group elements. The value must not be above 1.001, which permits a rounding error. `sn3d` treats a negative value as it treats a negative `X_<nuclide>` value. |
| `X_<nuclide>` | all | The mass fraction of one nuclide at `t_model`, e.g. `X_Ni56`. The name has the element symbol and the mass number. The value must be in the range from 0 to 1. `sn3d` sets a negative value above -1e-6 to 0, and a value below that limit stops the run. |
| `Ye` or `cellYe` | all | The electron fraction. |
| `q` | all | The trapped radiation energy per mass [erg/g] at `t_model`. The value must not be negative. |

More rules for the cells:

- A cell has one line or two lines. If the first cell line has a value for each column, each cell has one line. If
  not, each cell has two lines: the first line holds the cell id, the position, and the density, and the second line
  holds the other columns.
- The cells of a 2D model go through the r index first and then through the z index. The z midpoints start at the
  negative boundary. `sn3d` calculates the midpoint of each cell from `vmax` and `t_model`, and it stops if a value
  of the file differs by more than 0.1 percent.
- The cells of a 3D model go through the x index first, then the y index, and then the z index. The three position
  values of a row can have the sequence x, y, z or the sequence z, y, x. `sn3d` compares the positions in both
  sequences with the values that it calculates from `vmax` and `t_model`. If each sequence has a difference of more
  than half of a cell width, the log has a warning.
- `sn3d` finds the first cell id from the first cell of `model.txt`. `abundances.txt` and `Ye.txt` must use the same
  ids. The log messages always give the 0-based cell index.
- The value of the column `q` must already include the adiabatic losses before `t_model`. `INITIAL_PACKETS_ON` uses
  the column for the initial temperature and for the initial packets.
- A column with a different name gets a warning in the log, and `sn3d` does not use it.

## abundances.txt

`abundances.txt` holds the mass fractions of the elements of each cell. These include the elements whose nuclide mass
fractions are not in `model.txt`. The file can have comment lines.

The file has one row for each cell of `model.txt`, also for an empty cell, in the same sequence:

```text
<cell id> <X of Z=1> <X of Z=2> <X of Z=3> ...
```

The values start at hydrogen and continue through the atomic numbers. A row can end after the last element that the
model needs. A value must not be negative. `sn3d` sets a negative value above -1e-6 to 0, and a value at or below
that limit stops the run.

For a 1D or a 2D model, `sn3d` divides the values of each cell by their sum, so the values can also be proportional to
the mass fractions, e.g. element densities. For a 3D model, `sn3d` uses the values without normalisation. A 3D cell
whose values do not sum to 1 within 2 percent gets a warning in the log.

The mass fraction of an element must not be less than the sum of the mass fractions of its nuclides in `model.txt`.
`sn3d` permits a difference of 0.01 for the rounding of the values.
A cell with a density, but with no mass in the elements of `compositiondata.txt` that have ions, becomes an empty
cell. An element with `nions` 0 does not keep the cell. The log gets a warning that names the elements of that cell.

## compositiondata.txt

`compositiondata.txt` selects the elements, the ion stages, and the levels that the code reads from the atomic data
files. The file can have comment lines, and the code ignores the text to the right of a `#` character.

```text
<number of elements>
0
0
<Z> <nions> <lowest ion stage> <highest ion stage> <nlevelsmax> <abundance> <mass>
```

The number of elements must be 1 or more. The second value and the third value are unused and must be 0. Then the
file has one row for each element:

| Column | Meaning |
| --- | --- |
| `Z` | The atomic number, which must be in the range from 1 to 118. |
| `nions` | The number of ion stages, which must not be negative. With the value 0, the element has no ions, and the code keeps only its nuclide abundances. |
| lowest ion stage | The lowest ion stage of the run. The ion stage 1 is the neutral atom, so the value must be 1 or more. |
| highest ion stage | The highest ion stage of the run. The value must not be above `Z` + 1, which is the bare nucleus. For an element with ions, `nions` must be equal to the highest ion stage - the lowest ion stage + 1. An element with `nions` 0 has no such rule. |
| `nlevelsmax` | The maximum number of levels that the code reads for each ion of the element. A negative value reads every level. The value 0 stops the run for an element with ions. With `SINGLE_LEVEL_TOP_ION`, the highest ion stage keeps one level. |
| abundance | Unused. The value must be 0 or more. |
| mass | The mean mass of the stable nuclides of the element in atomic mass units. The value must be above 0. |

## adata.txt

`adata.txt` holds the energy levels. It has one block for each ion. The blocks of the ions that the run uses must have
the sequence of `compositiondata.txt`, with the ion stages of an element in increasing order. `sn3d` goes past the
blocks of other ions.

The header of a block:

```text
<Z> <ion stage> <number of levels> <ionisation energy in eV>
```

The ionisation energy must not be negative.

Then one row for each level:

```text
<level index> <energy in eV> <statistical weight> <number of transitions> <level name>
```

- The energy of a level is relative to the ground level of the ion.
- The statistical weight of each level that the run keeps must be above 0.
- `sn3d` reads the number of transitions, and it does not use it. The level name is optional, because `sn3d` does
  not read the text after the four values.
- Give the levels of an ion in the sequence of their energy. A level with an energy below the energy of the level
  before gives a warning in the log, and the run continues.
- An ion that the run keeps must have one level or more.

The level indices start at 0 or at 1, and they must increase by one for each level of an ion. `sn3d` finds the
first index from the first level of `adata.txt`. `transitiondata.txt` and `phixsdata_v2.txt` must use the same
numbering. `phixsdata.txt` always numbers the levels
from 1. The log messages always give the 0-based level index.

## transitiondata.txt

`transitiondata.txt` holds the lines. It has one block for each ion, in the same sequence as `adata.txt`.

The header of a block:

```text
<Z> <ion stage> <number of transitions>
```

The number of transitions can be 0, and it must not be negative.

Then one row for each transition:

```text
<lower level index> <upper level index> <A value in 1/s> <collision strength> <forbidden flag>
```

- The lower level index must be less than the upper level index.
- A collision strength of zero or above is the effective collision strength of the transition. The value zero gives
  a collisional rate of zero. If the collision strength is not available, the value is -1.0 for a permitted transition
  and -2.0 for a forbidden transition.
- The forbidden flag is 1 for a forbidden transition and 0 for a permitted transition.
- An A value must not be negative.
- The rows of a block can have no comment lines between them.

The code also accepts an older format of four columns, without the collision strength and the forbidden flag:

```text
<transition index> <lower level index> <upper level index> <A value in 1/s>
```

The code finds the format of each block from the number of columns in its first row.

## phixsdata_v2.txt and phixsdata.txt

These files hold the photoionisation cross-sections. A run needs one of them, and it can read both. A level must
then have its table in only one of the two files. `phixsdata_v2.txt`
has the same frequency grid for every table, and `phixsdata.txt` has an energy table for each level.
`phixsdata.txt` can have comment lines. `phixsdata_v2.txt` can have them only before the header line of a table and
before a target line, and not before the grid values or between the cross-section values.

`phixsdata_v2.txt` starts with the grid:

```text
<NPHIXSPOINTS>
<NPHIXSNUINCREMENT>
```

Both values must be above 0. Point number i of a table, from 0, has the frequency nu_edge (1 + i
`NPHIXSNUINCREMENT`), where nu_edge is the threshold frequency of the level.

Then each table of `phixsdata_v2.txt` has this layout:

```text
<Z> <upper ion stage> <upper level index> <lower ion stage> <lower level index> <threshold energy in eV>
<number of target levels>
<target level index> <probability>
<NPHIXSPOINTS cross-sections in megabarns>
```

- The lower level is the level that absorbs the photon. The upper ion stage must be the lower ion stage + 1.
- A list of target levels must have one target level or more.
- An upper level index of 0 or more gives one target level. The table then has no row with the number of target
  levels and no rows of target levels.
- A negative upper level index, e.g. -1, gives a list of target levels. The probabilities of the target levels must
  sum to 1 within 0.01, and each target level must occur only once. Each probability must be above 0.
- Each target level must be a level that the run keeps for the upper ion. A target above the level limit of
  `compositiondata.txt` stops the run. An upper ion with one level from `SINGLE_LEVEL_TOP_ION` is the exception.
- `sn3d` reads the threshold energy and does not use it. The threshold comes from the level energies of `adata.txt`.

Each table of `phixsdata.txt` has this layout:

```text
<Z> <upper ion stage> <upper level index> <lower ion stage> <lower level index> <number of points>
<energy above the threshold> <cross-section in megabarns>
```

- The second line occurs once for each point, and a table must have one point or more. The energy is in units of
  13.6 eV.
- The first energy must be 0, which is the threshold, and the energies must increase. `sn3d` uses the first
  cross-section at the threshold, and it does not read the energy of that point.
- The level indices of `phixsdata.txt` always start at 1.
- A cross-section must not be negative, in `phixsdata.txt` and in `phixsdata_v2.txt`.
- Each table must have the ground level of the upper ion as its only target level.
- `sn3d` interpolates each table onto a grid of 100 points with an increment of 0.1. If the run also reads
  `phixsdata_v2.txt`, the grid of that file applies. Above the last point of a table, the cross-section decreases
  with the third power of the frequency.

An ion with a table for an excited level must also have a table for its ground level.

## vpkt.txt

Required when virtual packets are enabled (`VPKT_ON` in `artisoptions.h`). The file sets the observer directions, the
number of spectra for each observer, and the ranges of time and wavelength. `artistools makevpktinput` writes the
file for a model. See
[tests/classicmode_3d_inputfiles/vpkt.txt](https://github.com/artis-mcrt/artis/blob/main/tests/classicmode_3d_inputfiles/vpkt.txt)
for an example.

`sn3d` reads the values in sequence, and each space or line break separates two values. The file can have no comment
lines.

| Sequence | Values | Meaning |
| --- | --- | --- |
| 1 | `nobsdirections` | The number of observer directions, which must be 1 or more. |
| 2 | one value for each direction | The cosine of the polar angle theta of each direction, from -1 to 1. |
| 3 | one value for each direction | The azimuthal angle phi of each direction in degrees. |
| 4 | flag, and optional values | With the flag 1, the number of spectra for each observer follows, and then one opacity value for each spectrum. The number must be 1 or more. With a different flag, each observer gets one spectrum with the full opacity. |
| 5 | `flag tmin_days tmax_days` | With the flag 1, the virtual packets must arrive between these two times in days. With a different flag, the range is from `VSPEC_TIMEMIN` to `VSPEC_TIMEMAX` of `vpkt.h`, but the two values must stay. |
| 6 | flag, and optional values | With the flag 1, the number of wavelength ranges follows, and then the lowest and the highest wavelength of each range in Angstroms. The number must be 1 or more. With a different flag, the range is from `VSPEC_NUMIN` to `VSPEC_NUMAX` of `vpkt.h`. |
| 7 | `flag optical_depth_is_thick_vpkt` | With the flag 1, a cell with a grey optical depth above this value sends no virtual packets. With a different flag, the value of `input.txt` applies, but the value must stay. With the flag 1, the value must be above 0. |
| 8 | `tau_max_vpkt` | The code stops a virtual packet when its optical depth is larger than this value in every spectrum. The value must be above 0. |
| 9 | flag | A value that is not 0 switches on the velocity grid map. |
| 10 | `tmin_grid_days tmax_grid_days` | Only with the velocity grid map: the time range of the map in days. |
| 11 | number, and the ranges | Only with the velocity grid map: the number of wavelength ranges of the map, and then the lowest and the highest wavelength of each range in Angstroms. The number must be 1 or more. |

The opacity value of a spectrum removes one part of the opacity along the path of the virtual packets:

| Value | Opacity of the spectrum |
| --- | --- |
| 0 | the full opacity. The first spectrum must have this value. |
| -1 | no line opacity |
| -2 | no bound-free opacity |
| -3 | no free-free opacity |
| -4 | no electron scattering opacity |
| a positive number Z | no line opacity of the element with the atomic number Z. A run with `VPKT_USE_BINNED_OPACITIES` cannot use this value. |

The ranges have limits. The time range of the spectra must be inside the range from `VSPEC_TIMEMIN` to
`VSPEC_TIMEMAX` and inside the time range of the model. The wavelength ranges of the spectra must be inside the
range that the frequencies `VSPEC_NUMAX` and `VSPEC_NUMIN` of `vpkt.h` give, which is 3500 to 10000 Angstroms. The ranges of the velocity grid map must be inside the ranges of the
spectra.

## Ye.txt

An optional file with the electron fraction of the cells. It replaces the values of the `Ye` column of `model.txt`.
Each value must be in the range from 0 to 1, and `sn3d` accepts a value up to 1.001 for the rounding. With
`RPKT_GREY_TYPE = TANAKA2020_ELECTRONFRAC`, each cell with matter
needs an electron fraction above 0 from this file or from that column, or the run stops. The file can have no comment lines.

```text
<number of rows>
<cell id> <electron fraction>
```

The second line occurs once for each row. The cell ids are the ids of `model.txt`. The file can hold a subset of the
cells. The number of rows can be 0, and it must not be negative.

## recombrates.txt

An optional file with recombination rate coefficients. With this file, `sn3d` scales the photoionisation
cross-sections of each ion, so that its recombination rate coefficient at the temperature
`RECOMBCALIBRATION_T_ELEC` agrees with the file. The file can have comment lines.

The calibration has two limits, and the log gives the result for each ion:

- `sn3d` makes no calibration of an ion whose recombination rate coefficient from the atomic data is not positive.
- The calibration to `rrc_low_n` scales the cross-sections of all levels of the ion. `sn3d` skips this calibration if
  the necessary factor is below 0.05, or if it is 2.0 or more. It then still makes the calibration to the total rate
  coefficient.

The file has one block for each ion that it calibrates. An ion with no block gets no calibration:

```text
<Z> <upper ion stage> <number of rows>
<log10 of T_e in K> <rrc_low_n> <rrc_total>
```

- The upper ion stage is the ion stage that recombines.
- The number of rows can be 0, and it must not be negative.
- The second line occurs once for each row. `rrc_low_n` is the rate coefficient to the low-n levels, and `rrc_total`
  is the total rate coefficient, both in cm^3/s.
- `sn3d` interpolates each rate coefficient between the two rows around `RECOMBCALIBRATION_T_ELEC`. It makes no
  calibration with a result that is negative. For a rate coefficient that is not available, give a negative value in
  both of these rows. A negative value in only one row can give a positive result.
- The rows of an ion that the run uses must include a temperature below and a temperature above
  `RECOMBCALIBRATION_T_ELEC`.

## autoion.txt

An optional file with the autoionisation rates. The file can have comment lines. It has one row for each
transition:

```text
<Z> <upper ion stage> <upper level index> <lower ion stage> <lower level index> <A value in 1/s>
```

- The lower level is the level that autoionises. The upper level is the level of the ion that the transition makes.
  The upper ion stage must be higher than the lower ion stage.
- Put the transitions of each lower level on consecutive lines.
- `sn3d` finds the first level index of `autoion.txt` separately from `adata.txt`: the numbering starts at 0 if a
  level number in the file is 0, and at 1 otherwise.
- The autoionising levels of an ion must be one block of levels with no gap, and the block must end at the highest
  level of the ion.
- The autoionising levels of the file must all be {term}`NLTE` levels, or none of them must be an NLTE level.
- The A value must not be negative.

## The gamma-ray line files

The gamma-ray lines of each nuclide are in the `data` folder of the repository. The file of a nuclide has the name
`gamma_<nuclide>.txt` with lower-case letters, e.g. `gamma_ni56.txt`. The file can have comment lines.

```text
<number of lines>
<energy in MeV> <number of photons for each decay>
```

The number of lines must not be negative. The second line occurs once for each gamma-ray line. The energy and the
number of photons must not be negative. A nuclide without a file, but with a gamma-ray energy in its
decay data, gets one line at that energy. These nuclides are an exception, and the run stops when one of them is in
the run and has no file:

- Ni56 and Co56;
- Cr48 and V48;
- Ni57 and Co57;
- Sc44.

Fe52 and Mn52 have no file and no gamma-ray energy in the decay data. Without a file, they get one line each, at
0.86 MeV and at 3.415 MeV.
