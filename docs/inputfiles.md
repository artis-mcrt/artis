# Input files

Put these files in the model folder. The model folder must also contain the ARTIS source folder, or a symlink to it,
with the name `artis`. The code then finds the physics data that comes with it in `artis/data`, and you copy nothing
(see [Atomic and model data](data.md)).

A build with libzstd also reads a text input file in zstd compressed form, e.g. `model.txt.zst`,
`transitiondata.txt.zst`, or `packets/packets00_0000.out.zst` for `exspec`. The plain file has priority when both
exist.

`sn3d` writes the final packet files of a job into the folder `packets/`, and the virtual packet files into
`vpackets/`, `vspecpol/`, and `vpkt_grid/`. `exspec` reads the packet files there.

| File | Necessary | Content |
| --- | --- | --- |
| [input.txt](#inputtxt) | yes | the run-time parameters |
| [model.txt](#modeltxt) | yes | the grid, the densities, and the nuclide mass fractions |
| [abundances.txt](#abundancestxt) | yes | the element mass fractions of each cell |
| [compositiondata.txt](#compositiondatatxt) | yes | the elements, the ion stages, and the level counts of the run |
| [adata.txt](#adatatxt) | yes | the energy levels of each ion |
| [transitiondata.txt](#transitiondatatxt) | yes | the lines of each ion |
| [phixsdata_v2.txt or phixsdata.txt](#phixsdata_v2txt-and-phixsdatatxt) | one of them | the photoionisation cross-sections |
| [vpkt.txt](#vpkttxt) | with `VPKT_ON` | the observer directions of the virtual packets |
| [Ye.txt](#yetxt) | no | the electron fraction of each cell |
| [recombrates.txt](#recombratestxt) | no | the recombination rate coefficients for the calibration |
| [autoion.txt](#autoiontxt) | no | the autoionisation rates |
| [gridcontributions.txt](#gridcontributionstxt) | no | `sn3d` does not read this file |
| [gamma_*.txt](#the-gamma-ray-line-files) | in `artis/data` | the gamma-ray lines of each nuclide |

In the layouts of this page, `<name>` is a value that you supply. Where this page says that a file can have comment
lines, `sn3d` ignores each blank line and each line whose first character that is not a space is `#`.

## input.txt

`input.txt` holds the run-time parameters. The format is positional: each line has a fixed meaning, and the unused
lines must stay. The file can have comment lines. `sn3d` ignores the text after the values of a line, so each line
can end with a comment.

The line numbers of this table count only the lines that are not comment lines, and they start at 0:

| Line | Values | Meaning |
| --- | --- | --- |
| 0 | `pre_zseed` | The seed of the random number generator. A positive value is the seed. For zero or a negative value, `sn3d` makes a random seed. Reproducible results need a fixed seed, and a build with `REPRODUCIBLE=ON` stops without a positive value. |
| 1 | `ntimesteps` | The number of timesteps of the model. |
| 2 | `timestep_start timestep_finish` | The first timestep of this job, and the stop timestep. The job does not run the stop timestep. |
| 3 | `tmin_days tmax_days` | The start time and the end time of the model in days. |
| 4 | | Unused (`nusyn_min_mev nusyn_max_mev`). |
| 5 | | Unused (`nsyn_time`). |
| 6 | | Unused (the start and end times for the synthesis). |
| 7 | | Unused (the number of dimensions of `model.txt`). `sn3d` finds this number from `model.txt`. |
| 8 | | Unused (the mode of the light curve). |
| 9 | | Unused (`n_out_it`). |
| 10 | | Unused (the factor for the speed of light). |
| 11 | | Unused (`gamma_kappagrey`). |
| 12 | | Unused (`syn_dir`). The direction is always `0 0 1`. |
| 13 | | Unused (`opacity_case`). |
| 14 | | Unused (`rho_crit_para`). |
| 15 | | Unused (`debug_packet`). |
| 16 | `simulation_continued_from_saved` | 0 starts a new simulation. 1 continues from the restart files of a previous job. With `timestep_start` 0, the run is always a new simulation. |
| 17 | | Unused (`rfcut`). |
| 18 | `num_lte_timesteps` | The number of timesteps at the start of the model that use LTE populations. The value must be 1 or more. |
| 19 | `optical_depth_is_thick num_grey_timesteps` | The condition for the grey opacity treatment of a cell (see below). |
| 20 | | Unused (`max_bf_continua`). |
| 21 | `nprocs_exspec` | The number of packet files that `exspec` reads, one for each `sn3d` rank. `sn3d` sets this line itself. |
| 22 | | Unused (`do_emission_res`). |
| 23 | | Unused (the diffusion values of the k-packets). |

Sixteen of the 24 lines are unused. See
[tests/classicmode_3d_inputfiles/input-newrun.txt](https://github.com/artis-mcrt/artis/blob/main/tests/classicmode_3d_inputfiles/input-newrun.txt)
for a complete example with a comment on every line.

The timesteps:

- `timestep_start` must be 0 or more, and it must be less than `ntimesteps`. `timestep_finish` must not be less than
  `timestep_start`, and it must not be more than `ntimesteps`.
- Give the number of timesteps as the stop value for the whole model. Long simulations run in several jobs, and the
  wall-time mechanism (`sn3d -w`) advances the start timestep at each restart and keeps the stop value.
- A smaller stop value makes the simulation stop earlier, e.g. to examine the partial results. The simulation
  continues only after you edit `input.txt`. The last restart files of such a run hold the start of its last
  timestep. The continuation job therefore propagates that timestep again and writes its spectra and light curves
  again.
- `sn3d` writes the restart files at the start of each timestep after the first timestep of the job, so give a job at
  least two timesteps. A job of one timestep writes no restart files.
- `tmin_days` must be more than 0 and less than `tmax_days`. The snapshot time of `model.txt` must not be after
  `tmin_days`.

The LTE timesteps and the grey cells:

- With a value for `MULTIBIN_RADFIELD_FIRST_TIMESTEP`, `num_lte_timesteps` must not be larger than that value.
- With the Spencer-Fano solver (`NT_SCHEME`), the first solution comes at the timestep `num_lte_timesteps` + 1.
- `sn3d` calculates the grey optical depth from each cell to the outer boundary of the model. A cell with a depth of
  `optical_depth_is_thick` or more gets the grey opacity treatment. This applies to each timestep before the timestep
  `num_grey_timesteps`.

At the start of a new simulation, `sn3d` writes `input.txt` again with a comment on each line and with the value of
`nprocs_exspec`, and it keeps the same content in `input-newrun.txt`. If `input.txt` is absent at the start of a run,
rank 0 restores it from `input-newrun.txt`.

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

- The file can have comment lines before the number of cells and before `t_model`.
- `t_model` is the snapshot time of the model. The densities and the positions of the file apply at this time.
- `vmax` is the velocity of the outer boundary of the grid. Put this line directly after the line of `t_model`. A 1D
  model has no such line, and its `vmax` is the outer velocity of the last cell. `vmax` must be less than the speed
  of light.
- The corners of a 2D or a 3D grid move faster than `vmax`. The corner velocity is `vmax` times the square root of 2
  for a cylindrical grid, and `vmax` times the square root of 3 for a Cartesian grid. `sn3d` stops if that velocity
  is not less than the speed of light, unless the run has `FORCE_SPHERICAL_ESCAPE_SURFACE`.
- The number of cells of a 3D model must be a cube number, because the grid has the same number of cells on each axis.
- The header line names the columns. It is optional. The columns of the cell id, the position, and the density must
  have the names and the sequence of the layouts above. The column `velocity_outer` is a second name for
  `vel_r_max_kmps`. `sn3d` accepts each column name that starts with `pos_` as a position column.
- A file without a header line gets the columns of the layouts above. A row with two more values gets `X_Ni57` and
  `X_Co57` as the last two columns.

The columns of each cell:

| Column | Model | Meaning |
| --- | --- | --- |
| `#inputcellid` | all | The cell id. The ids start at 0 or at 1, and they increase by one for each cell. |
| `vel_r_max_kmps` | 1D | The velocity of the outer boundary of the shell [km/s]. Each value must be larger than the value of the cell before. |
| `logrho` | 1D | The base-10 logarithm of the density [g/cm^3] at `t_model`. A value of -90 or less gives an empty cell. |
| `pos_rcyl_mid`, `pos_z_mid` | 2D | The midpoint of the cell [cm] at `t_model`. |
| `pos_x_min`, `pos_y_min`, `pos_z_min` | 3D | The corner of the cell with the lowest coordinates [cm] at `t_model`. |
| `rho` | 2D, 3D | The density [g/cm^3] at `t_model`. A value of zero gives an empty cell. |
| `X_Fegroup` | all | The mass fraction of the iron-group elements. |
| `X_<nuclide>` | all | The mass fraction of one nuclide at `t_model`, e.g. `X_Ni56`. The name has the element symbol and the mass number. |
| `Ye` or `cellYe` | all | The electron fraction. |
| `q` | all | The trapped radiation energy per mass [erg/g] at `t_model`. |
| `tracercount` | all | `sn3d` reads this column and does not use it. |

More rules for the cells:

- A cell has one line or two lines. If the first cell line has a value for each column, each cell has one line. If
  not, each cell has two lines: the first line holds the cell id, the position, and the density, and the second line
  holds the other columns.
- The cells of a 2D model go through the r index first and then through the z index. The z midpoints start at the
  negative boundary. `sn3d` calculates the midpoint of each cell from `vmax` and `t_model`, and it stops if a value
  of the file differs by more than 0.1 percent.
- The cells of a 3D model go through the x index first, then the y index, and then the z index. `sn3d` compares the
  positions with the values that it calculates from `vmax` and `t_model`. A difference gives a warning in the log.
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
model needs.

For a 1D or a 2D model, `sn3d` divides the values of each cell by their sum, so the values can also be proportional to
the mass fractions, e.g. element densities. For a 3D model, `sn3d` uses the values without normalisation. A 3D cell
whose values do not sum to 1 within 2 percent gets a warning in the log.

The mass fraction of an element must not be less than the sum of the mass fractions of its nuclides in `model.txt`.
A cell with a density, but with no mass in the elements of `compositiondata.txt`, becomes an empty cell. The log gets
a warning that names the elements of that cell.

## compositiondata.txt

`compositiondata.txt` selects the elements, the ion stages, and the levels that the code reads from the atomic data
files. The file can have comment lines, and the code ignores the text to the right of a `#` character.

```text
<number of elements>
0
0
<Z> <nions> <lowest ion stage> <highest ion stage> <nlevelsmax> <abundance> <mass>
```

The second value and the third value are unused and must be 0. Then the file has one row for each element:

| Column | Meaning |
| --- | --- |
| `Z` | The atomic number. |
| `nions` | The number of ion stages. With the value 0, the element has no ions, and the code keeps only its nuclide abundances. |
| lowest ion stage | The lowest ion stage of the run. The ion stage 1 is the neutral atom. |
| highest ion stage | The highest ion stage of the run. `nions` must be equal to the highest ion stage - the lowest ion stage + 1. |
| `nlevelsmax` | The maximum number of levels that the code reads for each ion of the element. A negative value reads every level. |
| abundance | Unused. The value must be 0 or more. |
| mass | The mean mass of the stable nuclides of the element in atomic mass units. |

## adata.txt

`adata.txt` holds the energy levels. It has one block for each ion. The blocks of the ions that the run uses must have
the sequence of `compositiondata.txt`, with the ion stages of an element in increasing order. `sn3d` goes past the
blocks of other ions.

The header of a block:

```text
<Z> <ion stage> <number of levels> <ionisation energy in eV>
```

Then one row for each level:

```text
<level index> <energy in eV> <statistical weight> <number of transitions> <level name>
```

- The energy of a level is relative to the ground level of the ion.
- `sn3d` reads the number of transitions and the level name, and it does not use them.
- A level with an energy below the energy of the level before gives a warning in the log.

The level indices start at 0 or at 1. `sn3d` finds the first index from the first level of `adata.txt`.
`transitiondata.txt` and `phixsdata_v2.txt` must use the same numbering. `phixsdata.txt` always numbers the levels
from 1. The log messages always give the 0-based level index.

## transitiondata.txt

`transitiondata.txt` holds the lines. It has one block for each ion, in the same sequence as `adata.txt`.

The header of a block:

```text
<Z> <ion stage> <number of transitions>
```

Then one row for each transition:

```text
<lower level index> <upper level index> <A value in 1/s> <collision strength> <forbidden flag>
```

- The lower level index must be less than the upper level index.
- A positive collision strength is the effective collision strength of the transition. If the collision strength is
  not available, the value is -1.0 for a permitted transition and -2.0 for a forbidden transition.
- The forbidden flag is 1 for a forbidden transition and 0 for a permitted transition.
- The rows of a block can have no comment lines between them.

The code also accepts an older format of four columns, without the collision strength and the forbidden flag:

```text
<transition index> <lower level index> <upper level index> <A value in 1/s>
```

The code finds the format of each block from the number of columns in its first row.

## phixsdata_v2.txt and phixsdata.txt

These files hold the photoionisation cross-sections. A run needs one of them, and it can read both. `phixsdata_v2.txt`
has the same frequency grid for every table, and `phixsdata.txt` has an energy table for each level.
`phixsdata.txt` can have comment lines. `phixsdata_v2.txt` can have them only before the header line of a table and
before a target line, and not before the grid values or between the cross-section values.

`phixsdata_v2.txt` starts with the grid:

```text
<NPHIXSPOINTS>
<NPHIXSNUINCREMENT>
```

Point number i of a table, from 0, has the frequency nu_edge (1 + i `NPHIXSNUINCREMENT`), where nu_edge is the
threshold frequency of the level.

Then each table of `phixsdata_v2.txt` has this layout:

```text
<Z> <upper ion stage> <upper level index> <lower ion stage> <lower level index> <threshold energy in eV>
<number of target levels>
<target level index> <probability>
<NPHIXSPOINTS cross-sections in megabarns>
```

- The lower level is the level that absorbs the photon. The upper ion stage must be the lower ion stage + 1.
- An upper level index of 0 or more gives one target level. The table then has no row with the number of target
  levels and no rows of target levels.
- A negative upper level index, e.g. -1, gives a list of target levels. The probabilities of the target levels must
  sum to 1, and each target level must occur only once.
- `sn3d` reads the threshold energy and does not use it. The threshold comes from the level energies of `adata.txt`.

Each table of `phixsdata.txt` has this layout:

```text
<Z> <upper ion stage> <upper level index> <lower ion stage> <lower level index> <number of points>
<energy above the threshold> <cross-section in megabarns>
```

- The second line occurs once for each point. The energy is in units of 13.6 eV.
- The level indices of `phixsdata.txt` always start at 1.
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
| 1 | `nobsdirections` | The number of observer directions. |
| 2 | one value for each direction | The cosine of the polar angle theta of each direction, from -1 to 1. |
| 3 | one value for each direction | The azimuthal angle phi of each direction in degrees. |
| 4 | flag, and optional values | With the flag 1, the number of spectra for each observer follows, and then one opacity value for each spectrum. With a different flag, each observer gets one spectrum with the full opacity. |
| 5 | `flag tmin_days tmax_days` | With the flag 1, the virtual packets must arrive between these two times in days. With a different flag, the range is from `VSPEC_TIMEMIN` to `VSPEC_TIMEMAX` of `vpkt.h`, but the two values must stay. |
| 6 | flag, and optional values | With the flag 1, the number of wavelength ranges follows, and then the lowest and the highest wavelength of each range in Angstroms. With a different flag, the range is from `VSPEC_NUMIN` to `VSPEC_NUMAX` of `vpkt.h`. |
| 7 | `flag optical_depth_is_thick_vpkt` | With the flag 1, a cell with a grey optical depth above this value sends no virtual packets. With a different flag, the value of `input.txt` applies, but the value must stay. |
| 8 | `tau_max_vpkt` | The code stops a virtual packet when its optical depth is larger than this value in every spectrum. |
| 9 | flag | 1 switches on the velocity grid map. |
| 10 | `tmin_grid_days tmax_grid_days` | Only with the velocity grid map: the time range of the map in days. |
| 11 | number, and the ranges | Only with the velocity grid map: the number of wavelength ranges of the map, and then the lowest and the highest wavelength of each range in Angstroms. |

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
range from `VSPEC_NUMIN` to `VSPEC_NUMAX`. The ranges of the velocity grid map must be inside the ranges of the
spectra.

## Ye.txt

An optional file with the electron fraction of the cells. It replaces the values of the `Ye` column of `model.txt`.
The file can have no comment lines.

```text
<number of rows>
<cell id> <electron fraction>
```

The second line occurs once for each row. The cell ids are the ids of `model.txt`. The file can hold a subset of the
cells.

## recombrates.txt

An optional file with recombination rate coefficients. With this file, `sn3d` scales the photoionisation
cross-sections of each ion, so that its recombination rate coefficient at the temperature
`RECOMBCALIBRATION_T_ELEC` agrees with the file. The file can have comment lines.

The file has one block for each ion:

```text
<Z> <upper ion stage> <number of rows>
<log10 of T_e in K> <rrc_low_n> <rrc_total>
```

- The upper ion stage is the ion stage that recombines.
- The second line occurs once for each row. `rrc_low_n` is the rate coefficient to the low-n levels, and `rrc_total`
  is the total rate coefficient, both in cm^3/s.
- `sn3d` interpolates each rate coefficient between the two rows around `RECOMBCALIBRATION_T_ELEC`. It makes no
  calibration with a result that is negative. Give a negative value in both of these rows for a rate coefficient that
  is not available, because a negative value in only one row can give a positive result.
- The rows of an ion that the run uses must include a temperature below and a temperature above
  `RECOMBCALIBRATION_T_ELEC`.

## autoion.txt

An optional file with the autoionisation rates. The file can have comment lines. It has one row for each
transition:

```text
<Z> <upper ion stage> <upper level index> <lower ion stage> <lower level index> <A value in 1/s>
```

- The lower level is the level that autoionises. The upper level is the level of the ion that the transition makes.
- Put the transitions of each lower level on consecutive lines.
- `sn3d` finds the first level index of `autoion.txt` separately from `adata.txt`: the numbering starts at 0 if a
  level number in the file is 0, and at 1 otherwise.
- The highest autoionising level of each ion must be the highest level of the ion.

## gridcontributions.txt

`sn3d` does not read this file. [artistools](https://github.com/artis-mcrt/artistools) writes it together with a
kilonova model, and it gives the contribution of each tracer particle to each cell. The files of the tests have the
columns `particleid`, `cellindex`, `frac_of_cellmass`, and `frac_of_cellmass_includemissing`.

## The gamma-ray line files

The gamma-ray lines of each nuclide are in the `data` folder of the repository. The file of a nuclide has the name
`gamma_<nuclide>.txt` with lower-case letters, e.g. `gamma_ni56.txt`. The code also accepts the name
`<nuclide>_lines.txt`. The code looks for each file in the model folder first, then in `data/`, and then in
`artis/data/`. The file can have comment lines.

```text
<number of lines>
<energy in MeV> <number of photons for each decay>
```

The second line occurs once for each gamma-ray line. A nuclide without a file, but with a gamma-ray energy in its
decay data, gets one line at that energy.
