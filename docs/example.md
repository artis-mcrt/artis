# Example: a nebular spectrum

This page takes a small model through a complete non-LTE (NLTE) run in the nebular phase. It makes the model folder,
limits the levels of the atomic data, sets the options and the timesteps, follows the run, checks the results, and
plots a spectrum. Do the [Quickstart](quickstart.md) first, because it installs the compiler, the Message Passing
Interface (MPI) library, and artistools.

The simulation of this page took 75 seconds with four ranks on a laptop (Apple M4).

:::{note}
The example uses the model and the atomic data of the `nebular_1d_3dgrid` test, so that every step runs on a laptop.
That dataset holds only iron, cobalt, and nickel. Section 8 shows the effect of this limit on the results, and section
10 says what a production run changes.
:::

## 1. Make the model folder

Each simulation has its own model folder. Clone the source code into the subfolder `artis`, so that the options of
section 4 belong to this simulation only:

```sh
mkdir ~/nebular_w7
cd ~/nebular_w7
git clone https://github.com/artis-mcrt/artis.git
```

Copy the input files of the test into the model folder. Then download the atomic data of the test and unpack it:

```sh
rsync -a --exclude 'results_md5_*' --exclude 'input-resume.txt' artis/tests/nebular_1d_3dgrid_inputfiles/ ./
curl -fLO "https://github.com/artis-mcrt/artis/releases/download/$(cat artis/tests/atomic_data_release.txt)/atomicdata_feconi.tar.xz"
tar -xf atomicdata_feconi.tar.xz
cp input-newrun.txt input.txt
```

`artis/tests/atomic_data_release.txt` names the release that holds the archive. The model folder now holds these
files:

- `model.txt.zst` and `abundances.txt.zst`: the model of the ejecta;
- `adata.txt`, `transitiondata.txt`, `phixsdata_v2.txt`, and `compositiondata.txt`: the atomic data;
- `recombrates.txt.zst`: the recombination rates for the calibration of each ion;
- `input.txt`: the run-time parameters.

A build with libzstd reads the `.zst` files directly. [Input files](inputfiles.md) describes each file.

## 2. Look at the model

```sh
artistools inputmodel describe .
```

The command gives the grid, the mass, and the mass of each element. This model is 1D and has three spherical shells
with an outer velocity of 12415 km/s. The total mass is 0.635 solar masses, and the iron group has 49.67% of it. The
comment line of `model.txt` gives the source of the cells: cells 70, 100, and 120 of the W7 model.

| Cell | Outer velocity [km/s] | Largest mass fractions |
| --- | --- | --- |
| 0 | 7075 | Ni 0.975 |
| 1 | 10260 | Si 0.496, S 0.284, Ni 0.066, Fe 0.060 |
| 2 | 12415 | Si 0.467, S 0.253, O 0.190 |

The cell identifiers of `model.txt` start at 1. The log messages and artistools give the 0-based cell index, as this table
does.

## 3. Limit the levels of the atomic data

`compositiondata.txt` selects the elements, the ion stages, and the levels that the run reads from the atomic data:

```text
3
0
0
26  5  1  5  500 0.0 55.8450
27  3  2  4  500 0.0 58.9332
28  4  2  5  500 0.0 58.6934
```

The first line is the number of elements. The next two lines must be 0. Each later line is one element, with these
columns:

1. the atomic number;
2. the number of ion stages;
3. the lowest ion stage, where 1 is the neutral atom;
4. the highest ion stage;
5. the maximum number of levels of each ion, or -1 for all the levels of `adata.txt`;
6. a value that the code does not use, but that must be zero or positive;
7. the mean atomic mass of the element in atomic mass units.

The file of this dataset thus selects Fe I to Fe V, Co II to Co IV, and Ni II to Ni V, with at most 500 levels for
each ion. `adata.txt` holds up to 2823 levels for an ion (Fe II). A run reads only the lowest levels up to the limit,
and only the lines between these levels.

A lower limit gives a quicker run with less memory, and it also changes the results. The run then has fewer lines,
and it has no photoionisation from the levels above the limit. The run stops at the start if a kept level
photoionises to a level of the upper ion that is above the limit. This table gives the memory lines of the log (see section 7) and the wall time for three limits. Each row is a
run of this model with four ranks, 20 timesteps, and 2e6 packets.

| Column 5 | Line list | Transition lists | Cell cache | Wall time |
| --- | --- | --- | --- | --- |
| -1 | 81.1 MB | 106.5 MB | 257.4 MB | 33 s |
| 500 | 15.8 MB | 20.8 MB | 52.0 MB | 28 s |
| 150 | 3.0 MB | 4.0 MB | 10.4 MB | 8 s |

Keep the value 500 for this example. To change the limit, write the new number in column 5 of each element.

## 4. Select the options

Symlink the NLTE nebular preset:

```sh
cd artis
ln -s artisoptions_nltenebular.h artisoptions.h
cd ..
```

[Option presets](presets.md) describes the preset, and [Methods and references](methods.md) gives the method and the
reference of each part. The preset has a value of `NUM_PACKETS` for about 1000 MPI ranks. Change it in
`artis/artisoptions.h`:

```c++
constexpr std::int64_t NUM_PACKETS = 1e7;
```

Each of the four ranks then holds 2.5e6 packets, which take 572 MB of memory for each rank.

Two values of the preset are important for the timesteps of the next section:

- [MULTIBIN_RADFIELD_FIRST_TIMESTEP](options.md#multibin_radfield_first_timestep) is 12. From that timestep, the
  rates use the fit of each frequency bin of the radiation field.
- [DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP](options.md#detailed_bf_estimators_first_timestep) is 13. From that
  timestep, the photoionisation rates use the detailed bound-free estimators.

## 5. Set the timesteps

`input.txt` comes from the test, which runs 8 of 10 timesteps. Change three lines, so that the run has 50 timesteps
from 170 to 250 days and runs all of them:

```text
50                       #  1: ntimesteps: number of timesteps
000 050                  #  2: timestep_start timestep_finish: ...
170 250                  #  3: tmin_days tmax_days: start and end times [day]
```

Keep the other lines. Two of them set the first timesteps:

```text
4                        # 18: num_lte_timesteps
0.0 4                    # 19: optical_depth_is_thick num_grey_timesteps
```

The first four timesteps solve the gas state in LTE, and every cell is thick in these timesteps, because the limit
of the optical depth is 0. [Plan a simulation](planning.md) describes each of these values.

## 6. Build and run

```sh
cd artis && make && cd ..
mpirun -np 4 artis/sn3d &
tail -F output_0-0.txt
```

`sn3d` writes the name of the job folder to the standard output. `tail -F` then shows the log of rank 0. Press
Ctrl+C to stop following it.

## 7. Follow the run

The log of rank 0 gives the progress after each timestep:

```sh
grep "is complete" output_0-0.txt | tail -n 1
```

```text
timestep 23 is complete (24 of 50 timesteps of the model) after 0.00 hours. If each timestep takes this time, the remaining 26 requested timesteps take 0.01 hours.
```

The `mem_usage` lines give the memory of the packets and of each large array:

```sh
grep mem_usage output_0-0.txt
```

```text
[info] mem_usage: packets occupy 572.205 MB
[info] mem_usage: linelist occupies 15.848 MB (node shared memory)
[info] mem_usage: transition lists occupy 20.800 MB (node shared memory)
...
```

The last lines of the log give the end state and the cost of the job:

```text
No need for restart
sn3d finished (job: pktprop ts 0 to ts 49 grid-preprop, 0.021 wallclock hours * 4 processes * 1 threads = 0.083 core hours)
```

`artistools corehours .` gives the sum of the core hours of all the jobs of a model folder.

## 8. Check the results

```sh
artistools check .
```

The command gives one row for each check. For this run, one row has a warning:

```text
"temperature limits" warning "... over all 50 time steps, 66.7% of the values of Te are at a limit ..."
```

A cell with an electron temperature T_e at `MINTEMP` or at `MAXTEMP` has no thermal balance, so that T_e is not
physical. Plot T_e of the first cell against the timestep:

```sh
artistools plotestimators Te -x timestep --cell 0 .
```

:::{figure} images/nebular_te.png
:alt: The electron temperature of cell 0 against the timestep

The electron temperature of cell 0 at each timestep.
:::

The curve shows the stages of the run:

- The timesteps 0 to 3 are the LTE timesteps. T_e is at `MINTEMP` in the timesteps 1 to 4. The Spencer-Fano solver
  starts at timestep 5, which is `num_lte_timesteps` + 1.
- T_e changes at timestep 13, when the detailed bound-free estimators start.
- From approximately timestep 23, T_e stays between 8230 K and 8400 K.

The gas state of a timestep comes from the estimators of the timestep before it, so a change of the method needs some
timesteps before the solution is stable. Use the results after that point only.

Cells 1 and 2 are different. The same command with `--cell 1` and with `--cell 2` shows that cell 1 is at `MAXTEMP`
from timestep 8 to the end, and cell 2 from timestep 5. The logs of their ranks have this warning:

```text
[warning] call_T_e_finder: heating bigger than cooling over the whole T_e range [1000,30000] in modelcell 1 (...). T_e forced to be MAXTEMP
```

These two cells hold mostly silicon and sulphur, and the dataset of this example has no atomic data for these
elements. The thermal balance of these cells thus has the cooling of iron, cobalt, and nickel only. The results of
these two cells are not physical. Cell 0 holds the nickel, and at 220 days it gets 95% of the deposited energy of the
model.

## 9. Plot the spectrum

```sh
artistools plotspectra -t 220 --emissionabsorption .
```

:::{figure} images/nebular_emission.png
:alt: The spectrum at 220 days with the emission of each ion

The spectrum at 220 days. The colours above the axis give the emission of each ion, and the colours below the axis
give the absorption of each ion.
:::

Fe III gives the strongest feature, near 4700 Angstroms. Without `--emissionabsorption`, the command draws only the
spectrum.

The row `arrival times` of `artistools check` gives the times that have light from the whole model. For this run,
they are 177.04 to 232.07 days, although the run goes from 170 to 250 days. artistools plots only these times. The
time of 220 days is in that range, and it is after the timestep at which T_e of cell 0 becomes stable.

[Output files](outputfiles.md) lists more artistools commands.

## 10. What a production run changes

- **The model.** A production model has the full structure of the ejecta. [Atomic and model data](data.md) says how
  artistools makes a model.
- **The atomic data.** A production dataset holds each element that has a significant mass in a cell, so that each
  cell has its cooling lines. [Atomic and model data](data.md) names the dataset for the nebular phase of a Type Ia
  supernova.
- **The packets and the ranks.** The noise of a spectrum decreases with the number of packets. Section 4 changed
  `NUM_PACKETS` in the preset file to 1e7, so set the value for your run again.
  In the model folder, `git -C artis restore artisoptions_nltenebular.h` gives the value of the repository, which is
  for approximately 1000 MPI ranks on a cluster (see [Runs on a cluster](clusters.md)).
- **The timesteps.** [Plan a simulation](planning.md) says what the timestep count, the time range, and the first
  timesteps set, and how you measure the cost of a run.
- **The checks.** Examine each warning of `artistools check`. Make sure that T_e of each cell is between the limits
  and is stable at the times of the spectra that you use.
