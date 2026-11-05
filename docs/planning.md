# Plan a simulation

A small number of values set the cost, the memory, and the noise of a run. This page says what each value does in
the code, and how you measure the cost of your own run. It gives no recommended values, because they depend on the
model. [Example: a nebular spectrum](example.md) shows the effect of some of these values on a small model.

## The values

| Value | Where | What it sets |
| --- | --- | --- |
| `NUM_PACKETS` | `artisoptions.h` | the memory of the packets, the work of the propagation, and the noise |
| `ntimesteps`, `tmin_days`, `tmax_days` | `input.txt`, lines `# 1:` and `# 3:` | the time grid |
| `timestep_start`, `timestep_finish` | `input.txt`, line `# 2:` | the timesteps of this job |
| `num_lte_timesteps` | `input.txt`, line `# 18:` | the timesteps that solve the gas state in LTE |
| `optical_depth_is_thick`, `num_grey_timesteps` | `input.txt`, line `# 19:` | the grey treatment of the first timesteps |
| the rank count | `mpirun -np` or the job script | the share of the packets and of the cells of each rank |
| the level limits | `compositiondata.txt` and `ION_NLEVELS_EXCITED_NLTE` | the memory of the atomic data and the size of the NLTE solution |
| the grid | `model.txt` and `CUBOID_NCOORDGRID_X`, `_Y`, `_Z` | the number of cells |

[Input files](inputfiles.md) describes `input.txt`, and [Compile-time options](options.md) describes each option.

## The packets

Each rank gets an equal share of [NUM_PACKETS](options.md#num_packets), and it keeps that share in memory for the
whole run. The log of each rank gives the number and the memory of its packets:

```text
Simulation propagates 2500000 packets on this rank (total 1e+07 with nprocs 4)
[info] mem_usage: packets occupy 572.205 MB
```

These lines come from the run of [Example: a nebular spectrum](example.md), where one packet takes 240 bytes. The
size depends on the build, e.g. a GPU build adds the state of the random generator to each packet, so read the value
from your log.

At the start of each timestep after the first timestep of a job, each rank writes its packets to a restart file,
`packets_nnnn_ts*.tmp`. The file of a rank has approximately the size of its packets in memory. `sn3d` removes the
files of the previous timestep only after all the ranks have written the new files, so the model folder holds two
sets for a short time.

The spectrum of a timestep holds the escaped packets that arrive at the observer in that timestep, so their number
sets the noise of the spectrum. The relative noise of a sum of N packets of equal energy is 1/sqrt(N). The log of
each rank gives the packets that escaped during the propagation of a timestep. This count is only a guide, because
the light travel time moves a packet to the timestep of its arrival:

```text
During timestep 33 on MPI process 0, 47523 pellets decayed and 49411 packets escaped. (t=220.124 [d])
```

`artistools check --frompackets <model folder>` adds a check of the packet noise.

## The time grid

`tmin_days` and `tmax_days` give the time range, and `ntimesteps` divides it. With the default
[TIMESTEP_SIZE_METHOD](options.md#timestep_size_method), the timesteps are equal in the logarithm of the time, so
each timestep is wider than the timestep before it by the same factor. The other schemes use a constant width for
all or for a part of the range. `timesteps.out` lists the start, the midpoint, and the width of each timestep.

Each timestep has one grid update, one propagation of all the packets, and one set of output (see
[How a simulation works](simulation.md)). Thus the timestep count sets these properties:

- the number of grid updates and of writes of the restart files;
- the time resolution of `spec.out` and of `light_curve.out`, because `sn3d` writes one spectrum for each timestep;
- the number of iterations of the gas state. The grid update of a timestep uses the estimators of the timestep
  before it, so the gas state needs some timesteps to follow a change.

With [NLTE_TIME_DEPENDENT_FIRST_TIMESTEP](options.md#nlte_time_dependent_first_timestep), the error of the time term
is first order in the ratio of the width of a timestep to its midpoint time. `sn3d` writes a warning to the log when
that ratio is above 0.1.

### The time range

The time range also sets how the packets divide over the time:

- With [UNIFORM_PELLET_ENERGIES](options.md#uniform_pellet_energies) true, which every preset has, each pellet
  holds the same energy. Its decay time is a sample of the decay curve between the first and the last time of the
  simulation. With the value false, the decay times are uniform, and the energy of each pellet follows the decay
  power at its time. In both modes, a longer time range gives each interval of time a smaller share of the packets.
- [INITIAL_PACKETS_ON](options.md#initial_packets_on) gives a share of the packets to the energy that the ejecta
  hold at `tmin_days`. With this option, the decay times of the pellets start at the snapshot time of `model.txt`,
  and not at `tmin_days`. A pellet that decays before `tmin_days` becomes a packet of that initial energy.

The light of the far side of the ejecta reaches the observer later than the light of the near side. artistools
therefore accepts only the observer times that have light from the whole model, and that range is shorter than the
time range of the run. The row `arrival times` of `artistools check` gives it. For the model of
[Example: a nebular spectrum](example.md), which has a maximum velocity of 12415 km/s, a run from 170 to 250 days
gives the range from 177.04 to 232.07 days.

### The first timesteps

The first timesteps of a run use approximations:

- The first `num_lte_timesteps` timesteps solve the ionisation balance in LTE. The value must be 1 or more. With the
  Spencer-Fano solver, the first solution is at timestep `num_lte_timesteps` + 1.
- In the first `num_grey_timesteps` timesteps, a cell is thick when its grey optical depth to the outer boundary is
  at or above `optical_depth_is_thick`. A thick cell gets a grey opacity, and its k-packets emit from a Planck
  function (see [Methods and references](methods.md#the-grey-opacity-of-a-thick-cell)). The log of the rank names
  each thick cell at each timestep.
- [MULTIBIN_RADFIELD_FIRST_TIMESTEP](options.md#multibin_radfield_first_timestep) and
  [DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP](options.md#detailed_bf_estimators_first_timestep) give the first timestep
  of two methods. When `MULTIBIN_RADFIELD_FIRST_TIMESTEP` has a value, `num_lte_timesteps` must be at or below it,
  or the run stops.

The gas state changes after each of these points. The example page shows the electron temperature of a cell against
the timestep, and how many timesteps the solution needs to become stable there. Make the same plot for your model.

## The ranks

Each rank does two types of work in each timestep:

- **The grid update.** `sn3d` divides the model grid cells that contain matter equally among the ranks, in the order
  of the cells. Each rank solves the state of its own cells. `modelgridrankassignments.out` and the log of each
  rank give the cells of that rank. When the rank count is at or above the number of model grid cells, each of the
  first ranks gets one cell, and the other ranks update no cell.
- **The propagation.** Each rank moves its own share of the packets through the whole grid.

The ranks wait for each other after the grid update and after the propagation. The slowest rank therefore sets the
time of a timestep. `artistools plotlogfiles <model folder>` plots the time that each rank took.

The memory of a node has two parts. The `mem_usage` lines with the text `(node shared memory)` are arrays that the
ranks of one node share, so a node holds one copy of them. Each rank holds its own copy of the other arrays and of
its packets.

A resumed run needs the same rank count as the run that wrote the restart files. `NUM_PACKETS` must be at or above
the rank count. [Running sn3d](running.md) describes the threads.

## The grid

Each model grid cell that contains matter needs one solution of the gas state at each timestep. The packets move
through the propagation cells. With [GRID_TYPE_OVERRIDE](options.md#grid_type_override), a 1D or a 2D model goes on a
3D Cartesian grid, and [CUBOID_NCOORDGRID_X](options.md#cuboid_ncoordgrid_x) and the two other values set the size of
that grid. The number of solutions of the gas state stays the same, because more than one propagation cell then
shares one model grid cell. A model grid cell that holds the centre of no propagation cell is the exception. It
becomes an empty cell and gets no solution, so use a grid that is fine enough for the thinnest cell of the model.

## The atomic data and the NLTE levels

Column 5 of `compositiondata.txt` limits the levels that a run reads for each ion (see
[Example: a nebular spectrum](example.md)). The limit sets the number of lines, and thus the memory of the line list
and of the transition lists. The example page gives the memory and the wall time for three limits. The log lists the
levels, the lines, and the bound-free transitions of each ion after it reads the atomic data:

```text
[info]  in total 12 ions, 5263 levels (4262 ionising), 519292 lines, 14780 photoionisation transitions
```

[ION_NLEVELS_EXCITED_NLTE](options.md#ion_nlevels_excited_nlte) sets the number of levels of each ion in full NLTE.
The NLTE solver puts all these levels of all the ions of an element into one rate matrix, so the value sets the
size of that matrix. The log gives the sum over all ions:

```text
[info] Total NLTE levels: 1008, of which 11 are superlevels
```

## Measure the cost of your run

Measure the cost with a short run before you ask for the full allocation:

1. Set `timestep_finish` to a small value, so that the job stops early. Give a stop value that is after the first
   timesteps of the section above. A timestep in LTE solves no NLTE populations and no Spencer-Fano equation, so it
   does not show the time of a later timestep.
2. Use the same number of ranks on each node as the full run. The memory of a node depends on it. Each rank holds
   `NUM_PACKETS` divided by the total rank count, so a run with fewer ranks holds more packets in each rank. Use the
   total rank count of the full run, or scale the `mem_usage` line of the packets.
3. Run `sn3d`, and read the log of rank 0.

The log gives these values:

- The memory of the packets and of each large array:

    ```sh
    grep mem_usage output_0-0.txt
    ```

- The time of each timestep, and an estimate for the remaining requested timesteps from the time of the last
  timestep:

    ```text
    timestep 23 is complete (24 of 50 timesteps of the model) after 0.00 hours. If each timestep takes this time, the remaining 26 requested timesteps take 0.01 hours.
    ```

- The core hours of the job, in the last line of the log:

    ```text
    sn3d finished (job: pktprop ts 0 to ts 49 grid-preprop, 0.021 wallclock hours * 4 processes * 1 threads = 0.083 core hours)
    ```

`artistools corehours <model folder>` sums the core hours of all the jobs of a model folder.

To continue after the short run, set `timestep_finish` in `input.txt` to the number of timesteps, and start `sn3d`
again. `input.txt` already holds the start timestep and the continuation flag. The new job propagates the last
timestep of the short run again (see [Input files](inputfiles.md)).

For a run in several jobs, `sn3d -w WALLTIMELIMITHOURS` tests the remaining wall time at each timestep (see
[Running sn3d](running.md)).
A job then needs wall time for at least two timesteps, because `sn3d` writes the restart files from the second
timestep of a job.
