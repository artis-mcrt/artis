# Troubleshooting

This page gives the messages that stop a build or a run, and what each one needs.

## The build

### make cannot find artisoptions.h

The build needs an `artisoptions.h` in the repository root, and git does not track that file. Symlink a {term}`preset` to it:

```sh
ln -s artisoptions_classic.h artisoptions.h
```

[Option presets](presets.md) says what each preset is for.

### The build folder is large

Each set of make options and each content of `artisoptions.h` has its own folder in `build/` (see
[Make options](makeoptions.md#the-build-folders)). `make clean` removes all build folders and the program symlinks.

### make cannot run 'mpicxx --version'

The build needs an MPI library that supplies `mpicxx`. Install one, e.g. Open MPI, or load the module of your cluster.
The job scripts in `scripts/` name the modules of some clusters. With Open MPI, select the compiler with an
environment variable, e.g. `export OMPI_CXX=g++-16`.

### gcc is older than the minimum version

The build needs gcc 14 or a later version. Load a newer gcc. With Open MPI you can also point `OMPI_CXX` at it.

### The binary does not hold your change

The `sn3d` symlink of the repository root points at the most recent build of any configuration, e.g. at the build of a
test after `tests/runtest.sh`. Look at the symlink before you copy a binary, and copy it with `cp -L` from the
repository root (see [Make options](makeoptions.md#the-build-folders)).

## The run

### The photoionisation data needs phixsdata.txt or phixsdata_v2.txt

The {term}`model folder` holds no photoionisation file. [Atomic and model data](data.md) gives the datasets, and each dataset holds such a
file. A run also needs `adata.txt`, `transitiondata.txt`, and `compositiondata.txt`.

### Not enough slots available

Open MPI stops when the rank count is higher than the core count of the machine. Give `mpirun` the option
`--oversubscribe`, as `tests/runtest.sh` does.

### The log ends with RESTART_NEEDED

This is the clean end of a job, and not an error. `sn3d` writes that text when too little wall time remains for
another {term}`timestep`. It then writes the {term}`restart files <restart file>` and updates `input.txt`, so that the next job continues the
simulation. The job scripts of `scripts/` search the log for the text and submit the next job. See
[Running sn3d](running.md).

### A resumed run stops at a restart file

The restart files must come from the binary that wrote them. A resumed run therefore needs the same `artisoptions.h`
and the same source version, and also the same rank count, the same timestep count, and the same time range. The
message names the value of the file and the value of the run. The code has no fallback for a file of a different
build. Start a new simulation, or build the source version of that run again.

### The run starts again from the first timestep

`input.txt` holds the timesteps of the next job. A new simulation writes the file again and keeps a copy in
`input-newrun.txt`, and rank 0 restores `input.txt` from that copy when the file is absent. Keep `input.txt` after a
job, because it continues the simulation.

## The results

### artistools draws fewer times than the run holds

The light of the far side of the ejecta reaches the observer later than the light of the near side. artistools
therefore accepts only the observer times that have light from the whole model. The fastest ejecta set that range, so
it is shorter than the time range of the run. The row `arrival times` of `artistools check` gives it, and
`--plotinvalidpart` draws the other times.

## The log

[Running sn3d](running.md#the-logs) gives the names of the log files. Each line starts with the time in Coordinated
Universal Time (UTC). A {term}`rank` writes the messages of the cells that it updates, so the message about a cell can
be in the log of a rank other than rank 0.

### The tags of the log lines

Some lines have a tag after the time:

- `[info]` gives a fact of the setup, e.g. the memory of an array (`mem_usage`), the levels and the lines of each
  ion, or the number of NLTE levels. It needs no action.
- `[warning]` says that the code limited or replaced a value, or that an input file fails a check. The run
  continues.
- `[error]` gives the values of a fault. An assertion stops the run directly after most of these lines.

A message that stops the run has the rank, the source file, and the line. `sn3d` writes it to the log of the rank
and to the standard error, which a Slurm job puts into its output file. The message has one of these two forms:

```text
[rank 0] [error] <file>:<line> in <function>: <message>
[rank 0] <file>:<line>: failed assertion `<condition>` in function <function>
```

`artistools logsummary <model folder>` counts each kind of warning in the logs of rank 0, with the first and the last
timestep of each kind.

### Warnings that most runs give

The run continues after each of these warnings:

| Start of the message | Meaning |
| --- | --- |
| `Z=... ionstage ... nlevels_groundterm ... phixstargetlevels(ion-1) ...` | A check of the atomic data at the start. It compares the number of levels in the ground term of the ion with the number of target levels of the photoionisation from the ground level of the lower ion. The two numbers are different. |
| `Z=... ionstage ... shellnum ... en_auger_ev is ... Setting to zero.` | The Auger energy of a shell is negative or not finite after the correction of the data table, and the code uses zero. |
| `cell ... ts ... bin ...: W ... too high or non-finite ...; continuing with W ... at T_R_max` | The fit of one bin of the multibin radiation field gave a dilution factor that is above 1e4 or not finite. The code then uses the maximum temperature of the bins for that bin. |
| `timestep ... cell ...: Some bin T_R values were clamped.` | The fit of some bins gave a temperature at the lower or the upper limit of the bins. |
| `cell ... ts ...: pop inversion greater than factor ...` | A level has a population inversion relative to the ground level of its ion that is above [STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING](options.md#strict_population_checking_inversion_factor_printout_warning). The end of the message says whether the solver keeps the solution. |
| `frac_sum is ..., but should be 1.0` | The three fractions of the Spencer-Fano solution do not sum to 1. The code sets the heating fraction to 1 minus the two other fractions, with the limits 0 and 1. If the two other fractions sum to more than 1, the heating fraction is 0 and the sum stays above 1. |
| `cell ... timestep ...: B_planck(T_e=...) peak at ... is outside the radiation field frequency range` | The Planck function of the electron temperature has its peak outside the range from [NU_MIN_R](options.md#nu_min_r) to [NU_MAX_R](options.md#nu_max_r). |
| `calculate_ionfractions: ... non-finite ionfract set to zero` | An ion fraction of the ionisation balance was not finite, and the code uses zero for that ion. |
### Warnings of the temperature limits

The run also continues after the warnings below, but each of these warnings says that a
temperature is at [MINTEMP](options.md#mintemp) or at [MAXTEMP](options.md#maxtemp), and not at its solution:

- `get_T_J_from_J: ... Clamping to MINTEMP` or `Clamping to MAXTEMP`, and `temperature estimator T_R = ... Setting T_R
  = T_min!`: a radiation temperature is outside the range of the two limits.
- `call_T_e_finder: cooling bigger than heating at lower T_e boundary ... T_e forced to be MINTEMP`: the cell has no
  thermal balance in the range, because the cooling is larger than the heating at both limits.
- `call_T_e_finder: heating bigger than cooling over the whole T_e range ... T_e forced to be MAXTEMP`: the cell has no
  thermal balance in the range, because the heating is larger than the cooling at both limits.

An electron temperature at a limit is not physical. `artistools check <model folder>` counts the cells at a limit at
each timestep, and it gives a warning when more than 10% of the cells of a timestep are at a limit.
[Example: a nebular spectrum](example.md) shows a cell without enough cooling, because its elements have no atomic
data.

### Warnings that need a check of the input

- `... cells with no mass in the elements of compositiondata.txt became empty`: the cells have a density, but all of
  their mass is in elements without atomic data. `sn3d` removes the mass of these cells. The log names the atomic
  numbers of the elements, so add atomic data for them to keep the mass.
- `ignored model.txt columns not recognised ...`: `sn3d` does not use a column of `model.txt`. Check the name of the
  column.
- `adata.txt: ... level ... has energy ... below the previous level energy` and `transitiondata.txt: ... no transition
  has the ground state as its lower level`: check the order and the numbering of the levels in the atomic data files
  (see [Input files](inputfiles.md)).
- `Can't solve for NLTE populations in cell ... Attempting to use LTE solution instead`: the NLTE solution of one
  element failed in one cell, and that element gets LTE populations there. Examine the populations of the cell when
  the warning occurs at many timesteps.
- `timestep ... has a width of ... times its midpoint time`: the time-dependent ionisation and thermal balance need a
  ratio of 0.1 or less. The warning gives the first timestep above that ratio. Increase the timestep count, unless
  that timestep is in the constant segment of a hybrid scheme of
  [TIMESTEP_SIZE_METHOD](options.md#timestep_size_method). For such a timestep, decrease
  [TIMESTEP_SIZE_CONSTANT_MAX_WIDTH_DAYS](options.md#timestep_size_constant_max_width_days). The constant segment
  then takes more of the timestep count, so increase the timestep count also.
