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

The build needs gcc 14 or a later version. Load a newer gcc. With Open MPI you can also export `OMPI_CXX` with the name of the newer g++.

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

### The log has the line RESTART_NEEDED

This is the clean end of a job, and not an error. `sn3d` writes that text when the job stops because too little wall time remains
for another {term}`timestep`. The {term}`restart files <restart file>` and `input.txt` are then already current,
because `sn3d` writes them after the grid update of each timestep, except in the first timestep of a job. The next job continues the
simulation. The job scripts of `scripts/` search the log for the text and submit the next job. See
[Running sn3d](running.md).

### A resumed run stops at a restart file

The restart files must come from the binary that wrote them. A resumed run therefore needs the same `artisoptions.h`
and the same source version, and also the same rank count, the same timestep count, and the same time range. The
message names the value of the file and the value of the run. The code has no fallback for a file of a different
build. Start a new simulation, or build the source version of that run again.

### A job stopped with a crash or a memory kill

A crash or a kill stops the job before the last line of the log, `sn3d finished (job: ...)`. The job scripts then
submit no next job (see [Runs on a cluster](clusters.md#what-the-job-scripts-do)). Find the cause in these places:

- the output file of the Slurm job, `slurm-<job ID>.out`. It holds the standard error of `sn3d` and the messages of
  Slurm, e.g. about a kill at the memory limit.
- the logs of the ranks in the {term}`job folder`. The message that stops the run is in the log of the rank that found
  the fault (see [The tags of the log lines](#the-tags-of-the-log-lines)). This command lists the logs that have an
  error line or a failed assertion:

  ```sh
  grep -l -e '\[error\]' -e 'failed assertion' job_from_ts*/output_*.txt
  ```

- `sacct -j <job ID> --format=JobID%20,State%20,MaxRSS`. It gives the end state of each step of the job. It also
  gives the largest memory of one task of the step. The state `OUT_OF_MEMORY` is a kill at the memory limit.

`sn3d` writes the {term}`restart files <restart file>` and updates `input.txt` after the grid update of each
{term}`timestep`, except the first timestep of a job. It removes the restart files of the timestep before only after
each rank writes its new files. After a crash, `input.txt` thus names the timestep of the newest complete set of
restart files. After a crash in the first timestep of a job, `input.txt` and the restart files are the same as at the
start of that job.

Correct the cause. Then submit the job again, e.g. with `artis/scripts/submit.sh`. The next job continues from
the timestep that `input.txt` names. If that timestep is the start timestep of the crashed job, the next job uses the
same job folder and removes the logs of the crashed job. Copy these logs first if you need them.

For a memory kill, read the `mem_usage` lines of the log (see
[Plan a simulation](planning.md#measure-the-cost-of-your-run)). A resumed run needs the same rank count and the same
`NUM_PACKETS` as the run that wrote the restart files (see [Plan a simulation](planning.md#the-ranks)). Each rank
holds its share of `NUM_PACKETS`, and the ranks of a node share one copy of the large arrays. To resume, give the job
more nodes and fewer ranks on each node, with the same rank count. Change the `#SBATCH` lines of the job script in
`artis` for this. On Virgo, change `artis/scripts/artis-virgo-submit.sh`. An option on the command line of `sbatch`
changes only one job, because each job script submits the next job with its own values. A new simulation can
also use a smaller `NUM_PACKETS` or a lower level limit in `compositiondata.txt` (see
[Plan a simulation](planning.md#the-atomic-data-and-the-nlte-levels)).

### The code or the options changed between two jobs

A resumed run needs the same `artisoptions.h` and the same source version as the run that wrote the restart files.
Each job script runs `make` in `artis` before it starts `sn3d`, so a change of the clone between two jobs goes into the
next job. The checks of the restart files compare only some values (see
[A resumed run stops at a restart file](#a-resumed-run-stops-at-a-restart-file)), e.g.:

- the timestep count;
- the rank count;
- the count of the cells that contain matter;
- the count of the ions;
- the time range.

A change of the code that keeps these values can thus resume with no message.

Thus change the code or the options only between two simulations. Give each {term}`model folder` its own clone (see
[Runs on a cluster](clusters.md#make-the-model-folder)). A change for one simulation then does not reach another.
The `git version` and `git status` lines at the start of the log of each job give the code of that job.

### The run starts again from the first timestep

Keep the parameters for the saved timestep in the plain file `input.txt`.
If only `input.txt.zst` holds those parameters, a plain `input-newrun.txt` takes priority.
If the plain `input.txt` is absent, recover the current parameters with this command.

```sh
zstd -d input.txt.zst -o input.txt
```

Run the [input check](running.md#check-the-input-before-a-job) with the rank count of the saved run.
Read the selected mode and timestep in its log.

```sh
grep -E 'input: (resuming simulation|starting a new simulation)|sn3d stops before timestep' inputcheck/output_0-0.txt
```

For a restart, the log must show `resuming simulation from saved point` and the timestep of the restart files.
If it shows a new simulation, correct `input.txt` before you run `sn3d`.
A new simulation [removes the previous output and restart files](outputfiles.md#a-new-simulation-removes-the-files-of-the-previous-one).

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

An error in the command line of `sn3d` or `exspec` is different. Rank 0 writes `[error] <message>` and the usage to
the standard error, and the program stops before it reads a file.

`artistools logsummary <model folder>` counts each kind of warning in the logs of rank 0, with the first and the last
timestep of each kind. The line of each timestep in the log of rank 0 gives the number of warnings in the logs of all
ranks since the start of the job. The end of the log gives the total.

### Warnings that most runs give

The run continues after each of these warnings:

| Start of the message | Meaning |
| --- | --- |
| `Z=... ionstage ... nlevels_groundterm ... phixstargetlevels(ion-1) ...` | A check of the atomic data at the start. It compares the number of levels in the ground term of the ion with the number of target levels of the photoionisation from the ground level of the lower ion. The two numbers are different. |
| `Z=... ionstage ... shellnum ... en_auger_ev is ... Setting to zero.` | The Auger energy of a shell is negative or not finite after the correction of the data table, and the code uses zero. |
| `cell ... ts ... bin ...: W ... too high or non-finite ...; continuing with W ... at T_R_max` | The fit of one bin of the multibin radiation field gave a dilution factor that is above 1e4 or not finite. The code then uses the maximum temperature of the bins for that bin. |
| `timestep ... cell ...: Some bin T_R values were clamped.` | The fit of some bins gave a temperature at the lower or the upper limit of the bins. |
| `cell ... ts ...: pop inversion greater than factor ...` | A level has a population inversion relative to the ground level of its ion that is above [STRICT_POPULATION_CHECKING_INVERSION_FACTOR_PRINTOUT_WARNING](options.md#strict_population_checking_inversion_factor_printout_warning). The end of the message says whether the solver keeps the solution. |
| `cell ... timestep ...: B_planck(T_e=...) peak at ... is outside the radiation field frequency range` | The Planck function of the electron temperature has its peak outside the range from [NU_MIN_R](options.md#nu_min_r) to [NU_MAX_R](options.md#nu_max_r). |
| `calculate_ionfractions: ... non-finite ionfract set to zero` | An ion fraction of the ionisation balance was not finite, and the code uses zero for that ion. |
| `NLTE (Spencer-Fano/Te/pops) solver mgi ... failed to converge` | The iteration of T_e and the NLTE populations of a cell stopped at the limit [NLTE_TE_NNE_MAXITER](options.md#nlte_te_nne_maxiter) before it converged. When [NT_SCHEME](options.md#nt_scheme) is NT_SPENCERFANO, the iteration also includes the non-thermal solution. The code keeps the solution of the last iteration. |
### Warnings of the temperature limits

The run continues after these warnings. The code limits a radiation temperature to the range from
[MINTEMP](options.md#mintemp) to [MAXTEMP](options.md#maxtemp). For T_J the log gives
`get_T_J_from_J: ... Clamping to MINTEMP` or `Clamping to MAXTEMP`. For T_R it gives
`temperature estimator T_R = ... Setting T_R = T_min!` or `Setting T_R = T_max!`.

For the electron temperature, the solver tests heating minus cooling at both limits. When the values have the same
sign or either value is zero, the solver selects a limit. A solution can still lie between the limits.

- `call_T_e_finder: cooling bigger than heating at lower T_e boundary ... T_e forced to be MINTEMP` means cooling
  exceeds heating at the upper limit. The solver selects `MINTEMP`.
- `call_T_e_finder: heating bigger than cooling over the whole T_e range ... T_e forced to be MAXTEMP` means heating
  equals or exceeds cooling at the upper limit. The solver selects `MAXTEMP`.

The solver then limits the change in T_e for this iteration. The stored temperature can lie between the limits.
If the solver applies damping, the `T_e damping` message gives the temperature it stores.

`artistools check <model folder>` counts the cells at a limit at each timestep. It gives a warning when more than 10%
of the cells of a timestep are at a limit.
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

### A warning of the energy count of the Spencer-Fano solution

- `frac_sum is ..., but should be 1.0`: the heating, ionisation, and excitation fractions of the Spencer-Fano
  solution do not sum to 1 within a relative limit of 1e-5. The matrix and the analysis of the solution count each
  energy loss with the same weights, so a correct solution gives a sum of 1 to rounding. This warning thus shows that
  a term of the solver counts some energy two times or not at all. The run continues. Without this warning, the code
  uses the calculated heating fraction. With this warning, the code sets the heating fraction to 1 minus the two
  other fractions, with the limits 0 and 1. That value is the fraction of the packets that become k-packets. If the
  two other fractions sum to more than 1, the heating fraction is 0 and the sum stays above 1.

### Messages that stop a run

Each of these messages stops the run. The page in the last column gives more information about the cause.

| Start of the message | Cause | Page |
| --- | --- | --- |
| `Could not open file '...'. The search tried ... in the folder '...'.` | `sn3d` cannot read an input file in the model folder or in the data folders `data/` and `artis/data/`. The message names each path that `sn3d` tried, and the working folder. The file is absent, or it has no read permission. | [Input files](inputfiles.md) |
| `Could not open file '...' for mode '...'.` | The mode `rb` says that a resumed run cannot read a restart file of its first timestep, e.g. `gridsave_ts<timestep>.tmp`. The file is absent, or it has no read permission. The mode `wb` says that `sn3d` cannot write a restart file, e.g. because the file system is full or the model folder has no write permission. For a packet restart file, `sn3d` first writes this text with the tag `[error]` and no source line, and tries again. The next row gives the message that then stops the run. | [Output files](outputfiles.md#the-files-of-the-whole-simulation) |
| `The write of ... failed after ... attempts` | `sn3d` cannot write the packet restart file of a rank, e.g. because the file system is full. | [Output files](outputfiles.md#the-files-of-the-whole-simulation) |
| `Found '...', but this build has no libzstd.` | The input file is zstd compressed, and the build has no libzstd. | [Make options](makeoptions.md) |
| `... ends before line ...` | `input.txt` has fewer lines than the code reads. | [Input files](inputfiles.md#inputtxt) |
| `... has the unreadable text '...' at line ...` | A line of `input.txt` has a value of the wrong type. | [Input files](inputfiles.md#inputtxt) |
| `input.txt gives tmin_days ... and tmax_days ..., but they must obey 0 < tmin_days < tmax_days` | The time range of `input.txt` is not valid. | [Input files](inputfiles.md#inputtxt) |
| `The model snapshot time ... d is after tmin ... d` | `tmin_days` is before the snapshot time of `model.txt`. | [Make your own model](ownmodel.md#4-select-the-options-and-the-times) |
| `model.txt: found only ... cells instead of ... expected.` | `model.txt` has fewer cells than its first line gives. | [Input files](inputfiles.md#modeltxt) |
| `model.txt cell ... has vel_r_max_kmps ..., but it must be larger than ... of the cell before` | The outer velocities of a 1D model do not increase. | [Input files](inputfiles.md#modeltxt) |
| `abundances.txt ends before cell ... of the ... cells of model.txt` | `abundances.txt` has fewer cells than `model.txt`. | [Input files](inputfiles.md#abundancestxt) |
| `abundances.txt cell ... has element mass fractions that sum to ..., which is more than 1` | The element mass fractions of a cell sum to more than 1. | [Input files](inputfiles.md#abundancestxt) |
| `compositiondata.txt has ... elements, but the row of element ... is absent or unreadable` | `compositiondata.txt` has fewer element rows than its first line gives, or a row has a value of the wrong type. | [Input files](inputfiles.md#compositiondatatxt) |
| `... has ... timesteps and ... ranks, but this simulation has ... timesteps and ... ranks` | A resumed run has a different timestep count or rank count than the run that wrote the restart files. | [A resumed run stops at a restart file](#a-resumed-run-stops-at-a-restart-file) |
| `the restart files have ... ranks (nprocs_exspec in input.txt), but this simulation has ... ranks` | A resumed run has a different rank count than the run that wrote the restart files. The line `nprocs_exspec` of `input.txt` gives the rank count of that run. | [A resumed run stops at a restart file](#a-resumed-run-stops-at-a-restart-file) |
| `NUM_PACKETS = ... is less than the ... ranks, but every rank needs a packet` | The run has more ranks than packets. | [Plan a simulation](planning.md#the-ranks) |
| `invalid wall time hours '...' given with -w option` | The value of `-w` is not a positive number. | [Running sn3d](running.md#a-simulation-in-several-jobs) |
