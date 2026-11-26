# Running sn3d

Run `sn3d` in the {term}`model folder`. A local run with eight MPI {term}`ranks <rank>` looks like this:

```sh
mpirun -np 8 artis/sn3d &
tail -F output_0-0.txt
```

`tail -F` waits until `sn3d` makes the log file. Press Ctrl+C to stop following it.

## The logs

`sn3d` writes one line with the name of the {term}`job folder` to the standard output, and then nothing more unless a crash
occurs. Each MPI rank n writes a log file `output_n-t.txt` into the job folder, where t is the thread of that rank (see
[Output files](outputfiles.md)). A run without threads writes only the files of thread 0. `sn3d` keeps a symlink
`output_0-0.txt` in the model folder that points to the log of rank 0 and thread 0.

The log of rank 0 starts with these lines:

- the git version;
- the {term}`preset` that the symlink `artisoptions.h` points to;
- the compiler and the compile flags;
- the MPI library;
- the options of `artisoptions.h`.

After each {term}`timestep`, the log gives the progress of the model and an estimate of the
time for the remaining requested timesteps. The estimate uses the duration of the last timestep. The same line gives
the number of warnings in the logs of all ranks since the start of the job, and the end of the log gives the total.
Each rank writes the warnings of its own cells, so most warnings are in the logs of the other ranks and not in the log
of rank 0.

## Check the input before a job

Run `sn3d` with the option `-c` before you submit a job. The check reads the input files as a run does, and it stops
before the first timestep:

```sh
mpirun -np 4 artis/sn3d -c
```

An error in an input file then stops the check, and not a job that waited in the queue. When the check finds no
error, `sn3d` writes these lines to the standard output:

- the counts of the model and of its atomic data:
  - the model cells with matter and a propagation cell;
  - the elements of `compositiondata.txt`;
  - the ions, the levels, and the lines of their atomic data;
  - the excited NLTE levels and the superlevels;
- the largest peak memory of one rank in the check. For a new simulation, the line also gives the memory of one
  packet. Each rank of a run needs that memory for each of its packets;
- `input check: complete`.

The peak memory includes the arrays that the ranks of a node share, so it changes with the rank count of the check.

The check works in the folder `inputcheck` of the model folder, and it changes no other file of the model folder. The
folder holds the log of each rank, e.g. `inputcheck/output_0-0.txt`, with the `mem_usage` lines for the memory of the
arrays. It also holds the files that a simulation writes before its first timestep, e.g. `timesteps.out`. Each check
removes the folder of the check before it, so run only one check at a time in a model folder. A new simulation and
`scripts/clean.sh` remove the folder too.

A check differs from a run in these points:

- When a run would restore `input.txt` from `input-newrun.txt`, the check reads `input-newrun.txt` directly.
  See the [order of the parameter files](inputfiles.md#inputtxt).
- The check of a new simulation needs little memory for the packets, so a low rank count is sufficient, e.g. on a
  login node.
- The check of a resumed simulation needs the rank count and the memory of the simulation, because each rank reads
  its own {term}`restart files <restart file>`.

## A simulation in several jobs

To split a simulation across jobs, give `sn3d -w` the available wall time in hours.
The [cluster scripts](clusters.md#what-the-job-scripts-do) get this value from Slurm.
From the second timestep of a job, `sn3d` saves the grid and the packets before it tests the remaining wall time.
Read the log of rank 0 to decide what to do next.

| Text in the log | Next step |
| --- | --- |
| `RESTART_NEEDED` | The job has too little time for the next timestep. The cluster script submits the next job. For a manual run, use the current `input.txt` and the same rank count. |
| `Completed the requested timestep range` | The job reached `timestep_finish` before the end of the model. To continue, increase only `timestep_finish` in `input.txt`. Then submit the next job. |
| `No need for restart` | The job reached the end of the model. The spectra and light curves are ready for [analysis](analysis.md). |

For example, a new job with the range `000 009` propagates timesteps 0 to 8.
Its last restart files hold the state before packet propagation at timestep 8.
To continue a 20-timestep model, change `008 009` to `008 020` in the current `input.txt`.
Keep the continuation flag at 1. The next job propagates timestep 8 again.

Run `artis/sn3d -h` to list all command-line options. `artis/sn3d --version` gives the git version, the preset, and
the compiler of the build. An unknown option or argument stops `sn3d` before it reads a file.

## Check a run

Examine these points during and after a run:

- the progress line after each timestep, with its estimate of the remaining time (see
  [Example: a nebular spectrum](example.md#7-follow-the-run));
- the `mem_usage` lines, against the memory of a node (see
  [Plan a simulation](planning.md#measure-the-cost-of-your-run));
- the warnings of the logs. Most runs give some of them, and some need a change of the input (see
  [Troubleshooting](troubleshooting.md#warnings-that-most-runs-give));
- the electron temperature T_e of the cells after the first timesteps. It must be between `MINTEMP` and `MAXTEMP`, and
  stable at the times of the results that you use (see [Plan a simulation](planning.md#the-first-timesteps) and
  [Example: a nebular spectrum](example.md#8-check-the-results));
- the rows of `artistools check` (see [Analyse the results](analysis.md#start-with-the-checks));
- the packet noise of the spectra (see [Analyse the results](analysis.md#is-the-packet-count-sufficient)).

The log of a complete job ends with the line `sn3d finished (job: ...)`. A job without that line had a crash or a kill
(see [Troubleshooting](troubleshooting.md#a-job-stopped-with-a-crash-or-a-memory-kill)).

## The ranks and the threads

Each MPI rank holds an equal share of `NUM_PACKETS` in memory, and it updates its own range of the cells (see
[Plan a simulation](planning.md#the-ranks)).

A build with `OPENMP=ON` or with `STDPAR=ON` can add threads inside each rank, and `STDPAR=ON GPU=ON` runs the loops
on a GPU. [Make options](makeoptions.md) says which compiler gives threads with each option. A multithreaded run gives no reproducible results.
