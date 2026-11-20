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

- A run restores an absent `input.txt` from `input-newrun.txt`. The check restores no file. Without `input.txt`, it
  reads `input-newrun.txt` directly.
- The check of a new simulation needs little memory for the packets, so a low rank count is sufficient, e.g. on a
  login node.
- The check of a resumed simulation needs the rank count and the memory of the simulation, because each rank reads
  its own {term}`restart files <restart file>`.

## A simulation in several jobs

To split a long simulation across several queued jobs, run `sn3d` with `-w WALLTIMELIMITHOURS`. From the second
timestep of a job, `sn3d` writes the restart files and updates `input.txt` after the grid update of each timestep. Then it
tests the remaining wall time. When too little wall time remains to complete another timestep, the job stops cleanly, and the
log then ends with `RESTART_NEEDED`. The bundled cluster job scripts detect that text and submit a
continuation job, and they pass the remaining Slurm allocation time automatically (see
[Runs on a cluster](clusters.md)).

Run `artis/sn3d -h` to list all command-line options. `artis/sn3d --version` gives the git version, the preset, and
the compiler of the build. An unknown option or argument stops `sn3d` before it reads a file.

## The ranks and the threads

Each MPI rank holds an equal share of `NUM_PACKETS` in memory, and it updates its own range of the cells (see
[Plan a simulation](planning.md#the-ranks)).

A build with `OPENMP=ON` or with `STDPAR=ON` can add threads inside each rank, and `STDPAR=ON GPU=ON` runs the loops
on a GPU. [Make options](makeoptions.md) says which compiler gives threads with each option. A multithreaded run gives no reproducible results.
