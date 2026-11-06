# Running sn3d

Run `sn3d` in the model folder. A local run with eight MPI ranks looks like this:

```sh
mpirun -np 8 artis/sn3d &
tail -F output_0-0.txt
```

`tail -F` waits until `sn3d` makes the log file. Press Ctrl+C to stop following it.

## The logs

`sn3d` writes one line with the name of the job folder to the standard output, and then nothing more unless a crash
occurs. Each MPI rank n writes a log file `output_n-t.txt` into the job folder, where t is the thread of that rank (see
[Output files](outputfiles.md)). A run without threads writes only the files of thread 0. `sn3d` keeps a symlink
`output_0-0.txt` in the model folder that points to the log of rank 0 and thread 0.

The log of rank 0 starts with the git version and the content of `artisoptions.h` without the comment lines, the
preprocessor lines, and the blank lines. After each timestep, it gives the progress of the model and an estimate of the
time for the remaining requested timesteps. The estimate uses the duration of the last timestep.

## Check the input before a job

Run `sn3d` with the option `-c` before you submit a job. The check reads the input files as a run does, and it stops
before the first timestep:

```sh
mpirun -np 4 artis/sn3d -c
```

An error in an input file then stops the check, and not a job that waited in the queue. `sn3d` writes
`input check: complete` to the standard output when it finds no error.

The check works in the folder `inputcheck` of the model folder, and it changes no other file of the model folder. The
folder holds the log of each rank, e.g. `inputcheck/output_0-0.txt`, with the `mem_usage` lines for the memory of the
arrays. It also holds the files that a simulation writes before its first timestep, e.g. `timesteps.out`. Each check
removes the folder of the check before it, so run only one check at a time in a model folder. A new simulation and
`scripts/clean.sh` remove the folder too.

A check differs from a run in these points:

- The check needs `input.txt`. A run restores an absent `input.txt` from `input-newrun.txt`, but the check does not.
- The check of a new simulation places only one pellet for each rank. A low rank count is thus sufficient, e.g. on a
  login node. The `mem_usage` line of the packets gives the share of one rank for the rank count of the check. The
  arrays that the ranks of a node share have the size that they have in the run.
- The check of a resumed simulation needs the rank count of the simulation, because each rank has its own restart
  files. Each rank reads its restart files of the start timestep, so this check needs the memory of the run.
- With `VPKT_WRITE_CONTRIBS`, the check of a resumed simulation copies the `vpackets` restart file of each rank into
  the folder `inputcheck`.

## A simulation in several jobs

To split a long simulation across several queued jobs, run `sn3d` with `-w WALLTIMELIMITHOURS`. When too little wall
time remains to complete another timestep, the run finishes cleanly and prints `RESTART_NEEDED` into the log. It then
writes the restart files and updates `input.txt`. The bundled cluster job scripts detect that text and submit a
continuation job, and they pass the remaining Slurm allocation time automatically (see
[Runs on a cluster](clusters.md)).

Run `artis/sn3d -h` to list all command-line options.

## The ranks and the threads

Each MPI rank gets an equal share of `NUM_PACKETS`, and the first ranks also get one packet each of the remainder. Each
rank keeps its share in memory, so the rank count and `NUM_PACKETS` together set the memory of a run.

The ranks on one node share the large read-mostly arrays, e.g. the atomic data, through MPI shared-memory windows. A
run with `OPENMP=ON` adds OpenMP threads inside each rank, and a run with `STDPAR=ON` gives the loops to the
parallel algorithms of the standard library. Those algorithms use threads only with a compiler and a standard library
that supply a parallel backend (see [Make options](makeoptions.md)). A multithreaded run gives no reproducible
results, because each thread seeds its own random generator. A run with
`STDPAR=ON GPU=ON` is different. The parallel algorithms then run on a GPU and not on the threads of the rank, and
each packet holds its own generator state.
