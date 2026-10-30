# Running sn3d

Run `sn3d` in the model folder. A local run with eight MPI ranks looks like this:

```sh
mpirun -np 8 artis/sn3d &
tail -f output_0-0.txt
```

Press Ctrl+C to stop following the log file.

## The logs

`sn3d` writes one line with the name of the job folder to the standard output, and then nothing more unless a crash
occurs. Each MPI rank n writes a log file called `output_n-0.txt` into the job folder (see
[Output files](outputfiles.md)). `sn3d` keeps a symlink `output_0-0.txt` in the model folder that points to the log of
rank 0.

The log of rank 0 starts with the git version and the content of `artisoptions.h` without the comment lines, the
preprocessor lines, and the blank lines. After each timestep, it gives the progress of the model and an estimate of the
time for the remaining requested timesteps. The estimate uses the duration of the last timestep.

## A simulation in several jobs

To split a long simulation across several queued jobs, run `sn3d` with `-w WALLTIMELIMITHOURS`. When too little wall
time remains to complete another timestep, the run finishes cleanly and prints `RESTART_NEEDED` into the log. It then
writes the restart files and updates `input.txt`. The bundled cluster job scripts detect that text and submit a
continuation job, and they pass the remaining Slurm allocation time automatically (see
[Runs on a cluster](clusters.md)).

Run `./sn3d -h` to list all command-line options.

## The ranks and the threads

Each MPI rank gets an equal share of `NUM_PACKETS`, and the first ranks also get one packet each of the remainder. Each
rank keeps its share in memory, so the rank count and `NUM_PACKETS` together set the memory of a run.

The ranks on one node share the large read-mostly arrays, e.g. the atomic data, through MPI shared-memory windows. A
run with `OPENMP=ON` or with `STDPAR=ON` adds threads inside each rank (see [Make options](makeoptions.md)). A
multithreaded run gives no reproducible results, because each thread seeds its own random generator.
