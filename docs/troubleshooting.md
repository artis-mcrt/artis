# Troubleshooting

This page gives the messages that stop a build or a run, and what each one needs. The log file of rank 0 holds the git
version and the compile-time options of a run, so give both in an [issue](https://github.com/artis-mcrt/artis/issues).

## The build

### make cannot find artisoptions.h

The build needs an `artisoptions.h` in the repository root, and git does not track that file. Symlink a preset to it:

```sh
ln -s artisoptions_classic.h artisoptions.h
```

[Option presets](presets.md) says what each preset is for. Use a symlink and not a copy, because the symlink follows
each change of the preset file.

### A change of the preset compiles nothing again

Make compares the timestamps through a symlink, so a new symlink starts no rebuild. Give make the flag that reads the
timestamp of the symlink itself:

```sh
export MAKEFLAGS="--check-symlink-times --jobs=$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
```

The build warns you when `artisoptions.h` is a symlink and the flag is absent. A copy gets no such warning, which is
the second reason for the symlink. After a copy with the old timestamp, e.g. from `cp -p`, run `make clean`.

### make cannot run 'mpicxx --version'

The build needs an MPI library that supplies `mpicxx`. Install one, e.g. Open MPI, or load the module of your cluster.
The job scripts in `scripts/` name the modules of some clusters. With Open MPI, select the compiler with an
environment variable, e.g. `export OMPI_CXX=g++-16`.

### gcc is older than the minimum version

The build needs gcc 14 or a later version. Load a newer gcc. With Open MPI you can also point `OMPI_CXX` at it.

### The binary does not hold your change

`make` puts the objects in `build/<configuration>/` and makes a symlink in the repository root to each program. Each
option adds a suffix to the name of the build folder, so the builds of different options stay apart. The symlink
points at the most recent build of any configuration. Look at the symlink before you copy a binary, and copy it with
`cp -L` from the repository root. A build folder of an older compiler otherwise gives you the old binary.

## The run

### The photoionisation data needs phixsdata.txt or phixsdata_v2.txt

The model folder holds no photoionisation file. [Atomic and model data](data.md) gives the datasets and the download
command. A run also needs `adata.txt`, `transitiondata.txt`, and `compositiondata.txt`.

### Not enough slots available

Open MPI stops when the rank count is higher than the core count of the machine. Give `mpirun` the option
`--oversubscribe`, as `tests/runtest.sh` does.

### A run with TESTMODE=ON crashes in the MPI library

A build with `TESTMODE=ON` has the sanitizers, so a run needs two settings:

```sh
export ASAN_OPTIONS=detect_stack_use_after_return=1:detect_leaks=0
export OMPI_MCA_memory='^patcher'
```

`tests/runtest.sh` sets both when `TESTMODE=ON` is in the environment.

### The log ends with RESTART_NEEDED

This is the clean end of a job, and not an error. `sn3d` writes that text when too little wall time remains for
another timestep. It then writes the restart files and updates `input.txt`, so that the next job continues the
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

artistools accepts only the times at which the whole model can have sent light to the observer. The fastest ejecta set
that limit, so the limit is earlier than the last timestep of the run. `--plotinvalidpart` draws the later times.

### The checksums of a local test differ from the reference files

Continuous integration makes the reference checksums on an arm64 runner with one gcc version. Another CPU type or
another compiler can give other values, and no test covers that case. Push the branch and let CI decide. See
[Testing](testing.md).
