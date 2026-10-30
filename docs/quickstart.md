# Quickstart

This page gives a first simulation, a check of the results, and a plot. Everything here runs on a laptop. The test model takes a few minutes on four cores.

## 1. Install a compiler and MPI

On macOS:

```sh
brew install llvm open-mpi pkgconf zstd coreutils
```

`tests/runtest.sh` needs the `md5sum` command of coreutils for the next step.

On Linux, get a compiler and an MPI library from your package manager, from a module, or from pixi. [Installation](installation.md) gives the detail and the version that the build needs.

## 2. Clone the repository

```sh
git clone https://github.com/artis-mcrt/artis.git
cd artis
```

## 3. Run a test model

```sh
tests/runtest.sh kilonova_1d
```

This is the quickest complete simulation of the repository. The script does these steps:

1. It downloads the atomic data of the test (about 15 MB) into `tests/` and keeps the archive for a later run.
2. It makes the model folder `tests/kilonova_1d_testrun` with the input files and an `artisoptions.h` for the test.
3. It builds `sn3d` and `exspec`, and copies them into the model folder.
4. It runs a new simulation, resumes it from the restart files, and then runs `exspec`.
5. It compares the checksums of the output files with the stored reference checksums.

A different compiler or CPU type can give checksums that differ from the stored checksums. The simulation is still valid. [Testing](testing.md) describes the steps and the reference checksums.

The repository has eleven test models. `tests/runtest.sh` without arguments lists them.

## 4. Look at the results

The model folder `tests/kilonova_1d_testrun` now holds the output files, e.g. `light_curve.out.zst` and `spec.out.zst`. [Output files](outputfiles.md) describes each one.

[artistools](https://github.com/artis-mcrt/artistools) plots and analyses them. Install it with [uv](https://docs.astral.sh/uv/):

```sh
uv tool install artistools
```

Then describe the model, check the results, and plot them:

```sh
artistools describe tests/kilonova_1d_testrun
artistools check tests/kilonova_1d_testrun
artistools plotlightcurves tests/kilonova_1d_testrun
artistools plotspectra -t 1 tests/kilonova_1d_testrun
```

`artistools check` gives one row for each check of the results, so it is the quickest test that a run deserves trust.
It gives a warning for this test model, because the reduced atomic data and the low packet count limit its accuracy.

`plotspectra` needs a time in days (`-t`). The test runs only the first ten of twenty timesteps, so the spectra are
valid up to about 1.7 days.

## 5. Make your own model folder

Each simulation has its own model folder. Keep the source code in a subfolder `artis` of it. The code then finds its physics data in `artis/data`. The commands below make a model folder from the input files of the test model:

```sh
mkdir ~/mymodel
cd ~/mymodel
ln -s ~/artis artis
tar -xf artis/tests/atomicdata_feconi.tar.xz
rsync -a --exclude 'results_md5_*' artis/tests/kilonova_1d_inputfiles/ ./
cp input-newrun.txt input.txt
```

The archive and the input folder both hold a `compositiondata.txt`. The archive comes first, so the file of the test
model wins, as it does in the test. That file keeps ten elements, and the file of the archive keeps three.

For a production run, clone the repository into the model folder in place of the symlink (see [Runs on a cluster](clusters.md)).

The atomic data of the test covers only a few elements, with few levels. [Atomic and model data](data.md) says where a larger dataset and other models come from.

Now select the compile-time options and build:

```sh
cd artis
ln -s artisoptions_kilonova_lte.h artisoptions.h
make
cd ..
```

[Option presets](presets.md) says what each preset is for.

:::{note}
`artisoptions.h` is a symlink, so an editor changes the preset file of your clone. That is the intended workflow, because each model folder keeps its own clone of the source. Do not replace the symlink with a copy.
:::

The preset has a value of `NUM_PACKETS` for about 1000 MPI ranks. A laptop needs a much smaller value. Change it in `artis/artisoptions.h`:

```c++
constexpr std::int64_t NUM_PACKETS = 1e6;
```

Also set the range of the timesteps in `input.txt`. Line 2 holds the first and the last timestep of the job. The test runs only the first ten of twenty timesteps, so change that line to run the whole model:

```text
000 020                  #  2: timestep_start timestep_finish: ...
```

Then build again and run:

```sh
cd artis && make && cd ..
mpirun -np 4 artis/sn3d &
tail -f output_0-0.txt
```

`sn3d` writes one line with the name of the job folder to the standard output, and then nothing more unless a crash occurs. Each rank writes its own log file. `output_0-0.txt` in the model folder is a symlink to the log of rank 0. Press Ctrl+C to stop following the log file.

[Running sn3d](running.md) describes the command-line options, the logs, and the restart of a long run.

## Next steps

- [Option presets](presets.md) and [Compile-time options](options.md): the physics of the run.
- [Input files](inputfiles.md): the format of `model.txt`, `abundances.txt`, and the atomic data files.
- [Runs on a cluster](clusters.md): the Slurm job scripts and the wall-time mechanism.
