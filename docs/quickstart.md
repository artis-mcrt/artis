# Quickstart

This page gives a first simulation, a check of the results, and a plot. Everything here runs on a laptop. The test model takes a few minutes on four cores.

## 1. Install a compiler and MPI

On macOS:

```sh
brew install llvm open-mpi pkgconf zstd coreutils
export PATH="$(brew --prefix llvm)/bin:$(brew --prefix coreutils)/libexec/gnubin:$PATH"
```

On Debian, and on Ubuntu 24.04 or a later version:

```sh
sudo apt-get install g++-14 openmpi-bin libopenmpi-dev libzstd-dev zstd git make curl xz-utils rsync python3 perl
export OMPI_CXX=g++-14
```

On Fedora 40 or a later version:

```sh
sudo dnf install gcc-c++ openmpi-devel libzstd-devel zstd git make curl xz rsync python3 perl
module load mpi/openmpi-$(uname -m)
```

On a different Linux system, get the same packages from your package manager or from a module.
[Installation](installation.md) gives the versions and the commands that the test of step 3 calls.

## 2. Clone the repository

```sh
git clone https://github.com/artis-mcrt/artis.git
cd artis
```

## 3. Run a test model

```sh
tests/runtest.sh kilonova_1d
```

This is the quickest complete simulation of the repository. The script downloads the atomic data of the test (about
15 MB), builds `sn3d` and `exspec`, and makes the model folder `tests/kilonova_1d_testrun`. It then runs the
simulation in two jobs, and it compares the checksums of the output files with the stored checksums.

A different compiler or CPU type can give different checksums, and the script then stops at the comparison. The
output files of the first job are complete at that point, so the next section still works. [Testing](testing.md)
describes the steps and how you run one step by name.

`tests/runtest.sh` without arguments lists the other test models.

## 4. Look at the results

The model folder `tests/kilonova_1d_testrun` now holds the output files, e.g. `light_curve.out.zst` and `spec.out.zst`. [Output files](outputfiles.md) describes each one.

[artistools](https://github.com/artis-mcrt/artistools) plots and analyses them. It needs
[uv](https://docs.astral.sh/uv/), which section 1 does not install:

```sh
brew install uv                                  # macOS
curl -LsSf https://astral.sh/uv/install.sh | sh  # Linux
```

The Linux installer adds `uv` to the profile of your shell, but the current shell needs it too:

```sh
source "$HOME/.local/bin/env"  # Linux
```

Then install artistools:

```sh
uv tool install artistools
export PATH="$(uv tool dir --bin):$PATH"
```

The export gives the tool folder of `uv` to the current shell, and `uv tool update-shell` adds it to the profile of
your shell.

Then describe the model, check the results, and plot them:

```sh
artistools describe tests/kilonova_1d_testrun
artistools check tests/kilonova_1d_testrun
artistools plotlightcurves tests/kilonova_1d_testrun
artistools plotspectra -t 1 tests/kilonova_1d_testrun
```

:::{figure} images/lightcurve.png
:alt: The light curve of the test model

The light curve that `plotlightcurves` gives. The test model uses few packets, so the curve is noisy.
:::

:::{figure} images/spectrum.png
:alt: The spectrum of the test model at one day

The spectrum that `plotspectra -t 1` gives.
:::

`artistools check` gives one row for each check of the results, so it is the quickest test that a run deserves trust.
It gives a warning for this test model, because the reduced atomic data and the low packet count limit its accuracy.

`plotspectra` needs a time in days (`-t`). artistools plots only the times that have light from the whole model,
which end at about 1.7 days for this model (see [Troubleshooting](troubleshooting.md#the-results)).

## 5. Make your own model folder

Each simulation has its own model folder. Keep the source code in a subfolder `artis` of it. The code then finds its
physics data in `artis/data`. The commands below make a model folder from the input files of the test model. Run the
first command in the clone of step 2:

```sh
scripts/new-model-folder.sh --symlink ~/mymodel
cd ~/mymodel
tar -xf artis/tests/atomicdata_feconi.tar.xz
rsync -a --exclude 'results_md5_*' artis/tests/kilonova_1d_inputfiles/ ./
cp input-newrun.txt input.txt
```

The script makes the folder, a symlink `artis` to the clone, and an `input.txt` from `scripts/input-template.txt`.
The last command replaces that `input.txt` with the file of the test model.

For a production run, clone the branch `release` into the model folder in place of the symlink (see [Runs on a cluster](clusters.md)).

The atomic data of the test covers only a few elements, with few levels. [Atomic and model data](data.md) says where a larger dataset and other models come from.

Now select the compile-time options and build:

```sh
cd artis
ln -s artisoptions_kilonova_lte.h artisoptions.h
make
cd ..
```

[Option presets](presets.md) says what each {term}`preset` is for.

:::{note}
`artisoptions.h` is a symlink, so an edit changes the preset file of the clone, and each model folder that shares the
clone gets the same options. Give a model folder its own clone when it needs its own options.
:::

The preset has a value of `NUM_PACKETS` for about 1000 MPI {term}`ranks <rank>`. A laptop needs a much smaller value. Change it in
`artis/artisoptions.h`:

```c++
constexpr std::int64_t NUM_PACKETS = 1e6;
```

Also set the range of the {term}`timesteps <timestep>` in `input.txt`. The line with the comment `#  2: timestep_start
timestep_finish` holds the first timestep of the job and the stop timestep. The job does not run the stop timestep, so the test runs
the first nine of twenty timesteps. Change that line to run the whole model:

```text
000 020                  #  2: timestep_start timestep_finish: ...
```

Then build again and run:

```sh
cd artis && make && cd ..
mpirun -np 4 artis/sn3d &
tail -F output_0-0.txt
```

`output_0-0.txt` is the log of rank 0, and `tail -F` waits until `sn3d` makes it. Press Ctrl+C to stop following the
log file. [Running sn3d](running.md) describes the command-line options, the logs, and the restart of a long run.

## Next steps

- [Option presets](presets.md) and [Compile-time options](options.md): the physics of the run.
- [Input files](inputfiles.md): the format of `model.txt`, `abundances.txt`, and the atomic data files.
- [Runs on a cluster](clusters.md): the Slurm job scripts and the wall-time mechanism.
