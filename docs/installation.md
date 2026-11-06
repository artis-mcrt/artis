# Installation

## What the build needs

- git and make.
- A C++ compiler: gcc 14 or a later version, clang, nvc++, or hipcc.
- An MPI library that supplies an `mpicxx` command, e.g. Open MPI. The end-to-end tests need Open MPI itself, because `tests/runtest.sh` gives `mpirun` the options `--oversubscribe` and `--mca` (see [Testing](testing.md)).
- libzstd. A build with the library reads and writes zstd compressed files (see [Input files](inputfiles.md)). A build without it stops at the end-to-end tests, because `tests/runtest.sh` builds with `ZSTD=ON` and the test input files are compressed.
- The TBB library, but only for a build with `STDPAR=ON` and gcc. The Makefile then links `-ltbb`. The package is `libtbb-dev` on Debian and Ubuntu. The other compilers need no such library.
- The OpenMP runtime library, but only for a build with `OPENMP=ON` and clang. The Makefile then links `-lomp`. The
  package is `libomp` with homebrew and `libomp-dev` on Debian and Ubuntu.

The build needs no other program. [uv](https://docs.astral.sh/uv/) is separate from the build. It gives artistools and runs the Python helper scripts of `scripts/`. On an HPC cluster, a module or a spack command gives the compiler and the MPI library. For the clusters that we use, the top of the Slurm script in `scripts/artis-*.sh` names the modules that work.

With Open MPI, select the compiler with an environment variable:

```sh
export OMPI_CXX=g++
```

## pixi on a laptop

[pixi](https://pixi.sh) installs a compiler, Open MPI, libzstd, and each command that the end-to-end tests call. It
needs no other package manager, and it works on Linux (x86-64) and on macOS (Apple silicon). Install pixi, and then
start a new terminal:

```sh
curl -fsSL https://pixi.sh/install.sh | sh
```

In the clone of the repository, start a shell of the `laptop` environment:

```sh
pixi shell -e laptop
```

The first start downloads the packages into the `.pixi` folder of the clone. The shell then holds `mpicxx`, `mpirun`,
`make`, and the other commands, so the commands of this site work in it without a change.
`pixi run -e laptop <command>` runs one command in the environment, e.g.
`pixi run -e laptop tests/runtest.sh kilonova_1d`.

- The environment uses gcc on Linux and clang on macOS.
- macOS also needs the command line tools of Xcode, because they hold the system headers. `xcode-select --install`
  installs them.
- On Linux, keep the clone in a path without a space. The gcc of conda-forge does not work in a path that holds a
  space.
- Do not use the `laptop` environment on an HPC cluster. Its Open MPI does not have the network support of the
  cluster, so use the MPI library of the cluster there (see [Linux](#linux)).

## macOS

Without pixi, install the LLVM tools, Open MPI, and the other dependencies with homebrew:

```sh
brew install llvm lld libomp open-mpi pkgconf zstd coreutils prek
```

Homebrew keeps the LLVM tools and the GNU tools of coreutils out of the standard path, so put their two folders
first:

```sh
export PATH="$(brew --prefix llvm)/bin:$(brew --prefix coreutils)/libexec/gnubin:$PATH"
```

The homebrew LLVM gives clang, clang-format, clang-tidy, and the clangd language server. coreutils gives the `nproc`
command of the parallel build, and it gives `md5sum`, which the end-to-end tests call (see [Testing](testing.md)).
Homebrew installs that command as `gmd5sum`, and the `gnubin` folder of the export holds it under the plain name.

## Linux

Get the compiler, the MPI library, and libzstd from your package manager or, on a cluster, from a module. The default
environment of pixi is an alternative source of the compiler only. `pixi run make` builds with its gcc from
conda-forge (see `pixi.toml`):

- The MPI library must come from the system, e.g. from a module on a cluster.
- To use the gcc of pixi for the commands below, start `pixi shell` first.
- If the link stops with undefined references from the MPI library, export `LDFLAGS` with `-Wl,-rpath-link` and the library folder of the MPI library. The job scripts for Leonardo and Virgo show examples.

## Parallel builds

Add this line to your shell startup file:

```sh
export MAKEFLAGS="--jobs=$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
```

macOS has no `nproc` command, so the `sysctl` command gives the processor count there.

## Select a preset and build

The build needs an `artisoptions.h` in the repository root. That file is gitignored. Symlink one of the presets to it:

```sh
ln -s artisoptions_classic.h artisoptions.h
make
```

[Option presets](presets.md) says what each preset is for, and [Compile-time options](options.md) describes every option.

:::{warning}
Symlink the preset. Do not copy it. The symlink keeps `artisoptions.h` correct when the preset changes, e.g. after a change of branch.
:::

`make` builds `sn3d` and `exspec`. It puts the objects in `build/<configuration>/` and makes a symlink in the repository root to each program that it builds. `make unittests` builds the third program and makes its symlink (see [Testing](testing.md)).

Each make option of [Make options](makeoptions.md) adds a suffix to the name of the build folder, so different sets of make options coexist. `GCCTOOLCHAIN` is the exception, so run `make clean` after you change its value. The name also ends with a checksum of the content of `artisoptions.h`, so the builds of different presets coexist too. A change of the preset or of an option value compiles every object in a new build folder. A change back to an earlier preset compiles nothing again.

Each edit of an option keeps the build folder of the old value. `make clean` removes all build folders.

Each symlink points to the most recent build of any configuration. Look at the symlink before you copy a binary. Copy a binary with `cp -L` from the repository root.

`make help` prints one line for each target and each option of `make`, and [Make options](makeoptions.md) describes them.
