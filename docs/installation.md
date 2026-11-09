# Installation

## What the build needs

- git and make.
- A C++ compiler: gcc 14 or a later version, clang, nvc++, or hipcc.
- An MPI library that supplies an `mpicxx` command, e.g. Open MPI.
- libzstd (optional). A build with the library reads and writes zstd compressed files (see
  [Input files](inputfiles.md)).
- The TBB library, but only for a build with `STDPAR=ON` and gcc. The package is `libtbb-dev` on Debian and Ubuntu.
- The OpenMP runtime library, but only for a build with `OPENMP=ON` and clang. The package is `libomp` with homebrew
  and `libomp-dev` on Debian and Ubuntu.

The end-to-end tests of `tests/runtest.sh` (see [Testing](testing.md)) need more:

- Open MPI itself, because the script gives `mpirun` options that only Open MPI accepts;
- libzstd, which is the package `libzstd-dev` on Debian and Ubuntu, and `libzstd-devel` on Fedora;
- the commands `zstd`, `md5sum`, `curl`, `tar`, `xz`, `python3`, `perl`, and `rsync`.

[uv](https://docs.astral.sh/uv/) is separate from the build. It gives artistools and runs the Python helper scripts of `scripts/`. On an HPC cluster, a module or a spack command gives the compiler and the MPI library. For the clusters that we use, the top of the Slurm script in `scripts/artis-*.sh` names the modules that work.

With Open MPI, select the compiler with an environment variable:

```sh
export OMPI_CXX=g++
```

## macOS

Install the LLVM tools, Open MPI, and the other dependencies with homebrew:

```sh
brew install llvm lld libomp open-mpi pkgconf zstd coreutils prek
```

Homebrew keeps the LLVM tools and the GNU tools of coreutils out of the standard path, so put their two folders
first:

```sh
export PATH="$(brew --prefix llvm)/bin:$(brew --prefix coreutils)/libexec/gnubin:$PATH"
```

The homebrew LLVM gives clang, clang-format, clang-tidy, and the clangd language server. coreutils gives `nproc` and
`md5sum`.

## Linux

`pixi run make` builds with a gcc from conda-forge (see `pixi.toml`):

- The MPI library must come from the system, e.g. from a module on a cluster.
- To use the gcc of pixi for the commands below, start `pixi shell` first.
- If the link stops with undefined references from the MPI library, export `LDFLAGS` with `-Wl,-rpath-link` and the library folder of the MPI library. The job scripts for Leonardo and Virgo show examples.

## Parallel builds

Add this line to your shell startup file:

```sh
export MAKEFLAGS="--jobs=$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
```

## Select a preset and build

The build needs an `artisoptions.h` in the repository root. That file is gitignored. Symlink one of the presets to it:

```sh
ln -s artisoptions_classic.h artisoptions.h
make
```

[Option presets](presets.md) says what each preset is for and why it is a symlink, and [Compile-time options](options.md) describes every option.

`make` builds `sn3d` and `exspec`, and it makes a symlink in the repository root to each program. Each set of make options and each content of `artisoptions.h` has its own build folder (see [Make options](makeoptions.md)).

`make help` prints one line for each target and each option of `make`.
