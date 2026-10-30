# Installation

## What the build needs

- A C++ compiler: gcc 14 or a later version, clang, nvc++, or hipcc.
- An MPI library that supplies an `mpicxx` command, e.g. Open MPI.
- Optionally libzstd. A build with the library reads and writes zstd compressed files (see [Input files](inputfiles.md)).

The build needs no other program. On an HPC cluster, a module or a spack command gives the compiler and the MPI library. For the clusters that we use, the top of the Slurm script in `scripts/artis-*.sh` names the modules that work.

With Open MPI, select the compiler with an environment variable:

```sh
export OMPI_CXX=g++
```

## macOS

Install the LLVM tools, Open MPI, and the other dependencies with homebrew:

```sh
brew install llvm open-mpi pkgconf zstd coreutils prek
```

The homebrew LLVM gives clang, clang-format, clang-tidy, and the clangd language server. coreutils gives the `md5sum` command, which the end-to-end tests need (see [Testing](testing.md)).

## Linux

`pixi run make` builds with a gcc from conda-forge (see `pixi.toml`):

- The MPI library must come from the system, e.g. from a module on a cluster.
- To use the gcc of pixi for the commands below, start `pixi shell` first.
- If the link stops with undefined references from the MPI library, export `LDFLAGS` with `-Wl,-rpath-link` and the library folder of the MPI library. The job scripts for Leonardo and Virgo show examples.

## Parallel builds

Add this line to your shell startup file:

```sh
export MAKEFLAGS="--check-symlink-times --jobs=$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
```

macOS has no `nproc` command, so the `sysctl` command gives the processor count there. Without `--check-symlink-times`, make does not rebuild after a change of the symlink to a different preset.

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

`make` builds `sn3d` and `exspec`. It puts the objects in `build/<configuration>/` and makes `sn3d`, `exspec`, and `unittests` in the repository root symlinks to that folder. Each make option adds a suffix to the name of the build folder, so different sets of make options coexist. The name holds no preset, so a change of the preset must compile every object again. `--check-symlink-times` in `MAKEFLAGS` makes make see that change. Each symlink points to the most recent build of any configuration, so look at the symlink before you copy a binary. Copy a binary with `cp -L` from the repository root.

[Make options](makeoptions.md) lists the options of `make`.
