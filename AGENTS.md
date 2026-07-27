# AGENTS.md

Guidance for AI coding agents working on ARTIS. See [README.md](README.md) for a
full project description.

## Project overview

ARTIS is a 3D Monte Carlo radiative transfer code for supernovae and kilonovae,
written in modern C++ (C++26 with gcc/clang; code must remain compatible with
C++23 for the nvc++ and hipcc GPU builds) and parallelised with MPI (plus
optional OpenMP or C++ standard parallelism). Two executables are built from the
same sources:

- `sn3d`: the main simulation (radiative transfer over a time-step sequence)
- `exspec`: post-processes packet files into spectra and light curves

## Repository layout

- `*.cc` / `*.h` in the repository root: the entire C++ source. Each module pair
  (e.g. `grid.cc`/`grid.h`) exposes its interface in a matching namespace
  (`grid::`, `decay::`, `radfield::`, ...). `sn3d.cc` and `exspec.cc` contain the
  two `main()` functions; everything else is shared.
- `artisoptions_*.h`: compile-time option presets. The build requires an
  `artisoptions.h` in the repo root (gitignored), normally a copy or symlink of
  one of the presets. Options are documented in `artisoptions_doc.md`.
- `data/`: physics data bundled with the code (decay data, gamma spectra, etc.).
- `tests/`: integration test setups (`setup_*.sh` plus `*_inputfiles/` folders).
- `scripts/`: HPC cluster job scripts and helper utilities.
- `third_party/`: vendored Eigen and Boost. Never edit or lint these.

## Building

Requires a recent C++ compiler (gcc >= 14, clang, or nvc++) and an MPI library
providing `mpicxx`. With Open MPI, select the underlying compiler via
`export OMPI_CXX=g++`.

```sh
cp artisoptions_classic.h artisoptions.h   # or another preset; required, gitignored
make -j$(nproc)                            # builds sn3d and exspec
```

Objects go into a per-configuration `build/` subdirectory and the binaries are
symlinked into the repo root. `make OPTIMIZE=OFF` is the fastest way to check
that everything compiles (this is what the pre-commit hook runs).

Useful make options (see README for the full list): `TESTMODE=ON` (assertions +
address/UB sanitizers), `REPRODUCIBLE=ON`, `FASTMATH=OFF`, `OPTIMIZE=OFF`,
`EIGEN=OFF` (use GSL instead of Eigen), `OPENMP=ON`, `STDPAR=ON`, `GPU=ON`,
`MAX_NODE_SIZE=N`.

The build uses `-Werror` with extensive warnings — new warnings break CI.

## Testing

There are no unit tests. Testing is done through end-to-end regression runs
(`.github/workflows/ci.yml`): each test in `tests/` runs a small model with
`REPRODUCIBLE=ON FASTMATH=OFF MAX_NODE_SIZE=2` and compares `md5sum` checksums of
the output files against `tests/<testname>_inputfiles/results_md5_*.txt`. To run
one locally (downloads atomic data on first use):

```sh
cd tests
source ./setup_kilonova_1d.sh        # creates tests/kilonova_1d_testrun/
cd ..
make REPRODUCIBLE=ON MAX_NODE_SIZE=2 FASTMATH=OFF -j$(nproc) sn3d exspec
cp sn3d exspec tests/kilonova_1d_testrun/
cd tests/kilonova_1d_testrun
cp input-newrun.txt input.txt
mpirun -np 4 --oversubscribe ./sn3d   # logs go to output_0-0.txt, not stdout
```

If a change legitimately alters numerical results, the stored checksums must be
regenerated — maintainers do this via the "Update checksums" workflow
(`.github/workflows/updatechecksums.yml`), which commits the checksum files
produced by the last CI run of the branch. Say so in the PR/commit message when
results are expected to change.

CI also compiles every `artisoptions_*.h` preset with gcc, clang, Apple Clang,
nvc++, and hipcc — so a new compile-time option must be added to **all** preset
files and documented in `artisoptions_doc.md`, and code must compile under all
of those compilers (including `GPU=ON` paths).

## Linting and formatting

- clang-format (Google-based style, 120-column limit; config in `.clang-format`)
  is enforced on all C++ files.
- clang-tidy with an extensive check list (`.clang-tidy`); run via `make check`
  (needs `compile_commands.json`, e.g. `compiledb -n make TESTMODE=ON`).
- cppcheck runs in CI (`.github/workflows/ci-checks.yml`).
- Pre-commit hooks (`.pre-commit-config.yaml`, managed with `prek`/pre-commit)
  run clang-format, whitespace/EOF fixers, and a `make OPTIMIZE=OFF` compile.

Before committing, at minimum run clang-format on changed files and make sure
`make OPTIMIZE=OFF` succeeds.

## Conventions and gotchas

- Prefer modern C++: `constexpr`, `std::ranges`, `std::span`, structured
  bindings. Compile-time configuration lives in `artisoptions.h` as `constexpr`
  values rather than runtime flags.
- Global simulation state lives in `namespace globals` (`globals.h`). Large
  read-mostly arrays use MPI shared-memory windows (`MPI_shared_array`) so that
  ranks on a node share one copy — be careful about who writes to them and when.
- `sn3d` writes per-rank log files (`output_<rank>-0.txt`); stdout stays quiet
  unless there is a crash.
- The Makefile regenerates `version.h` from git metadata on every invocation;
  `version*.h` is gitignored.
- The default branch is `develop`; `release` is the production branch. CI runs
  on every push.
