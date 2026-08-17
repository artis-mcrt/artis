# AGENTS.md

Guidance for AI coding agents working on ARTIS. See [README.md](README.md) for a
full project description.

## Writing style

Write all English in ASD-STE100 (Simplified Technical English). This applies to
code comments, documentation, commit messages, pull request text, and replies to
the maintainers. It also applies to new log, warning, and error strings in the
code.

Obey these rules:

- Use the approved STE words, in the approved part of speech and the approved
  meaning. Technical names and technical verbs are permitted, e.g. `sn3d`,
  `MPI_shared_array`, "packet", "opacity", and "to sample".
- Use one term for one thing. Do not change between synonyms, e.g. between
  "cell" and "grid cell", or between "time step" and "timestep".
- Use the active voice. Write "the Makefile writes `version.h`" and not
  "`version.h` is written by the Makefile".
- Use the simple tenses (present, past, and future). Do not use the -ing form
  as a noun or as an adjective if a simple form is possible.
- Write short sentences. Use a maximum of 20 words in an instruction and a
  maximum of 25 words in descriptive text. Write a maximum of 6 sentences in a
  descriptive paragraph.
- Give one instruction in one sentence. Write the reason in a different
  sentence.
- Keep the articles "a", "an", and "the".
- Use a maximum of three words in a noun cluster. Write "the checksums of the
  output files" and not "output file checksum comparison".
- Write positive statements. Do not write a double negative.
- Do not use slang, idioms, or jokes. Do not use an abbreviation that the text
  does not define.
- Use a vertical list if the text contains more than three related items or
  conditions.

These rules do not apply to:

- Identifiers in the code, e.g. the names of variables, functions, and
  namespaces. Keep the conventions of the file that you change.
- The format and the content of the output files. `linestat.out` and the other
  `*.out` files must stay byte-compatible (see "Conventions and gotchas").
- Quoted text from an external source, e.g. a compiler message or a title of a
  publication.

Keep the British spellings that this repository uses, e.g. "parallelised" and
"normalise". STE controls the choice of words and the structure of the
sentences. It does not control the spelling variant here.

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
  `artisoptions.h` in the repo root (gitignored), normally a symlink to one of
  the presets. Options are documented in `artisoptions_doc.md`.
- `data/`: physics data bundled with the code (decay data, gamma spectra, etc.).
- `tests/`: integration test setups (`setup_*.sh` plus `*_inputfiles/` folders).
- `scripts/`: HPC cluster job scripts and helper utilities.
- `third_party/`: vendored Eigen (pruned to the modules ARTIS uses; see
  `third_party/eigenversion.txt`) and Random123. Never edit or lint these.

## Building

Requires a recent C++ compiler (gcc >= 14, clang, or nvc++) and an MPI library
providing `mpicxx`. With Open MPI, select the underlying compiler via
`export OMPI_CXX=g++`.

```sh
ln -s artisoptions_classic.h artisoptions.h   # or another preset; required, gitignored
make -j$(nproc)                               # builds sn3d and exspec
```

Symlinking (rather than copying) the preset keeps `artisoptions.h` in sync when
the preset changes, e.g. after switching branches. Add
`--check-symlink-times` to your make flags so that make notices preset edits
through the symlink.

Objects go into a per-configuration `build/` subdirectory and the binaries are
symlinked into the repo root. `make OPTIMIZE=OFF` is the fastest way to check
that everything compiles (this is what the pre-commit hook runs).

Useful make options (see README for the full list): `TESTMODE=ON` (assertions +
address/UB sanitizers), `REPRODUCIBLE=ON`, `FASTMATH=OFF`, `OPTIMIZE=OFF`,
`OPENMP=ON`, `STDPAR=ON`, `GPU=ON`, `MAX_NODE_SIZE=N`.

The build uses `-Werror` with extensive warnings — new warnings break CI.

## Testing

Unit tests for the pure numeric helpers live in `unittests.cc` (built with
`make unittests`, run as `./unittests`; CI runs them for the classic and nebular
presets), and compile-time checks of the `constexpr` helpers are `static_assert`s
next to their definitions. The main coverage is end-to-end regression runs
(`.github/workflows/ci.yml`): each test in `tests/` runs a small model with
`REPRODUCIBLE=ON FASTMATH=OFF MAX_NODE_SIZE=2` and compares `md5sum` checksums of
the output files against `tests/<testname>_inputfiles/results_md5_*.txt` (all
`*.out` files are checksummed; `output_*.txt` logs are not). Two separate smoke
test jobs run one model with `OPENMP=ON` and with `STDPAR=ON` (multicore CPU),
both without checksum comparison, to exercise the threaded code paths. To run a
regression test locally (downloads atomic data on first use):

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

CI also compiles every `artisoptions_*.h` preset with gcc and clang, and key
presets with Apple Clang, nvc++, and hipcc — so a new compile-time option must
be added to **all** preset files and documented in `artisoptions_doc.md`, and
code must compile under all of those compilers (including `GPU=ON` paths).

## Linting and formatting

- clang-format (Google-based style, 120-column limit; config in `.clang-format`)
  is enforced on all C++ files.
- clang-tidy with an extensive check list (`.clang-tidy`); run via `make check`
  (needs `compile_commands.json`, e.g. `compiledb -n --full-path make TESTMODE=ON`).
- cppcheck runs in CI (`.github/workflows/ci-checks.yml`).
- Pre-commit hooks (`.pre-commit-config.yaml`, managed with `prek`/pre-commit)
  run clang-format, whitespace/EOF fixers, and a `make OPTIMIZE=OFF` compile.

Before committing, at minimum run clang-format on changed files and make sure
`make OPTIMIZE=OFF` succeeds.

## Conventions and gotchas

- Output files (e.g. `linestat.out`, `*.out` spectra/light curve files) are
  parsed downstream by the artistools Python package. Do not change their
  formats or make them optional — `linestat.out` in particular must remain
  byte-compatible.
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
- The default branch is `main`; `release` is the production branch. CI runs
  on every push.
