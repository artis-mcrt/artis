# AGENTS.md

Guidance for the AI coding agents that work on ARTIS. See [README.md](README.md)
for the full description of the project.

This file names where a fact lives. It does not copy the fact. If an exact
value matters, read it from the Makefile, from `.github/workflows/`, or from
the source.

## Writing style

Write all English in ASD-STE100 (Simplified Technical English). This applies to:

- code comments;
- documentation;
- commit messages;
- pull request text;
- replies to the maintainers;
- new log, warning, and error strings in the code.

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
- The format and the content of the output files (see "Code conventions").
- A log string that a script reads. `sn3d` writes `RESTART_NEEDED`, and the job
  scripts in `scripts/` search the log for that text. Do not change such a
  string.
- Quoted text from an external source, e.g. a compiler message or a title of a
  publication.

Keep the British spellings that this repository uses, e.g. "parallelised" and
"normalise". STE controls the choice of words and the structure of the
sentences. It does not control the spelling variant here.

## Project overview

ARTIS is a 3D Monte Carlo radiative transfer code for supernovae and kilonovae.
The code uses modern C++ and the Message Passing Interface (MPI). You can also
add OpenMP threads or C++ standard parallelism. gcc and clang compile the code
as C++26. nvc++ and hipcc compile it as C++23, so the code must stay compatible
with C++23. The Makefile builds three programs from the same sources:

- `sn3d`: the main simulation. It does the radiative transfer over a sequence
  of timesteps.
- `exspec`: it combines the packet files into spectra and light curves.
- `unittests`: the unit tests. They cover the numeric and the parsing helpers,
  and also some physics functions, e.g. the Compton cross-section and the
  Spencer-Fano solver.

## Repository layout

- `*.cc` and `*.h` in the repository root: all of the C++ source. `sn3d.cc`,
  `exspec.cc`, and `unittests.cc` contain the three `main()` functions. The
  Makefile compiles every other `.cc` file into all three programs. That list
  is the variable `common_files` in the Makefile. Add a new source file to it,
  because the Makefile finds no file automatically.
- `artisoptions_*.h`: the presets of the compile-time options. The build needs
  an `artisoptions.h` in the repository root. That file is gitignored and is
  normally a symlink to one preset. `artisoptions_doc.md` documents the options
  in a single fenced C++ block.
- `data/`: the physics data that comes with the code, e.g. the decay data and
  the gamma spectra.
- `tests/`: the setups of the end-to-end tests (`setup_*.sh` and the related
  `*_inputfiles/` folders).
- `scripts/`: the job scripts for HPC clusters and some helper tools.
- `third_party/`: a pruned copy of Eigen (see `third_party/eigenversion.txt`).
  Do not edit it and do not lint it.
- `pixi.toml`: a conda-forge environment with Open MPI and gcc. `pixi run make`
  is the quickest way to get a toolchain that works.

The repository has no `CONTRIBUTING.md`, no `docs/` folder, and no template for
a pull request or an issue.

## Build

The build needs gcc 14 or later, clang, nvc++, or hipcc, and an MPI library
that supplies `mpicxx`. With Open MPI, select the compiler with
`export OMPI_CXX=g++`. The build needs no other program.

```sh
export MAKEFLAGS="--check-symlink-times --jobs=$(nproc --all)"
ln -s artisoptions_classic.h artisoptions.h   # required, gitignored
make                                          # builds sn3d and exspec
```

Symlink the preset instead of a copy. The symlink keeps `artisoptions.h`
correct when the preset changes, e.g. after a change of branch. Make reads the
timestamp through the symlink without `--check-symlink-times`, so a change of
preset then starts no rebuild. The Makefile warns you if `artisoptions.h` is a
symlink and the flag is absent. A copied `artisoptions.h` gets no such warning.

README.md lists the make options, and the Makefile is the authority. `make
OPTIMIZE=OFF` is the quickest test that everything compiles. The other targets
are `unittests`, `check` (clang-tidy), `sn3dwhole`, and `clean`.

`make` puts the objects in `build/<configuration>/` and makes `sn3d`, `exspec`,
and `unittests` in the repository root symlinks to that folder. Each symlink
points to the most recent build of any configuration. Look at the symlink
before you copy a binary or start a run. Copy a binary with `cp -L` from the
repository root. A stale build folder of an older compiler otherwise gives a
benchmark that does not measure your change.

Each build also writes `compile_commands.json` (see "Linting and formatting").

Points that surprise a new agent:

- Each option adds a suffix to the name of the build folder, so different
  option sets coexist. A new option that changes `CXXFLAGS` must also add a
  suffix to `BUILD_DIR`.
- The name of the build folder does not contain the preset. A change of preset
  therefore needs `make clean`, and CI does this before nearly every preset
  compile.
- The Makefile decides how the options interact, e.g. `REPRODUCIBLE=ON` with
  `FASTMATH`, and `OPENMP=ON` with `STDPAR=ON`. Read the Makefile for the rule.
- The build uses `-Werror` for all compilers except nvc++.
- clang adds `-Wunsafe-buffer-usage`, which gcc does not have. A build that is
  clean with gcc can still fail with clang.
- `TESTMODE=ON` adds the sanitizers, the extra assertions, and the hardened
  modes of the standard library. A run then needs
  `ASAN_OPTIONS=detect_stack_use_after_return=1:detect_leaks=0` and, for Open
  MPI, `OMPI_MCA_memory=^patcher`.
- `GCCTOOLCHAIN` selects the gcc toolchain of nvc++. The Makefile has fallbacks
  for some known hosts, so this option is not always necessary. CI sets it in
  the environment.
- The Makefile writes `version.cc` from the git metadata, but only when the
  content changes. `version.cc` is gitignored, and `version.h` is a tracked
  file that declares the three strings. Keep the definitions in `version.cc`.
  The git metadata changes at almost every git operation, so a definition in
  the header makes `sn3d.cc` and `exspec.cc` compile again each time.

## Tests

`make OPTIMIZE=OFF unittests && ./unittests` builds and runs `unittests.cc`.
The `OPTIMIZE=OFF` build folder already holds the objects of the compile test,
so only `unittests.cc` compiles again. The harness is hand-written, so there is
no gtest and no catch. Read `unittests.cc` to see if it calls a function that
you changed. The `constexpr` helpers also have `static_assert` checks next to
their definitions.

The main coverage is the end-to-end tests in `tests/`. Each test does a new run
and then a resume run, and then post-processes the packets with `exspec`.
Continuous integration (CI) compares md5 checksums of the output files with the
reference files in `tests/<testname>_inputfiles/`. Each test has two of them:
`results_md5_job0.txt` and `results_md5_final.txt`. CI is the test bed. Commit
and push early, and let CI run the full matrix.

To repeat one test locally, do the steps of the test job in
`.github/workflows/ci.yml`. Use `nebular_1d_3dgrid` as the default test. It
builds the nebular preset with the NLTE populations, the non-thermal solver,
and the detailed bound-free estimators on a 3D grid, so it covers the physics
that most changes touch. The `kilonova_1d` test is faster, but it uses the LTE
preset. The setup script downloads about 15 MB of atomic data from a GitHub
release, so this step needs the network. The script keeps the archive in
`tests/`, and a later run of the script uses that copy.

Three steps differ from a plain build and are easy to miss:

- Build from the `artisoptions.h` of the run folder, not from the preset. Each
  setup script copies a preset and then changes some option values with `sed`.
  Remove your `artisoptions.h` before the copy, because `cp` writes through a
  symlink and replaces the content of the tracked preset.
- Copy `input-resume.txt` to `input.txt` before the resume run. `sn3d` restores
  an absent `input.txt` from `input-newrun.txt`, but the resume run needs the
  restart state that the first run wrote.
- Remove the `*.tmp` files before `exspec`, and the `*_res_*.out` files after
  `mergeangleres.py`, as the workflow does.

CI writes `results_md5_job0.txt` from
`md5sum *.out job0/*.out speclc_angle_res/*.*` and `results_md5_final.txt` from
the same command with `job1/`. The log files stay outside both sets, because
their names do not match.

CI makes the reference checksums on an arm64 runner with g++-15. Local x86-64
builds with gcc 14 reproduced all of them for `kilonova_1d`, on two different
CPU types. `REPRODUCIBLE=ON` therefore gives portable results, but only these
combinations have a test. Examine a local mismatch as a real change of the
results first. Let CI give the decision.

A change that must not alter the results must give identical checksums in CI.
If a change alters the numerical results for a good reason, the stored
checksums must be regenerated. Only a maintainer can do this, with the
`updatechecksums.yml` workflow. That workflow takes the checksums from a
finished CI run. The run does not have to be of the head commit, if the later
commits do not change the results. Say in the commit message and in the pull
request that you expect the results to change.

A new end-to-end test must run in CI in about ten minutes or less, because
every push runs the full matrix.

## Continuous integration

The workflows in `.github/workflows/`:

- `ci.yml` runs each model of its `testname` matrix on arm64, and compares the
  checksums. Two more jobs run one model with `OPENMP=ON` and with
  `STDPAR=ON`. These two jobs compare no checksums (see "Reproducible
  results").
- `cislowtestmode.yml` calls `ci.yml` again with `TESTMODE=ON`. The sanitizers
  and the extra assertions make this run slow, so `ci.yml` gives it a longer
  timeout. This workflow does not enforce the checksums, because both checksum
  steps then continue after an error.
- `ci-checks.yml` runs the pre-commit hooks, clang-tidy, and cppcheck. It then
  compiles the code with each compiler of its matrix, and on macOS, with
  hipcc, and with nvc++. The GPU compilers build the classic and the nebular
  presets with `STDPAR=ON GPU=ON`. The gcc, the clang, and the macOS jobs
  compile every remaining preset. The gcc and the clang jobs also build and run
  the unit tests, for the classic and for the nebular preset.
- `updatechecksums.yml` writes the reference checksums (see "Tests").
- `depapprove.yml` enables auto-merge for the pull requests of Dependabot and
  of pre-commit-ci.
- `copilot-setup-steps.yml` runs only when you change that file.

The first three run on each push, except on a `classic*` branch. Your code
must therefore compile with each of these compilers, and also on the paths for
`GPU=ON`. A new compile-time option must go into every `artisoptions_*.h`
preset, because CI finds the presets with a glob.

CI also installs the artistools Python package and plots the output of each
test (see "Input and output files").

## Linting and formatting

- clang-format uses a Google-based style with a limit of 120 columns
  (`.clang-format`). It applies to all C++ files except `third_party/`.
- clang-tidy has a long list of checks (`.clang-tidy`) and makes almost every
  diagnostic an error. Make the compile database with the Makefile, then name
  the files that you changed:

  ```sh
  make TESTMODE=ON compile_commands.json
  run-clang-tidy grid.cc
  ```

  Each build of a program and the `make` pre-commit hook also write the
  database. The recipe needs python3 and replaces the file only when the
  content changes. A host without python3 gets a message and a normal build.
  The database always has the flags of `TESTMODE=ON`, also when the build does
  not, so that clangd and clang-tidy examine the body of each
  `assert_testmodeonly()`. `make TESTMODE=ON check` runs the same check as CI
  over all three programs and takes many minutes. Never make the database with
  `compiledb`, because it drops `mpicxx` from each command and every diagnostic
  is then false.
- cppcheck runs in CI and stops the job on an error. Read the command in
  `ci-checks.yml`: it excludes `third_party`, and it suppresses some message
  types that the code does not correct. Write a suppression as an inline
  `// cppcheck-suppress <id>` comment. That comment works because the command
  has `--inline-suppr`.
- The pre-commit hooks (`.pre-commit-config.yaml`, run with `prek`) apply
  clang-format, about a dozen whitespace and file checks, and a
  `make OPTIMIZE=OFF` compile. `prek run --all-files` runs all of them. CI runs
  the same hooks but skips the compile. One check finds a destroyed symlink,
  which is important here: `CLAUDE.md` is a symlink to `AGENTS.md`.

## Code conventions

### Names

The readers of this code are scientists, not only programmers. Give each
function, variable, and type a name that says which physical quantity it holds.
Do not name a thing after its role in an abstract algorithm.

- Name the quantity, not the position in a procedure. Write
  `get_log_te_nne_ionpops()` and not `get_outer_state()`. The first name says
  that the result holds the logarithms of T_e, of nne, and of the ion
  populations. The second name says only that some state belongs to some outer
  loop, which a reader cannot check against the physics.
- Use the symbol of the quantity if the code already uses it, e.g. `T_e`, `nne`,
  `nnion`, `phi`, and `Gamma`. Keep the case of the existing symbol.
- A single letter is permitted only for a mathematical object that has no
  physical value, e.g. the index `i` of a loop. Do not give a single letter to a
  vector of physical quantities.
- Give an index variable the name of its index space (see "Cell indices").
- Expand an abbreviation that this repository does not already use. `mgi`,
  `nne`, and `bf` are established here. A new short form is not.
- The name must stay correct for a reader who does not know the algorithm. A
  comment that explains the name is a sign that the name is wrong.

### Comments

- A comment explains the physics, the algorithm, or the reason for a choice.
  Cite the journal article where one exists. Do not repeat the code.
- A better name is better than a comment. Rename the variable or the function
  instead.
- Describe the current code only. The version control system holds the history.

### Cell indices

ARTIS uses three different index spaces for the cells, and all three are plain
integers. A mix of two of them compiles cleanly and reads the wrong cell. The
comment above `get_propcell_modelgridindex()` in `grid.h` defines them:

- `cellindex` selects a cell of the propagation grid.
- `modelgridindex` (`mgi`) selects a cell of the input model.
- `nonemptymgi` selects a model cell that contains matter. The per-cell physics
  arrays have only these entries, so the physics modules use this index.

Convert between the last two with `grid::get_nonemptymgi_of_mgi()` and
`grid::get_mgi_of_nonemptymgi()`. Do not call
`grid::get_nonemptymgi_of_mgi()` for an empty model cell. It has an assertion
for this, which only a `TESTMODE=ON` build finds.

`get_propcell_modelgridindex()` and `get_propcell_nonemptymgi()` give -1 for an
empty cell. Test for -1 before you use the result as an index. These two
functions are the only correct test for an empty cell.

Give each variable the name of its index space, e.g. `nonemptymgi` for a
`nonemptymgi`. The compiler cannot find a mix of two index spaces, so the names
are the only defence.

### Logs and assertions

- Write log lines with `printlog()` and `printlnlog()` from `mpi_logging.h`.
  They take a `std::format` string. Do not use `printf`, `std::cout`, or
  `std::cerr`.
- Every log line gets a timestamp and a flush, so a log line in a hot loop is
  expensive.
- Do not call the loggers in a `DEVICE_FUNC`. Their device branch does not
  compile with nvc++ or hipcc, and a local gcc or clang build does not find
  this. Put a diagnostic of a device function inside `MY_IF_HOST()`.
- `assert_always()` stays active in an optimised build, and its reporter works
  on the host and on the device. `assert_testmodeonly()` becomes nothing unless
  `TESTMODE=ON`, so its expression must have no side effect. A parameter that
  only an `assert_testmodeonly()` uses needs `[[maybe_unused]]`.
- Many `assert_always()` calls hold a call with a side effect, e.g. a read
  from a file. Do not change such an assertion to `assert_testmodeonly()` and
  do not delete it.
- A fatal error is an `assert_always()` on the condition. Write an optional
  `printlnlog("[error] ...")` with the values before it. Do not call
  `std::abort()` directly. A direct abort is very hard to debug, because the
  assertion reporter adds the file, the line, the expression, and the function,
  and it does this on the host and on the device.
- An invalid input stops the run. Do not give a warning for it.
- Two numeric helpers are an exception: `toms748.h` and `gausskronrod.h` throw
  `std::domain_error` on the host, inside a guard that returns a NaN for a GPU
  build. No code catches these exceptions. Keep the guards, because device code
  permits no exception. Both files are bit-exact extractions of Boost.Math. Do
  not change their floating-point expressions or their table literals.

### Global and shared state

- `namespace globals` holds the global simulation state. Two headers open it:
  `globals.h` and `mpi_logging.h`, which has the MPI ranks and communicators.
- Large read-mostly arrays are `MPI_shared_array` members, so all ranks on a
  node share one copy. Search for `MPI_shared_array` to find them. Some are
  outside `namespace globals`.
- The construction of a shared array is collective over the node. Every rank on
  the node must call it. Only `rank_in_node == 0` asks for the memory, but the
  other ranks must take part, or the run stops. Do not put a construction
  inside `if (globals::rank_in_node == 0)`.
- Each rank then fills its own slice of the new array, so that the operating
  system puts the pages near that rank. The constructor ends with a barrier
  over the node, so you do not add one.
- In the grid update, each rank changes only the cells of its own range. See
  `update_grid.cc` and the assignment of the ranks in `grid.cc`. After the grid
  update, `sn3d.cc` broadcasts the per-cell state to the other nodes. A new
  per-cell array that the packets read must join that broadcast, or the other
  nodes propagate with stale values and give no error.
- The propagation of the packets is different. Each rank adds to the estimators
  of every cell that its packets enter. A new estimator therefore needs a sum
  over the ranks with `MPI_Allreduce_safe()`, as `radfield.cc` does. Without
  that sum the code loses the results of the other ranks and gives no error.
- Use the wrappers `MPI_Allreduce_safe()`, `MPI_Bcast_safe()`, and
  `MPI_Reduce_safe()`. They map the types and divide a large message into
  chunks, because MPI counts are 32-bit.

### Reproducible results

- Sort with the macro `SORT_OR_STABLE_SORT` if the order of two equal elements
  can change a result. The macro gives a stable sort when `REPRODUCIBLE=ON`.
  Some sorts use `std::ranges::sort` directly, because their key is unique.
  Keep them as they are.
- A multithreaded run is never reproducible, because each thread seeds its own
  generator and the accumulation order changes.
- Do not regroup a floating-point sum for elegance. Some operation orders are
  frozen on purpose, e.g. the panel structure in `solve_upper_triangular()` in
  `nonthermal.cc`.
- The random generator is a Xoshiro128PP in `random.h`. Take numbers with
  `rng_uniform(get_rngstate(pkt))` and `rng_uniform_pos(get_rngstate(pkt))`.
  `get_rngstate()` hides the difference between the builds. The host keeps one
  generator for each thread. A GPU build keeps the state inside each `Packet`,
  because device code cannot use `thread_local`.

### Portability to GPUs

The code must compile with nvc++ and with hipcc, also with `STDPAR=ON GPU=ON`.

- Mark a function that device code calls with `DEVICE_FUNC`.
- Separate a host path from a device path with `MY_IF_HOST()` and
  `MY_IF_DEVICE()`.
- Device code permits no exception and no `thread_local` variable. Use
  `THREADLOCALONHOST` for a variable that is thread-local on the host only.
- Give a parallel algorithm the macro `EXEC_PAR` or `EXEC_PAR_UNSEQ`.
- Accumulate into shared counters with `atomicadd()`, `atomicload()`, and
  `atomicstore()`. Add `ALIGNAS_AVOID_FALSE_SHARING` to a shared counter. For a
  lock, use `PaddedMutex` and `ScopedMutex`, not `std::mutex`.
- Some workarounds look like unnecessary complexity but are necessary. An
  example is `lowest_set_bit()` in `constants.h`, because nvc++ cannot
  translate `std::countr_zero` for a device.

### Input and output files

- The artistools Python package reads the output files. Do not change their
  formats. Do not make them optional. `linestat.out` must stay
  byte-compatible. CI plots each test with artistools, so a change of a format
  fails CI there.
- The numeric values of some enumerations are also part of that interface, e.g.
  `packet_type`, `absorption_type`, and the `EMTYPE_*` constants in `packet.h`.
  Do not renumber them.
- `sn3d` writes one log file for each rank and thread
  (`output_<rank>-<thread>.txt`). The option `-o JOBFOLDER` moves the per-job
  files into a subfolder. The run-level files, e.g. the restart files, stay in
  the run folder, together with a symlink to the log of rank 0. The standard
  output stays quiet unless there is a crash.
- The restart files (`gridsave_ts*.tmp` and the packet files) must only be
  consistent with the binary that wrote them. A resumed run uses the same
  `artisoptions.h` and the same source version as the run that wrote the
  files. Do not add a format version or a fallback for a restart file of a
  different option set or of an older binary.
- `input.txt` is positional. Each line has its fixed meaning, and some lines
  are unused placeholders that must stay. `read_parameterfile()` reads a line
  with `get_noncommentline()` and then takes the numbers with a
  `std::istringstream`. Follow that pattern for a new line. A new run writes
  a commented copy of `input.txt` to `input-newrun.txt`. If `input.txt` is
  absent at the start of a run, rank 0 restores it from that copy.
- The model files use a different helper. `model.txt`, `abundances.txt`, and
  `transitiondata.txt` take each number with `parse_next_token()`, which
  advances a `std::string_view`.
- Open a file with `fopen_required()` or `fstream_required()`. For a read they
  look in `./`, `data/`, and `artis/data/`. For a write they use the name that
  you give.

### C++ style

- Prefer modern C++: `constexpr`, `std::ranges`, `std::span`, and structured
  bindings. Put a compile-time option into `artisoptions.h` as a `constexpr`
  value, not into a runtime flag. Prefer a boolean option when only one other
  value is useful.
- Write a trailing return type: `auto f(...) -> T`. clang-tidy requires it.
- Give each `.cc` file an anonymous namespace for its internal helpers. Close
  it with the comment `}  // anonymous namespace`. Every tracked `.cc` file
  does this.
- Some modules put their interface in a namespace, e.g. `grid::` and `decay::`.
  Other modules declare their functions at global scope. Follow the header of
  the module that you change.
- Use `#ifndef` header guards, not `#pragma once`.
- clang-format sorts the includes. The header of the module comes first, then
  the system headers, then the project headers.
- clang gives `-Wunsafe-buffer-usage` for some system headers and for the
  construction of a `std::span` from a raw pointer. Enclose the `#include` or
  the construction in `#pragma clang unsafe_buffer_usage begin` and `end`. Most
  of the existing pragmas enclose an `#include`.

## Pull requests

- Make one pull request for one work item. Do not combine an unrelated fix
  with a feature.
- Commit and push early. CI is the test bed for the compilers and the
  checksums that a local build does not cover.
- Say in the commit message and in the pull request when you expect the
  numerical results to change. Otherwise the checksums must stay identical.

## Before you commit

Do the text edits first. They can make a finished compile out of date.

1. Add a new compile-time option to every `artisoptions_*.h` preset and to
   `artisoptions_doc.md`.
2. Search `tests/setup_*.sh` for the name of an option that you renamed or
   reformatted. Those scripts change option lines with `sed` and an exact text
   match. The pattern also contains the type, e.g. `constexpr int`. A pattern
   that matches nothing gives no error, so the test then runs with the default
   value of the preset and the checksums drift.
3. Run `prek run --all-files` and correct every message.
4. Run `make OPTIMIZE=OFF unittests && ./unittests` if `unittests.cc` calls a
   function that you changed.
5. Run clang-tidy on the files that you changed. CI stops on any diagnostic,
   and a gcc build does not find it.
6. Check that no `DEVICE_FUNC` calls a logger (see "Logs and assertions").

The default branch is `main`, and `release` is the production branch.
