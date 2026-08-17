# AGENTS.md

Guidance for the AI coding agents that work on ARTIS. See [README.md](README.md)
for the full description of the project.

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
  `*.out` files must stay byte-compatible (see "Code conventions").
- Quoted text from an external source, e.g. a compiler message or a title of a
  publication.

Keep the British spellings that this repository uses, e.g. "parallelised" and
"normalise". STE controls the choice of words and the structure of the
sentences. It does not control the spelling variant here.

## Project overview

ARTIS is a 3D Monte Carlo radiative transfer code for supernovae and kilonovae.
It is written in modern C++ and parallelised with MPI. OpenMP threads or C++
standard parallelism can be added. gcc and clang compile the code as C++26.
nvc++ and hipcc compile it as C++23, so the code must stay compatible with
C++23. Three programs are built from the same sources:

- `sn3d`: the main simulation. It does the radiative transfer over a sequence
  of timesteps.
- `exspec`: it combines the packet files into spectra and light curves.
- `unittests`: the unit tests of the numeric and the parsing helpers.

## Repository layout

- `*.cc` and `*.h` in the repository root: all of the C++ source. `sn3d.cc`,
  `exspec.cc`, and `unittests.cc` contain the three `main()` functions. The
  Makefile compiles every other `.cc` file into all three programs, so a new
  file in the root is picked up automatically.
- `artisoptions_*.h`: six presets of the compile-time options. The build needs
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
`export OMPI_CXX=g++`.

```sh
ln -s artisoptions_classic.h artisoptions.h   # required, gitignored
make -j$(nproc)                               # builds sn3d and exspec
```

Symlink the preset instead of a copy. The symlink keeps `artisoptions.h`
correct when the preset changes, e.g. after a change of branch. Add
`--check-symlink-times` to your make flags, because make reads the timestamp
through the symlink without it. The Makefile prints a warning if you forget
this flag.

README.md lists all of the make options. `make OPTIMIZE=OFF` is the quickest
test that everything compiles, and the pre-commit hook runs it. The other
targets are `unittests`, `check` (clang-tidy), `sn3dwhole`, and `clean`.

Points that surprise a new agent:

- Each option adds a suffix to the name of the build folder. Different
  configurations therefore coexist, and `make clean` is rarely necessary. A new
  option that changes `CXXFLAGS` must also add a suffix to `BUILD_DIR`.
- `REPRODUCIBLE=ON` sets `FASTMATH=OFF` by itself.
- `OPENMP=ON` and `STDPAR=ON` together are an error.
- The build uses `-Werror` for all compilers except nvc++.
- clang adds `-Wunsafe-buffer-usage`, which gcc does not have. A build that is
  clean with gcc can still fail with clang.
- `TESTMODE=ON` adds the sanitizers, the hardened standard library modes, and
  the extra assertions. A run then needs
  `ASAN_OPTIONS=detect_stack_use_after_return=1:detect_leaks=0` and, for Open
  MPI, `OMPI_MCA_memory=^patcher`.
- nvc++ needs `GCCTOOLCHAIN=<path to g++>`. README.md does not list this
  option.
- The Makefile writes `version.h` from the git metadata, but only when the
  content changes. `version*.h` is gitignored.

## Tests

`make unittests && ./unittests` builds and runs the unit tests of the numeric
and the parsing helpers in `unittests.cc`. The harness is hand-written, so
there is no gtest and no catch. The `constexpr` helpers also have
`static_assert` checks next to their definitions.

The main coverage is the end-to-end tests in `tests/`. Each test does a new run
and then a resume run, and then post-processes the packets with `exspec`. CI
compares md5 checksums of the output files with the reference files
`results_md5_job0.txt` and `results_md5_final.txt`.

To repeat one test locally, do the same steps as `.github/workflows/ci.yml`:

```sh
cd tests
source ./setup_kilonova_1d.sh          # creates tests/kilonova_1d_testrun/
cd ..
cp tests/kilonova_1d_testrun/artisoptions.h .   # do not skip this step
make REPRODUCIBLE=ON MAX_NODE_SIZE=2 FASTMATH=OFF -j$(nproc) sn3d exspec
cp sn3d exspec tests/kilonova_1d_testrun/
cd tests/kilonova_1d_testrun
cp input-newrun.txt input.txt
mpirun -np 4 --oversubscribe ./sn3d -o job0    # logs go to job0/output_0-0.txt
md5sum -c results_md5_job0.txt
cp input-resume.txt input.txt
mpirun -np 4 --oversubscribe ./sn3d -o job1
rm *.tmp
mpirun -np 1 ./exspec
python3 ../../scripts/mergeangleres.py
rm -f light_curve_res_*.out spec_res_*.out specpol_res_*.out
md5sum -c results_md5_final.txt
```

The copy of `artisoptions.h` is necessary. Each setup script copies a preset
into the run folder and then changes some option values with `sed`. A build
from the preset alone uses different options and gives different results.

CI collects the checksums with
`md5sum *.out job0/*.out speclc_angle_res/*.*`. The log files stay outside this
set, because their names do not match.

CI makes the reference checksums on an arm64 runner with g++-15, but
`REPRODUCIBLE=ON` makes them portable. A local x86-64 build with gcc 14
reproduces them exactly. A mismatch is therefore a true change of the results,
and not an effect of your machine.

If a change alters the numerical results for a good reason, the stored
checksums must be regenerated. Only a maintainer can do this, with the "Update
checksums" workflow. That workflow needs a finished CI run for the exact head
commit. Say in the commit message and in the pull request that you expect the
results to change.

## Continuous integration

Every push runs three workflows. Pushes to `classic*` branches are the only
exception.

- `ci.yml` runs ten end-to-end test models on arm64 with g++-15, and compares
  the checksums. Two more jobs run one model with `OPENMP=ON` and with
  `STDPAR=ON`. These two jobs compare no checksums, because a multithreaded run
  is not reproducible.
- `cislowtestmode.yml` runs the same ten models again with `TESTMODE=ON`. This
  gives the sanitizers and the extra assertions, and it is slow. Its timeout is
  360 minutes.
- `ci-checks.yml` runs the pre-commit hooks, clang-tidy, and cppcheck. It then
  compiles the code with gcc 14, 15, and 16, with clang, with Apple Clang and
  LLVM clang on macOS, with hipcc, and with nvc++. It compiles the classic and
  the nebular presets with `OPENMP=ON`, with `STDPAR=ON`, and with
  `STDPAR=ON GPU=ON`, and then all of the other presets. It also builds and
  runs the unit tests.

A new compile-time option must therefore go into all six preset files and into
`artisoptions_doc.md`. The code must compile with all of these compilers,
including the paths for `GPU=ON`.

CI also installs the artistools Python package from PyPI and plots the output
of each test. A change to the format of an output file can break CI through
that package.

## Linting and formatting

- clang-format uses a Google-based style with a limit of 120 columns
  (`.clang-format`). It applies to all C++ files except `third_party/`.
- clang-tidy has a long list of checks (`.clang-tidy`) and treats every
  diagnostic as an error. Run it with `make check`, which needs a
  `compile_commands.json`, e.g. from
  `compiledb -n --full-path make TESTMODE=ON`.
- cppcheck runs in CI with `--enable=all --inconclusive --check-level=exhaustive`
  and stops the job on any message. Write a suppression as an inline
  `// cppcheck-suppress <id>` comment.
- The pre-commit hooks (`.pre-commit-config.yaml`, run with `prek`) apply
  clang-format, some whitespace and file checks, and a `make OPTIMIZE=OFF`
  compile. They also run after a checkout, a merge, and a rebase, so the
  compile hook starts again at those times. CI skips the compile hook.

## Code conventions

### Cell indices

ARTIS uses three different index spaces for the cells, and all three are plain
integers. A mix of two of them compiles cleanly and reads the wrong cell. The
comment above `get_propcell_modelgridindex()` in `grid.h` defines them:

- `cellindex` selects a cell of the propagation grid.
- `modelgridindex` (`mgi`) selects a cell of the input model.
- `nonemptymgi` selects a model cell that contains matter. The per-cell physics
  arrays have only these entries, so the physics modules use this index.

Convert between the last two with `grid::get_nonemptymgi_of_mgi()` and
`grid::get_mgi_of_nonemptymgi()`. `get_propcell_modelgridindex()` and
`get_propcell_nonemptymgi()` give -1 for an empty cell. Test for -1 before you
use the result as an index.

### Logs and assertions

- Write log lines with `printlog()` and `printlnlog()` from `mpi_logging.h`.
  They take a `std::format` string. Do not use `printf`, `std::cout`, or
  `std::cerr`. The function `printout()` of the older code does not exist any
  more, although `.clang-tidy` still names it.
- Every log line gets a timestamp and a flush, so a log line in a hot loop is
  expensive.
- `assert_always()` stays active in an optimised build. `assert_testmodeonly()`
  becomes nothing unless `TESTMODE=ON`, so its expression must have no side
  effect. A parameter that only an `assert_testmodeonly()` uses needs
  `[[maybe_unused]]`.
- A fatal error writes `printlnlog("[error] ...")` and then calls
  `std::abort()`. The host path has no exception-based error handling.

### Global and shared state

- `namespace globals` holds the global simulation state. It is split between
  `globals.h` and `mpi_logging.h`, which has the MPI ranks and communicators.
- Large read-mostly arrays are `MPI_shared_array` members, so all ranks on a
  node share one copy. Some of them are outside `namespace globals`, e.g. in
  `grid.h`, `kpkt.h`, `nltepop.h`, and `rpkt.h`.
- Only the rank with `rank_in_node == 0` allocates and fills a shared array.
  A barrier over the node must follow the write.
- In a timestep, each rank updates only the cells of its own range. See
  `update_grid.cc` and the assignment of the ranks in `grid.cc`.
- Use the wrappers `MPI_Allreduce_safe()`, `MPI_Bcast_safe()`, and
  `MPI_Reduce_safe()`. They map the types and divide a large message into
  chunks, because MPI counts are 32-bit.

### Reproducible results

- Sort with the macro `SORT_OR_STABLE_SORT` and not with `std::sort`. The macro
  gives a stable sort when `REPRODUCIBLE=ON`.
- A multithreaded run is never reproducible, because each thread seeds its own
  generator and the accumulation order changes.
- Do not regroup a floating-point sum for elegance. Some operation orders are
  frozen on purpose, e.g. the panel structure in `solve_upper_triangular()` in
  `nonthermal.cc`.
- The random generator is a Xoshiro128PP in `random.h`. Take numbers from
  `rng_uniform()` and `rng_uniform_pos()`. The host keeps one generator for
  each thread. A GPU build keeps the state inside each `Packet`, because device
  code cannot use `thread_local`.

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
  formats and do not make them optional. `linestat.out` must stay
  byte-compatible.
- The numeric values of some enumerations are also part of that interface, e.g.
  `packet_type`, `absorption_type`, and the `EMTYPE_*` constants in `packet.h`.
  Do not renumber them.
- `sn3d` writes one log file for each rank and thread
  (`output_<rank>-<thread>.txt`). The option `-o JOBFOLDER` moves these files
  into a subfolder. The standard output stays quiet unless there is a crash.
- `input.txt` is positional. Each line has its fixed meaning, and some lines
  are unused placeholders that must stay. Read lines with
  `get_noncommentline()` and numbers with `parse_next_token()`.
- Open a file with `fopen_required()` or `fstream_required()`. They look in
  `./`, `data/`, and `artis/data/`.

### C++ style

- Prefer modern C++: `constexpr`, `std::ranges`, `std::span`, and structured
  bindings. Put a compile-time option into `artisoptions.h` as a `constexpr`
  value, not into a runtime flag.
- Write a trailing return type: `auto f(...) -> T`. clang-tidy requires it.
- Give each `.cc` file an anonymous namespace for its internal helpers, and
  close it with the comment `}  // anonymous namespace`. Not every module has a
  namespace of its own: `decay`, `gammapkt`, `grid`, `kpkt`, `nonthermal`,
  `radfield`, `stats`, and `vpkt` have one, and the other modules declare their
  functions at global scope.
- Use `#ifndef` header guards, not `#pragma once`.
- clang-format sorts the includes. The header of the module comes first, then
  the system headers, then the project headers.
- Enclose raw pointer arithmetic in
  `#pragma clang unsafe_buffer_usage begin` and `end`.

## Before you commit

1. Run clang-format on each file that you changed.
2. Run `make OPTIMIZE=OFF` and correct all warnings.
3. Run `make unittests && ./unittests` if you changed a numeric or a parsing
   helper.
4. Add a new compile-time option to all six `artisoptions_*.h` presets and to
   `artisoptions_doc.md`.
5. Search `tests/setup_*.sh` for the name of an option that you renamed or
   reformatted. Those scripts change option lines with `sed` and an exact text
   match. A pattern that matches nothing gives no error, so the test then runs
   with the default value of the preset and the checksums drift.
6. Say in the commit message when you expect the numerical results to change.

The default branch is `main`, and `release` is the production branch.
