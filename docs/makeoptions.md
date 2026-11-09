# Make options

The Makefile is the authority for these options and for the way that they interact.

| Option | Effect |
| --- | --- |
| `TESTMODE=ON` | Enable additional assertions and the address and undefined behaviour sanitizers. |
| `FASTMATH=OFF` | Do not enable compiler transformations that affect round-off-level results, e.g. a\*(b\*c) = (a\*b)\*c. |
| `MAX_NODE_SIZE=N` | Artificially limit the MPI node size to N ranks. This is useful for testing and prevents MPI shared memory windows from crossing CPU sockets. |
| `REPRODUCIBLE=ON` | Use stable sorts and disable `FASTMATH`. |
| `GPU=ON` | Required to compile for GPUs. It works around incompatible function calls and uses a Simpson-rule integrator in place of Gauss-Kronrod. |
| `GCCTOOLCHAIN=path` | With nvc++, select the gcc toolchain at that path. The Makefile has fallbacks for some known hosts, so this option is not always necessary. |
| `GPUARCH=N` | With nvc++ and `GPU=ON`, compile for compute capability N, e.g. `GPUARCH=80`, in place of the GPU of the host. A host that has no GPU needs this option to get a defined target. |
| `OPENMP=ON` | Parallelise with OpenMP threads within each MPI rank. This cannot be combined with `STDPAR`. The Makefile gives the OpenMP flags to gcc, clang, and nvc++ only, so a build with hipcc gets no threads. |
| `STDPAR=ON` | Give the loops to the parallel algorithms of the standard library, on the CPU cores or, together with `GPU=ON`, on a GPU. The Makefile gives nvc++ the flag `-stdpar=multicore` and links `-ltbb` for gcc, so these two compilers use the CPU cores. With clang or hipcc, the threads depend on the parallel backend of the standard library, which can run the algorithms in one thread. gcc needs the TBB library (see [Installation](installation.md)). |
| `OPTIMIZE=OFF` | Compile without optimisation. This is the quickest way to check that everything still compiles. |
| `ZSTD=OFF` | Build without libzstd. |
| `PGO=GENERATE`, `PGO=USE` | Profile-guided optimisation with gcc or clang. Build with `GENERATE`, run a representative simulation to collect profile data, then build again with `USE`. |

The Makefile finds libzstd with pkg-config or with a link of `-lzstd`. A build with the library reads a zstd compressed
input file, e.g. `model.txt.zst`, when the plain file is absent, and it writes the output files compressed (see
[Input files](inputfiles.md)). `ZSTD=ON` stops the build if the library is absent.

Each option of the table above adds a suffix to the name of the build folder, so different option sets coexist in
`build/`. `GCCTOOLCHAIN` is the exception. It changes `CXXFLAGS` but not the name of the folder, so make keeps the
objects of the previous toolchain. Run `make clean` after you change its value.

The name of the build folder ends with `_options` and a checksum of the content of `artisoptions.h`, so the builds of
different {term}`presets <preset>` coexist too. A change of the preset or of an option value compiles every object in a new build
folder, and a change back to an earlier content compiles nothing again. The timestamp of `artisoptions.h` has no
effect. Each edit of an option keeps the build folder of the old value, so run `make clean` when `build/` becomes
large.

Do not change `artisoptions.h` during a build. Make then stops with an error and removes the object, because that
object does not belong in the build folder. Run `make` again.

## The targets

- `make` builds `sn3d` and `exspec`.
- `make unittests` builds the unit tests (see [Testing](testing.md)).
- `make check` runs clang-tidy over all three programs.
- `make compile_commands.json` writes the compilation database for the clang tools.
- `make sn3dwhole` builds `sn3d` with one compiler command over all the sources, in place of one object file for each source. The `sn3d` symlink then points at that binary.
- `make clean` removes all build folders and the three program symlinks.
- `make help` prints one line for each target and each option. It needs no `artisoptions.h` and no compiler.
