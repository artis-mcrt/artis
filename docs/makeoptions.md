# Make options

The Makefile is the authority for these options and for the way that they interact.

| Option | Effect |
| --- | --- |
| `TESTMODE=ON` | Enable additional assertions and the address and undefined behaviour sanitizers. |
| `FASTMATH=OFF` | Do not enable compiler transformations that affect round-off-level results, e.g. a\*(b\*c) = (a\*b)\*c. |
| `MAX_NODE_SIZE=N` | Artificially limit the MPI node size to N ranks. This is useful for testing and prevents MPI shared memory windows from crossing CPU sockets. |
| `REPRODUCIBLE=ON` | Use stable sorts and disable `FASTMATH`. |
| `GPU=ON` | Required to compile for GPUs. It works around incompatible function calls and uses a Simpson-rule integrator in place of Gauss-Kronrod. |
| `GPUARCH=N` | With nvc++ and `GPU=ON`, compile for compute capability N, e.g. `GPUARCH=80`, in place of the GPU of the host. A host that has no GPU needs this option to get a defined target. |
| `OPENMP=ON` | Parallelise with OpenMP threads within each MPI rank. This cannot be combined with `STDPAR`. |
| `STDPAR=ON` | Parallelise with C++ standard library parallel algorithms, either on multicore CPUs or, together with `GPU=ON`, on GPUs. |
| `OPTIMIZE=OFF` | Compile without optimisation. This is the quickest way to check that everything still compiles. |
| `ZSTD=OFF` | Build without libzstd. |
| `PGO=GENERATE`, `PGO=USE` | Profile-guided optimisation with gcc or clang. Build with `GENERATE`, run a representative simulation to collect profile data, then build again with `USE`. |

The Makefile finds libzstd with pkg-config or with a link of `-lzstd`. A build with the library reads a zstd compressed
input file, e.g. `model.txt.zst`, when the plain file is absent, and it writes the output files compressed (see
[Input files](inputfiles.md)). `ZSTD=ON` stops the build if the library is absent.

Each option adds a suffix to the name of the build folder, so different option sets coexist in `build/`.

## The targets

- `make` builds `sn3d` and `exspec`.
- `make unittests` builds the unit tests (see [Testing](testing.md)).
- `make check` runs clang-tidy over all three programs.
- `make compile_commands.json` writes the compilation database for the clang tools.
- `make clean` removes the build folder.
