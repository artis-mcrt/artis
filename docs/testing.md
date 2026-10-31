# Testing

## The unit tests

The unit tests cover the numeric and parsing helpers and some physics functions, e.g. the Compton cross-section and the
triangular solve of the Spencer-Fano matrix. Build and run them:

```sh
make unittests && ./unittests
```

CI runs them for the classic and the NLTE nebular presets. The harness is hand-written, so there is no gtest and no
catch.

## The end-to-end tests

The `tests` folder holds eleven small end-to-end test models. `tests/runtest.sh` runs one of them:

```sh
tests/runtest.sh kilonova_1d
```

Run the script without arguments to see the names of the tests and of the steps. The test job of
[ci.yml](https://github.com/artis-mcrt/artis/blob/main/.github/workflows/ci.yml) calls the same script. The script does
these steps:

1. It removes the model folder `tests/kilonova_1d_testrun` of an earlier run. Then it runs
   `tests/setup_kilonova_1d.sh`, which downloads the atomic data of the test and makes the model folder again with the
   input files and an `artisoptions.h` for the test.
2. It builds `sn3d` and `exspec` from that `artisoptions.h` with the make options of CI, and copies them into the model
   folder. During the build, `artisoptions.h` in the repository folder is a symlink to the options of the test. The
   script then moves your `artisoptions.h` back with its timestamp. A later build with other make options then compiles
   nothing again, and a build with the make options of the test compiles all objects again.
3. It starts a new simulation and compares the checksums of the output files with
   `tests/kilonova_1d_inputfiles/results_md5_job0.txt`.
4. It resumes the simulation from the restart files, and then runs `exspec`.
5. It compares the checksums of the output files with `results_md5_final.txt`.

The script stops at the first step that fails. To run only one step, give the name of the step after the name of the
test, e.g. `tests/runtest.sh kilonova_1d resume`. With `TESTMODE=ON` in the environment, the build has the sanitizers
and the extra assertions.

The script needs Open MPI and libzstd. It also calls these commands:

- `zstd` and `md5sum`, because the test input files are zstd compressed;
- `curl`, `tar`, and `xz`, because the setup step downloads and unpacks the atomic data;
- `python3`, for `scripts/mergeangleres.py` of the `exspec` step;
- `perl` and `rsync`, for the setup step.

The script writes the checksums of the run into the model folder, e.g.
`tests/kilonova_1d_testrun/results_md5_job0.txt`.

:::{note}
CI makes the reference checksums on an arm64 runner with g++-16, so a different compiler or CPU type can give
different checksums. Examine a local mismatch first as a real change of the results. Let CI give the decision. For a
change that must not alter the results, a mismatch in CI is an error in the change.
:::

## Continuous integration

CI runs on every push. It builds with `REPRODUCIBLE=ON ZSTD=ON FASTMATH=OFF MAX_NODE_SIZE=2` and compares md5 checksums of the
output files with the reference checksums in the `tests/*_inputfiles` folders. It also compiles every
`artisoptions_*.h` preset with gcc and clang, and the classic and the NLTE nebular presets additionally with Apple
Clang, nvc++, and hipcc, including the GPU code paths. CI installs artistools and plots the output of each test, so a
change of an output format fails CI there.

A change that must not alter the numerical results has to give identical checksums. For such a change, a mismatch is an
error in the change. A change that alters the results for a good reason needs new reference checksums, which the
"Update checksums" workflow writes.
