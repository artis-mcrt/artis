# Testing

## The unit tests

The unit tests cover the numeric and parsing helpers and some physics functions, e.g. the Compton cross-section and the
triangular solve of the Spencer-Fano matrix. Build and run them:

```sh
make unittests && ./unittests
```

CI runs them for the classic and the {term}`NLTE` nebular {term}`presets <preset>`. The harness is hand-written, so there is no gtest and no
catch.

## The end-to-end tests

The `tests` folder holds eleven small end-to-end test models. `tests/runtest.sh` runs one of them:

```sh
tests/runtest.sh kilonova_1d
```

Run the script without arguments to see the names of the tests and of the steps. The test job of
[ci.yml](https://github.com/artis-mcrt/artis/blob/main/.github/workflows/ci.yml) calls the same script. The script does
these steps:

1. It removes the {term}`model folder` `tests/kilonova_1d_testrun` of an earlier run. Then it runs
   `tests/setup_kilonova_1d.sh`, which downloads the atomic data of the test and makes the model folder again with the
   input files and an `artisoptions.h` for the test.
2. It builds `sn3d` and `exspec` from that `artisoptions.h` with the make options of CI, and copies them into the model
   folder. During the build, `artisoptions.h` of the repository folder is a symlink to the options of the test, so start
   no other build at that time. The script then moves your `artisoptions.h` back. The `sn3d` and `exspec` symlinks of
   the repository folder then point at the build of the test.
3. It checks the input with `sn3d -c` (see [Running sn3d](running.md)). Then it starts a new simulation and compares
   the checksums of the output files with `tests/kilonova_1d_inputfiles/results_md5_job0.txt`.
4. It checks the input and the {term}`restart files <restart file>` with `sn3d -c`, and it resumes the simulation from the restart files. It
   then removes the restart files and runs `exspec`.
5. It compares the checksums of the output files with `results_md5_final.txt`.

The script stops at the first step that fails. To run only one step, give the name of the step after the name of the
test, e.g. `tests/runtest.sh kilonova_1d resume`. With `TESTMODE=ON` in the environment, the build has the sanitizers
and the extra assertions. The script then sets the two variables that a run of such a build needs:

```sh
export ASAN_OPTIONS=detect_stack_use_after_return=1:detect_leaks=0
export OMPI_MCA_memory='^patcher'
```

Each test uses one preset and a reduced atomic dataset from a release of this repository.
`tests/atomic_data_release.txt` names the release, and `tests/setupfuncs.sh` holds the download command.

| Tests | Preset | Archive | Elements |
| --- | --- | --- | --- |
| `classicmode_1d_3dgrid`, `classicmode_1d_3dgrid_thermalinput`, `classicmode_3d` | `classic` | `atomicdata_classic.tar.xz` | the classic dataset, with cross-sections in the v1 format |
| `kilonova_1d`, `kilonova_2d*` | `kilonova_lte` | `atomicdata_feconi.tar.xz` | Fe, Co, Ni |
| `nebular_1d_3dgrid*` | `nltenebular` | `atomicdata_feconi.tar.xz` | Fe, Co, Ni |
| `nltephotospheric_dynamic_ion_range_1d_1dgrid` | `nltephotospheric` | `atomicdata_hefeconi_fe_i_to_vii.tar.xz` | He, Fe, Co, Ni, with Fe I to Fe VII |
| `kilonova_1d_timedepnlte` | `kilonova_nlte` | `atomicdata_sryzrlace.tar.zst` | Sr, Y, Zr, La, Ce |

`kilonova_1d` is the quickest test. `nebular_1d_3dgrid` covers the physics that most changes touch, so it is the
default test for a change of the code.

[Installation](installation.md) lists the library and the commands that the script needs. The script writes the checksums of the run into the model folder, e.g.
`tests/kilonova_1d_testrun/results_md5_job0.txt`.

:::{note}
CI makes the reference checksums on an arm64 runner with g++-16, so a different compiler or CPU type can give
different checksums. Examine a local mismatch first as a real change of the results. Let CI give the decision. For a
change that must not alter the results, a mismatch in CI is an error in the change.
:::

## Continuous integration

CI runs on every push, except on a branch whose name starts with `classic`. It builds with `REPRODUCIBLE=ON ZSTD=ON FASTMATH=OFF MAX_NODE_SIZE=2` and compares md5 checksums of the
output files with the reference checksums in the `tests/*_inputfiles` folders. It also compiles every
`artisoptions_*.h` preset with gcc and clang, and the classic and the NLTE nebular presets additionally with Apple
Clang, nvc++, and hipcc, including the GPU code paths. CI installs artistools and plots the output of each test, so a
change of an output format fails CI there.

Once each week, a separate workflow builds the site with the `linkcheck` builder of Sphinx and tests every link.

A change that must not alter the numerical results has to give identical checksums. For such a change, a mismatch is an
error in the change. A change that alters the results for a good reason needs new reference checksums, which the
"Update checksums" workflow writes.
