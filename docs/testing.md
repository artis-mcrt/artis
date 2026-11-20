# Testing

## The unit tests

The unit tests cover the numeric and parsing helpers and some physics functions, e.g. the Compton cross-section and the
triangular solve of the Spencer-Fano matrix. Build and run them:

```sh
make unittests && ./unittests
```

Continuous integration (CI) runs them for the classic and the {term}`NLTE` nebular {term}`presets <preset>`. The
harness is hand-written, so there is no gtest and no catch.

## The end-to-end tests

The `tests` folder holds small end-to-end test models. `tests/runtest.sh` runs one of them:

```sh
tests/runtest.sh kilonova_1d
```

`tests/runtest.sh --help` shows the names of the tests and of the steps, and `tests/runtest.sh --list` shows only the
names of the tests. The test job of
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
   the checksums of the output files with `tests/kilonova_1d_inputfiles/results_md5_job0.txt`. For the test
   `classicmode_3d_formatversion2`, the step `formatcheck` comes between the two. It runs `sn3d -c` on the small 1D
   and 2D models of format version 2 in the folder `formatcheck_models` of the input files. It also runs `sn3d -c`
   on damaged copies of them. `sn3d` must read each model and must stop with the correct message for each damaged
   copy.
4. It checks the input and the {term}`restart files <restart file>` with `sn3d -c`, and it resumes the simulation from the restart files. It
   then removes the restart files and runs `exspec`.
5. It compares the checksums of the output files with `results_md5_final.txt`.

The script stops at the first step that fails, except a checksum step that finds different checksums. The other steps
then run, and at the end the script lists the files that differ and exits with an error status. To run only one step,
give the name of the step after the name of the test, e.g. `tests/runtest.sh kilonova_1d resume`. With `TESTMODE=ON`
in the environment, the build has the sanitizers and the extra assertions. The script then sets the two variables that
a run of such a build needs:

```sh
export ASAN_OPTIONS=detect_stack_use_after_return=1:detect_leaks=0
export OMPI_MCA_memory='^patcher'
```

Each test uses one preset and a reduced atomic dataset from a release of this repository.
`tests/atomic_data_release.txt` names the release, and `tests/setupfuncs.sh` holds the download command.

| Tests | Preset | Archive | Elements |
| --- | --- | --- | --- |
| `classicmode_1d_3dgrid`, `classicmode_3d*` | `classic` | `atomicdata_classic.tar.xz` | the classic dataset, with cross-sections in the v1 format |
| `kilonova_1d`, `kilonova_2d*`, `kilonova_3dspherical*` | `kilonova_lte` | `atomicdata_feconi.tar.xz` | Fe, Co, Ni |
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

## New reference checksums

Start the "Update checksums" workflow (`.github/workflows/updatechecksums.yml`) after a change that alters the results
for a good reason. Also start it after you add a test. The checksum steps of a new test compare nothing and give only a
warning, so CI passes. In all other cases, if each test of CI passes, the reference checksums are correct, and the
workflow is not necessary. The workflow needs write access to the repository, and the branch must be in
`artis-mcrt/artis`.

1. Wait until the CI run of the head commit of the branch is complete.
2. Make sure that each failed test job failed only at a checksum step, with no crash and no failed assertion.
3. Make sure that only the tests that the change must affect have different checksums.
4. Start the workflow with the branch as the input. With gh 2.87.0 or later, the first command prints the URL of the new
   run, and the URL ends with the run ID. With an older gh, `gh run list --workflow updatechecksums.yml` shows the run
   ID. The second command waits for the run, and it fails if the run fails.

   ```sh
   gh workflow run updatechecksums.yml -f branch=<branch>
   gh run watch <run ID> --exit-status
   ```

5. If the run fails, read its log and correct the cause.
6. Read the log of the run. Make sure that it selected the artifact of the head commit.
7. Make sure that the CI run of the new commit passes.
8. Pull the new commit before your next push to the branch.

The workflow gets a token with write access, so it runs only from `main`, and the first command needs no `--ref`. The
workflow takes the checksums from the CI artifact of the nearest commit in the history of the branch. If these
checksums are different from the reference checksums, it pushes a commit with them. A failed step of a test job, e.g. a
crash, makes the artifact incomplete. The workflow then takes an older commit and gives only a warning. An older
commit is correct only if the later commits do not change the results.

## A new end-to-end test

A new test needs these parts:

- a setup script `tests/setup_<name>.sh`;
- a folder `tests/<name>_inputfiles/` with the input files of the model, e.g. `input-newrun.txt` and
  `input-resume.txt`, or with only the files that are different from the folder of another test;
- the entry `<name>` in the `testname` matrix of `.github/workflows/ci.yml`;
- the atomic data archive in the release that `tests/atomic_data_release.txt` names;
- a row in the table of the tests on this page, if no row names the test.

The setup script does these steps:

1. It gets the atomic data archive with `getatomicdata`.
2. It makes the model folder `tests/<name>_testrun`.
3. In the model folder, it extracts the archive and copies the input files.
4. It makes the symlink `artis` to the repository folder. `sn3d` finds the files of `data/` through this symlink.
5. It copies a preset from `artis/` to `artisoptions.h`.
6. It changes some options with `setoption`.

`tests/setupfuncs.sh` holds `getatomicdata` and `setoption`. The setup script of a variant first copies the input
folder of the other test. Then it copies its own input folder with `rsync --ignore-times`. The input folder must hold
at least one file. Git keeps no empty folder.

`ci.yml` keeps each archive in a cache. Each cache step selects the tests with an `if` condition on their names. The
step of `atomicdata_feconi.tar.xz` takes each name that the other steps do not select. Make sure that these conditions
give the new test the cache of its archive. A new archive also needs its own cache step.

The plot steps of `ci.yml` also select a time from the start of the test name, e.g. `-t 200` for `nebular`. artistools
stops with an error if this time is outside the valid time range of the test. Make sure that the time of each plot
step is in the valid time range of the new test.

`tests/runtest.sh` finds the test from its setup script. Commit the test with no reference checksums. A checksum step
without a reference file compares nothing and gives a warning. After the first complete CI run, the "Update checksums"
workflow writes the reference checksums.

A new test must run in CI in about ten minutes or less. Every push runs the full matrix.
