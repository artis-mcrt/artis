# Contributing to ARTIS

Thank you for your interest in ARTIS. This page gives the steps of a change. The
[documentation site](https://artis-mcrt.readthedocs.io/) describes the code, and
[Development setup](https://artis-mcrt.readthedocs.io/en/latest/development.html) gives the tools.

No one is employed specifically to maintain ARTIS. The developers improve the code when their own research needs it,
so an answer can take time.

## Before you start

- Open an [issue](https://github.com/artis-mcrt/artis/issues) for a large change. The maintainers then say whether it
  fits the code.
- Make one pull request for one work item. Do not put an unrelated fix in the same request.

## Make the change

1. Make a branch from `main`.
2. Install the pre-commit hooks with `prek install`.
3. Symlink a preset to `artisoptions.h`, e.g. `ln -s artisoptions_classic.h artisoptions.h`. Do not copy the preset,
   because the symlink follows each change of the preset file.
4. Build with `make`. Run the unit tests with `make OPTIMIZE=OFF unittests && ./unittests` when you change a function
   that `unittests.cc` calls.
5. Run one end-to-end test that your change can affect, e.g. `tests/runtest.sh nebular_1d_3dgrid`. Do not run the
   whole matrix on your own machine. [Testing](https://artis-mcrt.readthedocs.io/en/latest/testing.html) describes the
   tests.
6. Run `prek run --all-files`, and run clang-tidy on each file that you changed with `run-clang-tidy <file>.cc`.

A new compile-time option goes into every `artisoptions_*.h` preset and into `artisoptions_doc.md`, under the same
heading and at the same position in each file. `scripts/make_options_docs.py` then writes `docs/options.md` again, and
a pre-commit hook stops a commit that leaves the page out of date.

## The style of the code

- Give each function, variable, and type a name that says which physical quantity it holds. The readers of this code
  are scientists, so a name must explain itself.
- Three index spaces select a cell, and all of them are plain integers. Give each variable the name of its index
  space, which is `cellindex`, `modelgridindex`, or `nonemptymgi`.
- Write a log line with `printlnlog()`, and a fatal error with `assert_always()` or with `fatal_crash()`.
- Give the full reference of each paper that you name: the authors, the year, the journal, the volume, the page
  range, and the DOI.
- Write the documentation, the comments, and the commit messages in Simplified Technical English (ASD-STE100).

[AGENTS.md](AGENTS.md) holds these conventions in full, together with the rules for the shared state, the
reproducible results, and the portability to GPUs.

## The pull request

- Push early. Continuous integration compiles the code with each supported compiler and runs the end-to-end tests, so
  it finds more than a local build. A pull request from a fork gets no automatic run.
- Say in the commit message and in the pull request when you expect the numerical results to change, and name the
  tests that change. The checksums must stay the same for every other change.
- After an intended change of the results, start the `updatechecksums.yml` workflow. It runs from `main` and takes
  your branch as an input, so give no `--ref`:

  ```sh
  gh workflow run updatechecksums.yml -f branch=<branch>
  ```

  The workflow needs write access to the repository, and the branch must be in `artis-mcrt/artis`. Ask a maintainer
  when you have neither.
- Add an entry to the first section of `docs/releasenotes.md` for a change that a user sees. Such a change alters the
  numerical results, a file format, an option, a preset, a program, or a script.
