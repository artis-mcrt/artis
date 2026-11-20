# Development setup

Clone the source code repository and check out the default branch `main`.

```sh
git clone https://github.com/artis-mcrt/artis.git
cd artis
```

[Installation](installation.md) gives the compiler and the MPI library that the build needs. The hooks below need
`prek`, which the homebrew command of that page installs. On Linux, install it with your package manager or with
[uv](https://docs.astral.sh/uv/):

```sh
uv tool install prek
export PATH="$(uv tool dir --bin):$PATH"
```

[Quickstart](quickstart.md#4-look-at-the-results) says why the export is necessary.

The hooks also need python3, version 3.11 or a later version, because `scripts/make_options_docs.py` runs as one of
them.

Then install the pre-commit hooks, select a {term}`preset`, and make a compilation database for the clang tools:

```sh
prek install
ln -s artisoptions_nltenebular.h artisoptions.h
make TESTMODE=ON compile_commands.json
```

Every build writes the database again, so it stays current. It always uses the `TESTMODE=ON` flags, so that the clang
tools check the assertions of the test mode. Use the clangd language server to edit the code, e.g. with the
[VS Code plugin](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd).

## A new source file

The Makefile compiles every `.cc` file of the variable `common_files` into the three programs. Add a new source file
to it, because the Makefile finds no file automatically. [Code map](codemap.md) names the purpose of each file.

## Linting and formatting

Run these commands before a commit:

```sh
prek run --all-files
make TESTMODE=ON compile_commands.json
run-clang-tidy grid.cc
```

In the last command, name the `.cc` files that you changed. For a changed header, name the `.cc` files that include
it. The database holds only the `.cc` files. Continuous integration (CI) runs the same tools, and an error of a
tool makes the CI job of that tool fail:

- The pre-commit hooks apply clang-format with the settings of `.clang-format` to the C++ files outside `third_party/`.
  They also do the checks of the whitespace and of the files. They run the two documentation scripts and the checks of
  the Python scripts. CI skips the `make` hook, which compiles the code.
- clang-tidy reads `compile_commands.json` and makes almost every diagnostic an error (see `.clang-tidy`). The database
  holds the flags of the compiler that `mpicxx` calls, so only a database of clang has `-Wunsafe-buffer-usage`. In CI,
  `mpicxx` calls clang++ 23, and clang-tidy 23 checks the sources with the `nltenebular` preset and the libc++ headers
  of LLVM 23. `make check` runs clang-tidy over all sources with the tools, the libraries, and the preset of your
  computer. Its diagnostics can therefore be different from the diagnostics of CI.
- cppcheck checks all sources except `third_party/`. `.github/workflows/ci-checks.yml` holds its command and the
  message types that it ignores. To ignore one message, put the comment `// cppcheck-suppress <id>` on the line above
  the code. This works because the command has the option `--inline-suppr`.

Make the compilation database with a build or with `make TESTMODE=ON compile_commands.json`. Do not use `compiledb`.
That tool removes `mpicxx` from each command, so each diagnostic of clang-tidy is then false.

## Compiler warnings

The build adds `-Werror` for every compiler except nvc++, so almost every warning stops the build. clang also adds
`-Wunsafe-buffer-usage`, which gcc does not have. A change that compiles cleanly with gcc can therefore fail with clang.
CI compiles the code with both. Some system headers and the construction of a `std::span` from a raw pointer give this
warning. Put such an `#include` or such a construction between `#pragma clang unsafe_buffer_usage begin` and
`#pragma clang unsafe_buffer_usage end`.

## The Python scripts

The scripts of `scripts/` use only the standard library of Python. The hooks check them with ruff, pyrefly, and ty, and
`scripts/ruff.toml`, `scripts/pyrefly.toml`, and `scripts/ty.toml` hold the settings. These settings agree with
`pyproject.toml` of artistools, except some ruff rules that apply only to artistools. `.pre-commit-config.yaml` holds
the hooks. Add each new package of a script to `additional_dependencies` of the `pyrefly` hook and of the `ty` hook.
The type checkers otherwise cannot resolve the imports of the script.

## Documentation

This site comes from the `docs/` folder. The pages are [MyST](https://myst-parser.readthedocs.io/en/latest/)
Markdown, and [Read the Docs](https://about.readthedocs.com/) builds them with
[Sphinx](https://www.sphinx-doc.org/en/master/) on each push.
`docs/conf.py` and `.readthedocs.yaml` hold that configuration. To build the site locally:

```sh
uv run --python 3.14 --with-requirements docs/requirements.txt sphinx-build -W --keep-going docs docs/_build/html
```

`-W` makes a warning an error, as CI and Read the Docs do. Open `docs/_build/html/index.html` to read the result.

Read the Docs builds one version of the site for each branch and for each tag that a maintainer activates on the
dashboard. `latest` comes from the `main` branch. A local build shows no version selector.

`scripts/make_options_docs.py` writes `docs/options.md` and the comparison table of `docs/presets.md` from
`artisoptions_doc.md` and the `artisoptions_*.h` presets, and a pre-commit hook runs the script. Change the comment of
the option in `artisoptions_doc.md`, not the generated text.

## A new compile-time option

A new option must go into every `artisoptions_*.h` preset and into `artisoptions_doc.md`, under the same heading and at
the same position in each file. CI finds the presets with a glob, so a missing preset stops the build.
`scripts/make_options_docs.py` stops when a preset holds an option that `artisoptions_doc.md` does not describe.

When you rename or remove an option, or change its type or its form, search the repository for its name:

```sh
git grep -w <name>
```

Change each match outside the text that `scripts/make_options_docs.py` writes. Most matches are in `tests/setup_*.sh`,
in `.github/workflows/ci.yml`, and in the pages of `docs/`. Also search
[artistools](https://github.com/artis-mcrt/artistools) for the name. artistools reads some options from the
`artisoptions.h` of a model folder.

## A new make option

A new make option needs these changes:

- a suffix of `BUILD_DIR` in the Makefile, if the option changes `CXXFLAGS`;
- an entry in the `foreach` list of the Makefile, if the values of the option are `ON` and `OFF`;
- a fixed value in the `make` command of `tests/runtest.sh`;
- a line in the help text of the Makefile, and a row in the table of [Make options](makeoptions.md).

Without the suffix, two values of the option share one build folder, and make can then link objects that it compiled
with the other value. The `foreach` list stops make for a bad value, e.g. `on`. The fixed value keeps a value of your
environment out of the test build. `scripts/check_docs.py` compares the help text with the table.

## More for contributors

[CONTRIBUTING.md](https://github.com/artis-mcrt/artis/blob/main/CONTRIBUTING.md) gives the steps of a change, from
the branch to the pull request.

[AGENTS.md](https://github.com/artis-mcrt/artis/blob/main/AGENTS.md) holds the conventions of the code: the names, the
comments, the three cell index spaces, the logs and the assertions, the shared state, the reproducible results, and the
portability to GPUs. It also describes the workflow of a pull request.
