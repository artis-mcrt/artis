# Development setup

Clone the source code repository and check out the default branch `main`.

```sh
git clone https://github.com/artis-mcrt/artis.git
cd artis
```

[Installation](installation.md) gives the compiler and the MPI library that the build needs.

Then install the pre-commit hooks, select a preset, and make a compilation database for the clang tools:

```sh
prek install
ln -s artisoptions_nltenebular.h artisoptions.h
make TESTMODE=ON compile_commands.json
```

Every build writes the database again, so it stays current. It always uses the `TESTMODE=ON` flags, so that the clang
tools check the assertions of the test mode. For editing, the clangd language server is recommended, e.g. with the
[VS Code plugin](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd).

## The three programs

The Makefile builds three programs from the same sources. `sn3d.cc`, `exspec.cc`, and `unittests.cc` hold the three
`main()` functions, and the Makefile compiles every other `.cc` file into all three programs. That list is the variable
`common_files` in the Makefile. Add a new source file to it, because the Makefile finds no file automatically.

## Linting and formatting

- clang-format uses a Google-based style with a limit of 120 columns (`.clang-format`). It applies to all C++ files
  except `third_party/`.
- clang-tidy has a long list of checks (`.clang-tidy`) and makes almost every diagnostic an error. Make the compilation
  database with the Makefile, then name the files that you changed:

    ```sh
    make TESTMODE=ON compile_commands.json
    run-clang-tidy grid.cc
    ```

    Never make the database with `compiledb`, because it drops `mpicxx` from each command and every diagnostic is then
    false.

- cppcheck runs in CI and stops the job on an error. Read the command in `.github/workflows/ci-checks.yml`.
- The pre-commit hooks (`.pre-commit-config.yaml`, run with `prek`) apply clang-format, more than a dozen whitespace
  and file checks, and a `make OPTIMIZE=OFF` compile. `prek run --all-files` runs all of them.
- The Python scripts in `scripts/` get ruff, pyrefly, and ty. `scripts/ruff.toml`, `scripts/pyrefly.toml`, and
  `scripts/ty.toml` hold the settings. The scripts use only the standard library.

## Documentation

This site comes from the `docs/` folder. The pages are [MyST](https://myst-parser.readthedocs.io/en/latest/)
Markdown, and [Read the Docs](https://about.readthedocs.com/) builds them with
[Sphinx](https://www.sphinx-doc.org/en/master/) on each push.
`docs/conf.py` and `.readthedocs.yaml` hold that configuration. To build the site locally:

```sh
uv run --with-requirements docs/requirements.txt sphinx-build -W --keep-going docs docs/_build/html
```

`-W` makes a warning an error, as CI and Read the Docs do. Open `docs/_build/html/index.html` to read the result.

Read the Docs builds one version of the site for each branch and for each tag that you activate. `latest` comes from
the `main` branch, and a version of a release comes from its branch or from its tag. The version selector under the
name of the project goes from one version to the other. `docs/_static/version-selector.js` fills that list with the
data of Read the Docs, so a local build shows no selector. A new version needs two steps. The branch or the tag must
hold the `docs/` folder and `.readthedocs.yaml`, and a maintainer must then activate the version on the dashboard.

`scripts/make_options_docs.py` writes `docs/options.md` from `artisoptions_doc.md` and the `artisoptions_*.h`
presets, and a pre-commit hook runs the script. Change the comment of the option in `artisoptions_doc.md`, not the
generated page.

## A new compile-time option

A new option must go into every `artisoptions_*.h` preset and into `artisoptions_doc.md`, under the same heading and at
the same position in each file. CI finds the presets with a glob, so a missing preset stops the build.
`scripts/make_options_docs.py` stops when a preset holds an option that `artisoptions_doc.md` does not describe.

## More for contributors

[AGENTS.md](https://github.com/artis-mcrt/artis/blob/main/AGENTS.md) holds the conventions of the code: the names, the
comments, the three cell index spaces, the logs and the assertions, the shared state, the reproducible results, and the
portability to GPUs. It also describes the workflow of a pull request.
