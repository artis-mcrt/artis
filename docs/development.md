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

Name the files that you changed in the last command. The section "Linting and formatting" of
[AGENTS.md](https://github.com/artis-mcrt/artis/blob/main/AGENTS.md) describes each tool and its settings.

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

`scripts/make_options_docs.py` writes `docs/options.md` from `artisoptions_doc.md` and the `artisoptions_*.h`
presets, and a pre-commit hook runs the script. Change the comment of the option in `artisoptions_doc.md`, not the
generated page.

## A new compile-time option

A new option must go into every `artisoptions_*.h` preset and into `artisoptions_doc.md`, under the same heading and at
the same position in each file. CI finds the presets with a glob, so a missing preset stops the build.
`scripts/make_options_docs.py` stops when a preset holds an option that `artisoptions_doc.md` does not describe.

## More for contributors

[CONTRIBUTING.md](https://github.com/artis-mcrt/artis/blob/main/CONTRIBUTING.md) gives the steps of a change, from
the branch to the pull request.

[AGENTS.md](https://github.com/artis-mcrt/artis/blob/main/AGENTS.md) holds the conventions of the code: the names, the
comments, the three cell index spaces, the logs and the assertions, the shared state, the reproducible results, and the
portability to GPUs. It also describes the workflow of a pull request.
