# ARTIS

[![DOI](https://zenodo.org/badge/11279591.svg)](https://zenodo.org/badge/latestdoi/11279591)
[![CI](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml/badge.svg)](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml)

[ARTIS](https://github.com/artis-mcrt/artis) ([Sim 2007](https://ui.adsabs.harvard.edu/abs/2007MNRAS.375..154S/abstract); [Kromer & Sim 2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1809K/abstract)) is a 3D radiative transfer code for Type Ia supernovae using the Monte Carlo method with indivisible energy packets ([Lucy 2002](https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract)). The latest version incorporates polarisation and virtual packets ([Bulla et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450..967B/abstract)), non-LTE physics appropriate for the nebular phase of Type Ia supernovae ([Shingles et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.2029S/abstract)), and alpha- and beta-decays with time-dependent thermalisation ([Shingles et al. 2023](https://ui.adsabs.harvard.edu/abs/2023ApJ...954L..41S/abstract)).

The code is modern C++23 and scales to thousands of CPU cores across multiple node using MPI with shared memory windows on each node. Experimental support is also provided for OpenMP and C++ standard parallelism (for multicore CPU and upcoming GPU targets).

## Why is this code available?
The ARTIS source code is available because it forms part of the method used to obtain published scientific results. Those interested in understanding the numerical techniques (in greater detail than the published descriptions) now have access to this code. We anticipate that some developers might find our code useful when building similar simulation codes, and in this case we ask that authors of any derivative works acknowledge and cite the ARTIS collaboration. This is in addition to the legal requirements of attribution and preservation of copyright notices on any substantial copies under the BSD 3-Clause licence.

## Can you help me to run the code?
We do not have the resources to support users of the code outside our team of direct collaborators.

## Installation and development
Clone the source code repository:
```sh
git clone https://github.com/artis-mcrt/artis.git
cd artis
```
For macOS, it's recommend that you install homebrew llvm to get clang-format, clang-tidy, clangd language server and the clang C++ compiler.
```sh
brew install llvm
```
Install an MPI library that provides an mpicxx command, e.g., on macOS:
```sh
brew install open-mpi
```
Install the pre-commit hooks:
```sh
pip install pre-commit
pre-commit install
```
For editing, the clangd language server is recommended (e.g., with the [VS Code plugin](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd)).

## License

Distributed under the BSD 3-Clause license. See [LICENSE](https://github.com/artis-mcrt/artis/blob/nebular/LICENSE) for more information.
