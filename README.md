# ARTIS

[![DOI](https://zenodo.org/badge/11279591.svg)](https://zenodo.org/badge/latestdoi/11279591)
[![CI](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml/badge.svg)](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml)

[ARTIS](https://github.com/artis-mcrt/artis) ([Sim 2007](https://ui.adsabs.harvard.edu/abs/2007MNRAS.375..154S/abstract); [Kromer & Sim 2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1809K/abstract)) is a 3D radiative transfer code for Type Ia supernovae using the Monte Carlo method with indivisible energy packets ([Lucy 2002](https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract)). The latest version incorporates polarisation and virtual packets ([Bulla et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450..967B/abstract)), non-LTE physics appropriate for the nebular phase of Type Ia supernovae ([Shingles et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.2029S/abstract)), and alpha- and beta-decays with time-dependent thermalisation ([Shingles et al. 2023](https://ui.adsabs.harvard.edu/abs/2023ApJ...954L..41S/abstract)).

The code is modern C++23 and scales to thousands of CPU cores across multiple node using MPI with shared memory windows on each node. Experimental support is also provided for OpenMP and C++ standard parallelism (for multicore CPU and upcoming GPU targets).

## Why is this code available?
The ARTIS source code is available because it forms part of the method used to obtain published scientific results. Those interested in understanding the numerical techniques in greater detail than the published descriptions have full access to the underlying code. We anticipate that some developers might find our code useful when building similar simulation codes, and in this case we ask that authors of any derivative works acknowledge and cite the ARTIS collaboration. This is in addition to the legal requirements of attribution and preservation of copyright notices on any substantial copies under the BSD 3-Clause licence.

## Can you help me to run the code?
We do not have the resources to support users of the code outside our team of direct collaborators.

## Installation of release version for production runs on Linux
We recommended that you retain the full source code and Git version metadata within each simulation folder  for future reference (i.e. don't just copy the executables).

Clone the source code repository from the release branch:
```sh
git clone --branch release https://github.com/artis-mcrt/artis.git
cd artis
```

To compile and run artis, you require a recent C++ compiler (gcc, clang, or nvc++), the GNU Scientific Library, and an MPI library with a wrapper command `mpicxx'. Typically these are made available on an HPC system by running module or spack commands. For systems we use frequently, look at the top of the relevant SLURM script in scripts/artis-*.sh to find compatible modules to load.

With the requirements met, select an options preset and compile with `make'. For example:
```sh
ln -s artisoptions_classic.h artisoptions.h
make
```

You will likely want to change the number of packets per rank (MPKTS), and the GRID_TYPE (SPHERICAL1D, CYLINDRICAL2D, or CARTESIAN3D) using a text editor. A CARTESIAN3D grid can be used with any 1D, 2D, or 3D input model (as was historically done) but there will be a loss of accuracy due to the mismatch between model cells and propagation cells.

Next, go up to the model folder and symlink the executables and data folder:
```sh
cd ..
ln -s artis/sn3d
ln -s artis/exspec
ln -s artis/data
```

The next steps are to ensure a full set of snapshot files (model.txt and abundances.txt) and an atomic database are present, and configure the timesteps in input.txt. Then, queue the relevant job script, with a command such as:
```sh
sbatch artis/scripts/artis-juwels.sh
```

## Development
Clone the source code repository and checkout the default branch:
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
pip install prek
prek install
```
For editing, the clangd language server is recommended (e.g., with the [VS Code plugin](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd)).

## License

Distributed under the BSD 3-Clause license. See [LICENSE](https://github.com/artis-mcrt/artis/blob/nebular/LICENSE) for more information.
