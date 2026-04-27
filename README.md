# ARTIS

[![DOI](https://zenodo.org/badge/11279591.svg)](https://zenodo.org/badge/latestdoi/11279591)
[![CI](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml/badge.svg)](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/artis-mcrt/artis)](https://github.com/artis-mcrt/artis/blob/develop/LICENSE)

ARTIS is a 3D radiative transfer code that uses Monte Carlo methods with indivisible energy packets ([Lucy 2002](https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract)) for ejecta in homologous (ballistic) expansion.

## Key Features

### Multi-dimensional geometry
ARTIS simulates ejecta in 1D spherical, 2D cylindrical, and 3D Cartesian coordinates, enabling the study of asymmetric explosions and geometry-dependent observables such as polarisation that are inaccessible to 1D codes.

### Physics fidelity
- **Line-by-line opacities**: In the default mode, ARTIS uses individual line opacities (Sobolev), with optional support for binned expansion opacities.
- **Macroatom radiative transfer**: the Lucy macroatom scheme self-consistently propagates packets through absorption, fluorescence, and multi-level de-excitation, capturing line-to-line energy redistribution without simplification. This can also be disabled in favour of simplified thermalisation/scattering ratio treatment.
- **Non-LTE level populations**: statistical-equilibrium non-LTE populations are solved alongside a multibin radiation field model and trajectory-based photoionisation estimators for accurate ionisation and excitation.
- **Non-thermal physics**: a detailed Spencer-Fano solver tracks the thermalization of fast electrons from radioactive decays and their contributions to ionisation and excitation rates.
- **Nuclear decay network**: alpha, beta, and fission decays are handled natively, including time-dependent Monte Carlo particle thermalisation, making ARTIS well suited for kilonova and r-process transient modelling.
- **Polarisation**: full Stokes-parameter polarised radiative transfer via both real and virtual packets enables direct comparison with spectropolarimetric observations of asymmetric ejecta.
- **Full-phase coverage**: a single simulation framework spans both the photospheric and nebular phases of a transient, eliminating the need to switch between specialised codes.

### High-performance computing
- **Modern C++23**: the codebase uses current language standards, including `constexpr`, concepts, ranges, and structured bindings, resulting in expressive and maintainable code that benefits from the full optimisation pipeline of modern compilers.
- **Distributed memory parallelism with MPI**: the simulation scales to thousands of CPU cores across multiple nodes. Intra-node communication uses MPI shared-memory windows to avoid redundant data copies between ranks on the same host.
- **Cache-friendly data layout**: data structures are organised for spatial locality so that packet updates and cell lookups achieve high cache hit rates, reducing memory-bandwidth bottlenecks on modern CPU architectures.
- **Cell-batched packet updates**: packets are processed in cell-ordered batches — analogous to ray-coherence methods used in production rendering engines — further improving instruction and data cache reuse.

### Software engineering practices
- **Continuous integration**: every pull request is validated automatically via [GitHub Actions CI](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml), catching regressions before they reach the main branch.
- **Static analysis and linting**: the project is configured for clang-tidy, clang-format, and cppcheck, enforcing consistent style and catching common C++ pitfalls at compile time rather than at runtime.
- **Address and undefined-behaviour sanitizers**: the `TESTMODE=ON` build flag enables ASan and UBSan, providing runtime detection of memory errors and undefined behaviour during development.

## Citing ARTIS
We maintain a list of [papers that use ARTIS](https://ui.adsabs.harvard.edu/user/libraries/g5NyA9gKT5KdDFLY6SixWg) and [papers that use ARTIS with non-LTE enabled](https://ui.adsabs.harvard.edu/user/libraries/CX8fnPInSu2q1rAE4wWQrg).

If you use ARTIS, please cite it using the [DOI from Zenodo](https://zenodo.org/records/18670358).

An early version of the code is described in [Sim (2007)](https://ui.adsabs.harvard.edu/abs/2007MNRAS.375..154S/abstract) and [Kromer & Sim (2009)](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1809K/abstract). Some specific features are described in:
- Polarisation and virtual packets: [Bulla et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450..967B/abstract)
- Non-LTE level populations, multibin radiation field model, trajectory-based photoionisation estimators, and non-thermal ionisation/excitation: [Shingles et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.2029S/abstract)
- Alpha, beta, and fission decay, and time-dependent particle thermalisation for kilonovae: [Shingles et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...954L..41S/abstract)
- Expansion opacities and parameterised scattering/thermalisation ratio (instead of default line-by-line opacity and macroatom): Shingles et al. (in prep)

## Source code availability and license
The ARTIS source code is available under a [BSD 3-Clause license](https://github.com/artis-mcrt/artis/blob/develop/LICENSE), which requires attribution and preservation of copyright notices on any substantial copies. If you find the ARTIS code useful in any way, we request that you cite us as described above and star the repository to help show impact in funding proposals.

## Setting up production runs on Linux
We recommended retaining the exact source code and Git metadata within each simulation folder for future reference (i.e., don't just copy the executables).

Clone the source code repository from the release branch:
```sh
git clone --branch release https://github.com/artis-mcrt/artis.git
cd artis
```

To compile and run artis, you will need a recent C++ compiler (g++, Clang, or nvc++) and an MPI library (e.g., Open MPI) that provides an `mpicxx` command. Usually, these are available on HPC clusters using module or spack commands. For systems that we use, look at the top of the relevant SLURM script in scripts/artis-*.sh to find compatible modules specifications. For Open MPI, set the C++ compiler using `export OMPI_CXX=g++`.

Next, select an options preset. For example:
```sh
ln -s artisoptions_classic.h artisoptions.h
```

You will likely want to change the number of packets per rank (MPKTS) using a text editor, e.g. `vim artisoptions.h`. The options are explained in [artisoptions_doc.md](https://github.com/artis-mcrt/artis/blob/develop/artisoptions_doc.md).

Next, compile with `make` and go up a level to the model folder:
```sh
make
cd ..
```

The next steps are to ensure a full set of snapshot files (model.txt and abundances.txt) and an atomic database are present, and to configure the timesteps in input.txt. Then, queue the relevant job script with a command such as:
```sh
sbatch artis/scripts/artis-juwels.sh
```

## Setting up for development
Clone the source code repository and checkout the default branch:
```sh
git clone https://github.com/artis-mcrt/artis.git
cd artis
```
On macOS, it's recommended that you install homebrew llvm to get clang-format, clang-tidy, clangd language server, and the clang C++ compiler. You can install this and other dependencies with:
```sh
brew install llvm open-mpi gsl prek compiledb
```
Install the pre-commit hooks and generate a compilation database for clang tools:
```sh
prek install
make clean && compiledb -n make TESTMODE=ON
```
For editing, the clangd language server is recommended (e.g., with the [VS Code plugin](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd)).

### Running
sn3d will not write to the standard output (unless a crash occurs) but each MPI rank n will produce a log file called output_n-0.txt. A local run might look something like this:
```bash
mpirun -np 8 ./sn3d&
tail -f output_0-0.txt
```
Press Ctrl+C to stop following the log file.

## Bundled Scripts
- clean.sh: Remove all output files while keeping input files and resetting the simulation to the beginning.
- movefiles.sh [DIRNAME]: Move artis output files from the simulation folder into another folder. Usually called automatically by the job scripts, but should be run manually if the simulation crashes or is terminated early.
- sumcorehourslogs.py: Calculate the summed core hours of all jobs using the timing information in the last line of the output*.txt log files. This cannot include runs where the job was terminated early.
- sumcorehoursslurm.py: Calculate the summed core hours of all jobs from the slurm job output files.

## make options
- TESTMODE=ON: Enable additional assertions and the address and undefined behaviour sanitizers.
- FASTMATH=OFF: Don't enable compiler transformations that affect round-off-level results (e.g. a*(b\*c) = (a\*b)*c).
- EIGEN=OFF: Use GSL (not Eigen) for matrix-vector solving (Spencer-Fano and NLTE pops).
- MAX_NODE_SIZE=N: Artificially limit MPI node size to N ranks. Useful for testing and preventing MPI shared memory windows from crossing CPU sockets.
- REPRODUCIBLE=ON: Use stable sorts and disable FASTMATH.
- GPU=ON: Required to compile for GPUs. Works around incompatible function calls and uses a Simpson-rule integrator in place of Gauss-Kronrod.

## Input files
### input.txt
Run-time configuration with:
- number of timesteps
- the start and end time in days
- number of pure LTE timesteps,
- optically-thick condition that switches cells to a grey opacity treatment.

### model.txt
Grid parameters, cell densities and nuclear composition.

### abundances.txt
Per-cell elemental mass fractions (in case the set of isotopic abundances in model.txt is not complete).

### adata.txt
One block per ion consisting of:

Header:<br/>
`[atomic number] [ionisation stage] [level count] [ionisation energy in eV]`

Level count rows of:<br/>
`[level index] [energy in eV] [statistical weight] [level transition count] [level name or description]`

### transitiondata.txt
One block per ion consisting of:

Header:<br/>
`[atomic number] [ionisation stage] [transition count]`

Transition count rows of:<br/>
`[lower level index] [upper level index] [A value] [collision strength] [forbidden flag]`

In case the collision strength is not available, it will be -1.0 for permitted and -2.0 for forbidden transitions.

### compositiondata.txt
Sets a filter on the elements, ion stages, and energy levels that are read in from atomic data files.

### phixsdata.txt or phixsdata_v2.txt
Photoionisation cross sections in v1 format (arbitrary energy tables) or v2 format (regularly spaced energy grid).

### gamma_*.txt
Gamma decay spectra containing energies [MeV] and probabilities.
