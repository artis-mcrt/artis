# ARTIS

[![DOI](https://zenodo.org/badge/11279591.svg)](https://zenodo.org/badge/latestdoi/11279591)
[![CI](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml/badge.svg)](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml)

[ARTIS](https://github.com/artis-mcrt/artis) ([Sim 2007](https://ui.adsabs.harvard.edu/abs/2007MNRAS.375..154S/abstract); [Kromer & Sim 2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1809K/abstract)) is a 3D radiative transfer code for Type Ia supernovae using the Monte Carlo method with indivisible energy packets ([Lucy 2002](https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract)). The latest version incorporates polarisation and virtual packets ([Bulla et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450..967B/abstract)), non-LTE physics appropriate for the nebular phase of Type Ia supernovae ([Shingles et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.2029S/abstract)), and alpha- and beta-decays with time-dependent thermalisation ([Shingles et al. 2023](https://ui.adsabs.harvard.edu/abs/2023ApJ...954L..41S/abstract)).

The code is modern C++23 and scales to thousands of CPU cores across multiple node using MPI with shared memory windows on each node. Experimental support is also provided for OpenMP and C++ standard parallelism (for multicore CPU and upcoming GPU targets).

## Why is this code available?
The ARTIS code forms part of the method used to obtain published scientific results. Those interested in understanding the numerical techniques in greater detail than the published descriptions have full access to the underlying code. We anticipate that some developers building similar simulation codes might find our code useful, and in this case we ask that authors of any derivative works acknowledge and cite the ARTIS collaboration. This is in addition to the legal requirements of attribution and preservation of copyright notices on any substantial copies under the BSD 3-Clause licence.

## Can you help me to run the code?
We do not have the resources to support users of the code outside our team of direct collaborators.

## Setting up production runs on Linux
We recommended that you retain the full source code and Git version metadata within each simulation folder for future reference (i.e. don't just copy the executables).

Clone the source code repository from the release branch:
```sh
git clone --branch release https://github.com/artis-mcrt/artis.git
cd artis
```

To compile and run artis, you will need a recent C++ compiler (g++, LLVM Clang, Apple Clang, or nvc++), the GNU Scientific Library (or Boost and Eigen headers), and an MPI library with a wrapper command `mpicxx'. Typically these are made available on HPC systems using module or spack commands. For systems that we use, look at the top of the relevant SLURM script in scripts/artis-*.sh to find compatible modules specifications.

Next, select an options preset and compile with `make'. For example:
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
On macOS, it's recommend that you install homebrew llvm to get clang-format, clang-tidy, clangd language server, and the clang C++ compiler. You can install this and other dependencies with:
```sh
brew install llvm open-mpi gsl boost eigen prek
```
Install the pre-commit hooks:
```sh
prek install
```
For editing, the clangd language server is recommended (e.g., with the [VS Code plugin](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd)).

### Running
sn3d will not write to the standard output (unless a crash occurs) but each MPI rank n will produce a log file called output_n-0.txt. A local run might look something like this:
```bash
mpirun -np 8 ./sn3d&
tail -f output_0-0.txt
```
Press Ctrl+C to stop following the log file.

## make options
- TESTMODE=ON: Enable additional assertions and compile with address sanitizer.
- FASTMATH=OFF: Don't enable compiler transformations that affect round-off-level results (e.g. a*(b+c).
- BOOST=ON: Use Boost library (not GSL) for root finding (classic nne solver, T_e, and binned radfield) and adaptive integration.
- EIGEN=ON: Use Eigen library (not GSL) for matrix-vector solving (Spencer-Fano and NLTE pops).
- BOOST=ON EIGEN=ON: Do not use GSL at all.
- MAX_NODE_SIZE=N: Artificially limit MPI node size to N ranks. Useful for testing or preventing MPI shared memory windows from crossing CPU sockets.
- REPRODUCIBLE=ON: Use stable sorts and disable FASTMATH.
- GPU=ON: Required to compile for GPUs. Avoids incompatible GSL and std::random calls.

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

## License

Distributed under the BSD 3-Clause license. See [LICENSE](https://github.com/artis-mcrt/artis/blob/nebular/LICENSE) for more information.
