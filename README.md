# ARTIS

[![DOI](https://zenodo.org/badge/11279591.svg)](https://zenodo.org/badge/latestdoi/11279591)
[![CI](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml/badge.svg)](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/artis-mcrt/artis)](https://github.com/artis-mcrt/artis/blob/main/LICENSE)

<img src="data/artislogo.png" alt="ARTIS logo" width="280" />

ARTIS is a 3D radiative transfer code that uses Monte Carlo methods with indivisible energy packets ([Lucy 2002, A&A, 384, 725-735, doi:10.1051/0004-6361:20011756](https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract)) for ejecta in homologous (ballistic) expansion such as supernovae and kilonovae. The code is designed for high performance on modern HPC clusters, with a focus on physics fidelity and multi-dimensional geometry.

## Key features

### Multi-dimensional geometry
ARTIS simulates ejecta in 1D spherical, 2D cylindrical, and 3D Cartesian coordinates, enabling the study of asymmetric explosions and geometry-dependent observables such as polarisation that are inaccessible to 1D codes.

### Physics fidelity
- **Line-by-line opacities**: In the default mode, ARTIS uses individual line opacities (Sobolev), with optional support for binned expansion opacities.
- **Macroatom radiative transfer**: the Lucy macroatom scheme self-consistently propagates packets through absorption, fluorescence, and multi-level de-excitation, capturing line-to-line energy redistribution without simplification. This can also be disabled in favour of simplified thermalisation/scattering ratio treatment.
- **Non-LTE level populations**: statistical-equilibrium non-LTE populations are solved alongside a multibin radiation field model and trajectory-based photoionisation estimators for accurate ionisation and excitation.
- **Non-thermal physics**: a detailed Spencer-Fano solver ([Spencer & Fano 1954, Phys. Rev., 93, 1172-1181, doi:10.1103/PhysRev.93.1172](https://ui.adsabs.harvard.edu/abs/1954PhRv...93.1172S/abstract)), following [Kozma & Fransson (1992), ApJ, 390, 602-621, doi:10.1086/171311](https://ui.adsabs.harvard.edu/abs/1992ApJ...390..602K/abstract), tracks the thermalization of fast electrons from radioactive decays and their contributions to ionisation and excitation rates.
- **Nuclear decay network**: alpha, beta, and fission decays are handled natively, including time-dependent Monte Carlo particle thermalisation, making ARTIS well suited for kilonova and r-process transient modelling.
- **Polarisation**: full Stokes-parameter polarised radiative transfer via both real and virtual packets enables direct comparison with spectropolarimetric observations of asymmetric ejecta.
- **Full-phase coverage**: a single simulation framework spans both the photospheric and nebular phases of a transient, eliminating the need to switch between specialised codes.

### High-performance computing
- **Modern C++**: the codebase uses current language standards, including `constexpr`, concepts, ranges, and structured bindings, resulting in expressive and maintainable code that benefits from the full optimisation pipeline of modern compilers. Builds use C++26 with gcc, clang, and nvc++, and C++23 with the hipcc GPU compiler.
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

An early version of the code is described in [Sim (2007), MNRAS, 375, 154-162, doi:10.1111/j.1365-2966.2006.11271.x](https://ui.adsabs.harvard.edu/abs/2007MNRAS.375..154S/abstract) and [Kromer & Sim (2009), MNRAS, 398, 1809-1826, doi:10.1111/j.1365-2966.2009.15256.x](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1809K/abstract). Some specific features are described in:
- Polarisation and virtual packets: [Bulla, Sim & Kromer (2015), MNRAS, 450, 967-981, doi:10.1093/mnras/stv657](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450..967B/abstract)
- Non-LTE level populations, multibin radiation field model, trajectory-based photoionisation estimators, and non-thermal ionisation/excitation: [Shingles et al. (2020), MNRAS, 492, 2029-2043, doi:10.1093/mnras/stz3412](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.2029S/abstract)
- Alpha, beta, and fission decay, and time-dependent particle thermalisation for kilonovae: [Shingles et al. (2023), ApJL, 954, L41, doi:10.3847/2041-8213/acf29a](https://ui.adsabs.harvard.edu/abs/2023ApJ...954L..41S/abstract)
- Expansion opacities and parameterised scattering/thermalisation ratio (instead of default line-by-line opacity and macroatom): Shingles et al. (in prep)

## Source code availability and license
The ARTIS source code is available under a [BSD 3-Clause license](https://github.com/artis-mcrt/artis/blob/main/LICENSE), which requires attribution and preservation of copyright notices on any substantial copies. If you find the ARTIS code useful in any way, we request that you cite us as described above and star the repository to help show impact in funding proposals.

## Setting up production runs on Linux
We recommend retaining the exact source code and Git metadata within each simulation folder for future reference (i.e., don't just copy the executables).

Clone the source code repository from the release branch:
```sh
git clone --branch release https://github.com/artis-mcrt/artis.git
cd artis
```

To compile and run ARTIS, you will need a recent C++ compiler (gcc 14 or newer, Clang, nvc++, or hipcc) and an MPI library (e.g., Open MPI) that provides an `mpicxx` command. Usually, these are available on HPC clusters using module or spack commands. For systems that we use, look at the top of the relevant SLURM script in scripts/artis-*.sh to find compatible modules specifications. For Open MPI, set the C++ compiler using `export OMPI_CXX=g++`.

Next, select an options preset. For example:
```sh
ln -s artisoptions_classic.h artisoptions.h
```

You will likely want to change the number of packets of all ranks together (NUM_PACKETS). Use a text editor, e.g. `vim artisoptions.h`. The values in the presets are for production runs with approximately 1000 ranks. Each rank keeps its share of the packets in memory. Decrease NUM_PACKETS for a run with fewer ranks. The options are explained in [artisoptions_doc.md](https://github.com/artis-mcrt/artis/blob/main/artisoptions_doc.md).

Next, compile with `make` and go up a level to the model folder:
```sh
make
cd ..
```

The next steps are to ensure a full set of input model files (model.txt and abundances.txt) and an atomic database are present, and to configure the timesteps in input.txt. Then, edit the relevant job script (setting e.g. number of ranks and your email address) and queue it with `sbatch`. For example, on the Cosma8 cluster:
```sh
vim artis/scripts/artis-cosma8.sh
sbatch artis/scripts/artis-cosma8.sh
```

The job scripts resubmit themselves until the simulation is finished and then queue a post-processing job (scripts/exspec-zip-*.sh) that merges and compresses the output files. sn3d writes the spectra and light curves itself, so this job runs exspec only if emission.out is absent. The post-processing job also converts the packets and estimator files to parquet format with [artistools](https://github.com/artis-mcrt/artistools), which it installs automatically using uv. See [Post-processing with exspec](#post-processing-with-exspec) below for what exspec does.

## Setting up for development
> [!IMPORTANT]
> **The develop branch is now called main**: Run the following commands to update your git repo:
> ```sh
> git branch -m develop main
> git fetch origin
> git branch -u origin/main main
> git remote set-head origin -a
> ```

Clone the source code repository and check out the default branch `main` (`release` is the production branch):
```sh
git clone https://github.com/artis-mcrt/artis.git
cd artis
```
On macOS, it's recommended that you install homebrew llvm to get clang-format, clang-tidy, clangd language server, and the clang C++ compiler. You can install this and other dependencies with:
```sh
brew install llvm open-mpi prek
```
Install the pre-commit hooks and generate a compilation database for clang tools:
```sh
prek install
make TESTMODE=ON compile_commands.json
```
Every build writes the database again, so it stays current. It always uses the `TESTMODE=ON` flags, so that the clang tools check the assertions of the test mode.
For editing, the clangd language server is recommended (e.g., with the [VS Code plugin](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd)).

### Running
sn3d writes one line with the name of the job folder to the standard output, and then nothing more unless a crash occurs. Each MPI rank n writes a log file called output_n-0.txt into the job folder (see [Output files](#output-files)). sn3d keeps a symlink output_0-0.txt in the simulation folder that points to the log of rank 0. A local run might look something like this:
```bash
mpirun -np 8 ./sn3d&
tail -f output_0-0.txt
```
Press Ctrl+C to stop following the log file.

To split a long simulation across several queued jobs, run sn3d with `-w WALLTIMELIMITHOURS`. When too little wall time remains to complete another timestep, the run finishes cleanly (writing the restart files and updating input.txt) and prints RESTART_NEEDED into the log, which the bundled cluster job scripts detect to submit a continuation job. The scripts pass the remaining SLURM allocation time automatically. Run `./sn3d -h` to list all command-line options.

### Output files
Each job writes the following into its job folder, e.g. `job_from_ts0000`:
- output_n-0.txt: a log file for each MPI rank n.
- estimators_nnnn.out: the plasma conditions of each cell (temperatures, ionisation, heating and cooling rates) at each timestep.

A run writes the following into the simulation folder:
- packets00_nnnn.out: the Monte Carlo packets from each rank, which exspec can turn into spectra and light curves again.
- light_curve.out, spec.out, and the other spectrum files that [Post-processing with exspec](#post-processing-with-exspec) lists: sn3d writes the light curves and spectra at each timestep, and the emission, absorption, and direction-resolved files at the last requested timestep.
- deposition.out: the radioactive energy deposition rate as a function of time.
- gridsave_ts*.tmp and packets_*_ts*.tmp: restart files that allow a later job to continue from the end of a timestep.

sn3d writes the per-job output files (the rank log files and the estimators, nlte, radfield, and macroatom files) into a job folder. sn3d names the folder from the start timestep of the job, e.g. `job_from_ts0000` for a new simulation and `job_from_ts0008` for a job that resumes at timestep 8. Rank 0 writes the name of the job folder to the standard output, so the log of a Slurm job names its folder. A new simulation first removes the files of the previous simulation: the output files, the restart files, and the `job_from_ts*` and `*.slurm` folders. It removes the same files as `scripts/clean.sh`, except `slurm-*.out`, `machine.file.*`, and `core.*`, which can belong to the current job. sn3d writes the shared run-level files, including the restart files, to the simulation folder. It also keeps an `output_0-0.txt` symlink there that points to the rank-0 log of the current job.

### Post-processing with exspec
As well as sn3d, `make` builds exspec, which combines the packet files from all ranks into spectra and light curves. sn3d writes the same files itself, so exspec is necessary only to make them again from the packet files, e.g. after a change of MNUBINS or of the frequency range. Run it in the simulation folder with any number of ranks from one up to the number of packet files:
```bash
mpirun -np 8 ./exspec
```
exspec reads the same input.txt, model, and atomic data files as sn3d, so it must be run in the same folder. The nprocs_exspec line of input.txt gives the number of packet files, which is the number of sn3d ranks and not the number of exspec ranks. Each exspec rank reads a block of the packet files, and the ranks of one node share the spectra arrays. Only rank 0 writes a log, to exspec.txt.

It writes light_curve.out, spec.out, emission.out, emissiontrue.out, and absorption.out, plus gamma_light_curve.out and gamma_spec.out when KEEP_ESCAPED_GAMMAS is set, and specpol.out, emissionpol.out and absorptionpol.out when POL_ON is set. Direction-resolved versions of these go in the speclc_angle_res folder.

To plot and analyse the output, use [artistools](https://github.com/artis-mcrt/artistools), a companion Python package for working with ARTIS light curves, spectra, and estimators.

### Testing
Unit tests for the pure numeric and parsing helpers are built and run with `make unittests && ./unittests` (CI runs them for the classic and NLTE nebular presets).

The tests folder contains eleven small end-to-end test models. Each tests/setup_*.sh script downloads the atomic data it needs and assembles a folder that is ready to run:
```sh
cd tests
source ./setup_kilonova_1d.sh   # creates tests/kilonova_1d_testrun/
```
[CI](.github/workflows/ci.yml) runs all of these on every push. It builds with `REPRODUCIBLE=ON FASTMATH=OFF MAX_NODE_SIZE=2` and compares md5 checksums of the output files against reference checksums stored in the tests/*_inputfiles folders. A change that legitimately alters the numerical results therefore needs new reference checksums, which maintainers regenerate using the "Update checksums" workflow. CI also compiles every artisoptions_*.h preset with gcc and clang, and the classic and NLTE nebular presets additionally with Apple Clang, nvc++, and hipcc (including the GPU code paths).

## Bundled scripts
- clean.sh: Remove all output files while keeping input files and resetting the simulation to the beginning. The script also removes the job_from_ts* folders of sn3d and the *.slurm folders of older versions.
- sumcorehourslogs.py: Sum the core hours of all jobs from the output_0-0.txt log of each job. The script reads the summary in the last line of the log. For a job that stopped early, it estimates the core hours from the first and the last timestamp of the log.
- sumcorehoursslurm.py: Calculate the summed core hours of all jobs from the slurm job output files.

## Make options
- TESTMODE=ON: Enable additional assertions and the address and undefined behaviour sanitizers.
- FASTMATH=OFF: Don't enable compiler transformations that affect round-off-level results (e.g. a*(b\*c) = (a\*b)*c).
- MAX_NODE_SIZE=N: Artificially limit MPI node size to N ranks. Useful for testing and preventing MPI shared memory windows from crossing CPU sockets.
- REPRODUCIBLE=ON: Use stable sorts and disable FASTMATH.
- GPU=ON: Required to compile for GPUs. Works around incompatible function calls and uses a Simpson-rule integrator in place of Gauss-Kronrod.
- OPENMP=ON: Parallelise with OpenMP threads within each MPI rank. Cannot be combined with STDPAR.
- STDPAR=ON: Parallelise with C++ standard library parallel algorithms, either on multicore CPUs or, together with GPU=ON, on GPUs.

- GPUARCH=N: With nvc++ and GPU=ON, compile for compute capability N (e.g. GPUARCH=80) in place of the GPU of the host. A host that has no GPU needs this option to get a defined target.
- OPTIMIZE=OFF: Compile without optimisation. This is the quickest way to check that everything still compiles.
- ZSTD=OFF: Build without libzstd. The Makefile finds libzstd with pkg-config or with a link of -lzstd, and a build with the library reads a zstd compressed input file, e.g. model.txt.zst, when the plain file is absent. ZSTD=ON stops the build if the library is absent.
- PGO=GENERATE and PGO=USE: Profile-guided optimisation with gcc or clang. Build with GENERATE, run a representative simulation to collect profile data, then rebuild with USE.

## Input files
These files go in the simulation folder, which should always contain the ARTIS source folder (or a symlink to it) named artis. The physics data bundled with the code (nuclear decay data, gamma-ray spectra, and the collisional ionisation and binding energy tables) is then found automatically in artis/data and does not need to be copied.

A build with libzstd (see "Make options") also reads a zstd compressed text input file, e.g. model.txt.zst, transitiondata.txt.zst, or packets00_0000.out.zst for exspec. The plain file has priority when both exist. The exceptions are input.txt, which sn3d writes again at the start of a new simulation, and the restart files.

### input.txt
Run-time configuration with:
- the random number seed, which must be a fixed value for reproducible runs
- number of timesteps
- the first and last timestep of this job. Normally the last timestep should be set to the number of timesteps: long simulations are split over resubmitted jobs by the wall-time mechanism (sn3d -w), which advances the start timestep on each restart but never the last one. Setting an earlier last timestep makes the simulation stop there (e.g. to inspect or post-process partial results) until input.txt is edited to continue
- the start and end time in days
- whether the run continues from the restart files of a previous job
- number of pure LTE timesteps
- optically-thick condition that switches cells to a grey opacity treatment
- nprocs_exspec: the number of packet files that exspec will read, one for each sn3d rank. sn3d sets this line itself. At the start of a new simulation, sn3d writes input.txt again with a comment on each line and with this value, and it keeps the same content in input-newrun.txt.

The format is positional, so every line must be present and in order, including those that are no longer used. See [tests/classicmode_3d_inputfiles/input-newrun.txt](tests/classicmode_3d_inputfiles/input-newrun.txt) for a complete example with a comment on every line.

### vpkt.txt
Required when virtual packets are enabled (VPKT_ON in artisoptions.h). Sets the observer directions, the number of spectra per observer, and the time and wavelength ranges to record. See [tests/classicmode_3d_inputfiles/vpkt.txt](tests/classicmode_3d_inputfiles/vpkt.txt) for an example.

### model.txt
Grid parameters, cell densities and nuclear composition. The optional column `q` holds the trapped radiation energy per mass [erg/g] at the snapshot time of the model. The value must already include the adiabatic losses before the snapshot time. `INITIAL_PACKETS_ON` uses it for the initial temperature and for the initial packets.

### abundances.txt
Required file with the per-cell elemental mass fractions, which include elements whose isotopic abundances are not given in model.txt.

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

A legacy four-column format (`[transition index] [lower level index] [upper level index] [A value]`, with no collision strength or forbidden columns) is also accepted. The format is detected from the column count of the first row of each block.

### compositiondata.txt
Sets a filter on the elements, ion stages, and energy levels that are read in from atomic data files. The code ignores the text to the right of a `#` character.

### phixsdata.txt or phixsdata_v2.txt
Photoionisation cross sections in v1 format (arbitrary energy tables) or v2 format (regularly spaced energy grid).

### gamma_*.txt
Gamma decay spectra containing energies [MeV] and probabilities.
