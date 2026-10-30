# ARTIS

[![DOI](https://zenodo.org/badge/11279591.svg)](https://zenodo.org/badge/latestdoi/11279591)
[![CI](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml/badge.svg)](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml)
[![Documentation](https://readthedocs.org/projects/artis/badge/?version=latest)](https://artis.readthedocs.io/)
[![License](https://img.shields.io/github/license/artis-mcrt/artis)](https://github.com/artis-mcrt/artis/blob/main/LICENSE)

<img src="data/artislogo.png" alt="ARTIS logo" width="280" />

ARTIS is a 3D radiative transfer code that uses Monte Carlo methods with indivisible energy packets ([Lucy 2002, A&A, 384, 725-735, doi:10.1051/0004-6361:20011756](https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract)) for ejecta in homologous (ballistic) expansion such as supernovae and kilonovae. The code is designed for high performance on modern HPC clusters, with a focus on physics fidelity and multi-dimensional geometry.

## Documentation

The full documentation is at [artis.readthedocs.io](https://artis.readthedocs.io/). It comes from the [docs](docs) folder of this repository.

- [Quickstart](https://artis.readthedocs.io/en/latest/quickstart.html): a first simulation and a plot on a laptop.
- [Installation](https://artis.readthedocs.io/en/latest/installation.html): the compiler, the MPI library, and the build.
- [Option presets](https://artis.readthedocs.io/en/latest/presets.html) and [Compile-time options](https://artis.readthedocs.io/en/latest/options.html): the physics of the run.
- [Atomic and model data](https://artis.readthedocs.io/en/latest/data.html): the data that a simulation reads.
- [Runs on a cluster](https://artis.readthedocs.io/en/latest/clusters.html): the Slurm job scripts and the wall-time mechanism.
- [Development setup](https://artis.readthedocs.io/en/latest/development.html) and [Testing](https://artis.readthedocs.io/en/latest/testing.html): the tools, the linters, and the tests.

## Quickstart

```sh
git clone https://github.com/artis-mcrt/artis.git
cd artis
tests/runtest.sh kilonova_1d
```

This builds the code and runs a small 1D kilonova model in two jobs, and then post-processes the packets with `exspec`. The [Quickstart](https://artis.readthedocs.io/en/latest/quickstart.html) page describes the steps and makes a model folder of your own. The build needs gcc 14 or a later version, clang, nvc++, or hipcc, and an MPI library that supplies `mpicxx`.

## Key features

### Multi-dimensional geometry
ARTIS simulates ejecta in 1D spherical, 2D cylindrical, and 3D Cartesian coordinates, enabling the study of asymmetric explosions and geometry-dependent observables such as polarisation that are inaccessible to 1D codes.

### Physics fidelity
- **Line-by-line opacities**: In the default mode, ARTIS uses individual line opacities (Sobolev), with optional support for binned opacities.
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
