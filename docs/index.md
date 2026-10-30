# ARTIS

ARTIS is a 3D radiative transfer code that uses Monte Carlo methods with indivisible energy packets ([Lucy 2002, A&A, 384, 725-735, doi:10.1051/0004-6361:20011756](https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract)) for ejecta in homologous (ballistic) expansion such as supernovae and kilonovae. The code is designed for high performance on modern HPC clusters, with a focus on physics fidelity and multi-dimensional geometry.

New users start with the [Quickstart](quickstart.md). It gives a first simulation and a plot on a laptop.

## Key features

### Multi-dimensional geometry

ARTIS simulates ejecta in 1D spherical, 2D cylindrical, and 3D Cartesian coordinates, enabling the study of asymmetric explosions and geometry-dependent observables such as polarisation that are inaccessible to 1D codes.

### Physics fidelity

- **Line-by-line opacities**: In the default mode, ARTIS uses individual line opacities (Sobolev), with optional support for binned opacities.
- **Macroatom radiative transfer**: the Lucy macroatom scheme self-consistently propagates packets through absorption, fluorescence, and multi-level de-excitation, capturing line-to-line energy redistribution without simplification. This can also be disabled in favour of simplified thermalisation/scattering ratio treatment.
- **Non-LTE level populations**: statistical-equilibrium non-LTE populations are solved alongside a multibin radiation field model and trajectory-based photoionisation estimators for accurate ionisation and excitation.
- **Non-thermal physics**: a detailed Spencer-Fano solver ([Spencer & Fano 1954, Phys. Rev., 93, 1172-1181, doi:10.1103/PhysRev.93.1172](https://ui.adsabs.harvard.edu/abs/1954PhRv...93.1172S/abstract)), following [Kozma & Fransson (1992), ApJ, 390, 602-621, doi:10.1086/171311](https://ui.adsabs.harvard.edu/abs/1992ApJ...390..602K/abstract), tracks the thermalisation of fast electrons from radioactive decays and their contributions to ionisation and excitation rates.
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

## The three programs

`make` builds three programs from the same sources:

- `sn3d` does the radiative transfer over a sequence of timesteps. It writes the spectra, the light curves, and the packet files.
- `exspec` combines the packet files into spectra and light curves again. [Post-processing with exspec](exspec.md) describes when you need it.
- `unittests` runs the unit tests (see [Testing](testing.md)).

[artistools](https://github.com/artis-mcrt/artistools) is a companion Python package that plots and analyses the output files.

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

## Support

The authors give ARTIS as it is. They do not promise to correct a defect, and they do not promise to answer a
question. They accept a [pull request](https://github.com/artis-mcrt/artis/pulls) sometimes.

You can open an [issue](https://github.com/artis-mcrt/artis/issues), but it can stay without an answer. The log file
of rank 0 holds the git version and the compile-time options of the run. Give both in an issue.

:::{toctree}
:caption: Getting started
:maxdepth: 1
:hidden:

quickstart
installation
presets
data
:::

:::{toctree}
:caption: Running a simulation
:maxdepth: 1
:hidden:

running
clusters
outputfiles
exspec
:::

:::{toctree}
:caption: Reference
:maxdepth: 1
:hidden:

options
inputfiles
makeoptions
scripts
:::

:::{toctree}
:caption: Development
:maxdepth: 1
:hidden:

development
testing
:::
