# ARTIS

ARTIS is a 3D Monte Carlo radiative transfer code for ejecta in homologous expansion, e.g. supernovae and kilonovae.
It follows indivisible energy packets through the ejecta, and it gives the light curves and the spectra.
[Lucy (2002), A&A, 384, 725-735, doi:10.1051/0004-6361:20011756](https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract)
describes the packets, and [How a simulation works](simulation.md) describes the method of ARTIS.

New users start with the [Quickstart](quickstart.md). It gives a first simulation and a plot on a laptop.

## Key features

### Geometry

ARTIS puts the ejecta on a 1D spherical, a 2D cylindrical, or a 3D Cartesian grid
([GRID_TYPE_OVERRIDE](options.md#grid_type_override)), so a model can be asymmetric. The output holds the light curve
and the spectrum of each observer direction (see [Output files](outputfiles.md)).

### Physics

- **Line-by-line opacity.** The default mode gives each line its own Sobolev optical depth. The binned opacities are
  an option ([RPKT_USE_BINNED_OPACITIES](options.md#rpkt_use_binned_opacities)).
- **Macro atoms.** After an absorption, a macro atom makes internal radiative and collisional transitions until it
  emits radiation or becomes a packet of thermal energy. A parameterised ratio of the scattering to the
  thermalisation replaces that treatment when you select it
  ([RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY](options.md#rpkt_boundbound_thermalisation_probability)).
- **Non-LTE populations.** The code solves the statistical equilibrium of the selected levels
  ([ION_NLEVELS_EXCITED_NLTE](options.md#ion_nlevels_excited_nlte)), together with a multibin radiation field model
  and detailed bound-free estimators
  ([DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP](options.md#detailed_bf_estimators_first_timestep)).
- **Non-thermal particles.** A Spencer-Fano solver ([Spencer & Fano (1954), Phys. Rev., 93, 1172-1181,
  doi:10.1103/PhysRev.93.1172](https://ui.adsabs.harvard.edu/abs/1954PhRv...93.1172S/abstract)) gives the heating, the
  ionisation, and the excitation by the fast leptons of the decays ([NT_SCHEME](options.md#nt_scheme)).
- **Nuclear decays.** The decay network holds the alpha, the beta, and the fission decays, and the particles
  thermalise over time ([PARTICLE_THERMALISATION_SCHEME](options.md#particle_thermalisation_scheme)).
- **Polarisation.** The code follows the Stokes parameters of the real packets and of the virtual packets
  ([POL_ON](options.md#pol_on), [VPKT_ON](options.md#vpkt_on)).
- **Both phases.** The presets cover the photospheric phase and the nebular phase of a transient (see
  [Option presets](presets.md)).

### Parallel runs

- **MPI.** A run can use many hundreds of ranks over many nodes. The ranks of one node share the large read-mostly
  arrays through MPI shared-memory windows, so each node holds one copy of the atomic data.
- **Threads and GPUs.** `OPENMP=ON` adds OpenMP threads inside each rank with gcc, clang, or nvc++, and `STDPAR=ON`
  gives the loops to the parallel algorithms of the standard library. `STDPAR=ON GPU=ON` builds for a GPU with nvc++
  or with hipcc (see [Make options](makeoptions.md)).
- **The order of the packets.** Before each timestep, the code sorts the packets by cell. The packets of one cell then
  run together and read the same cache entries.
- **Modern C++.** gcc, clang, and nvc++ compile the code as C++26, and hipcc compiles it as C++23.

### Development

- **Continuous integration.** Each push starts the
  [GitHub Actions CI](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml), except on a branch whose name
  starts with `classic`. CI compiles the code with each supported compiler and compares the checksums of the
  end-to-end tests, so a regression does not reach the main branch. A pull request from a fork gets no automatic run.
- **Static analysis.** clang-tidy, clang-format, and cppcheck run on each push, and almost every diagnostic is an
  error.
- **Sanitizers.** `TESTMODE=ON` adds the address sanitizer, the undefined behaviour sanitizer, and the extra
  assertions.

## The three programs

The Makefile builds three programs from the same sources. `make` builds the first two:

- `sn3d` does the radiative transfer over a sequence of timesteps. It writes the spectra, the light curves, and the
  packet files.
- `exspec` combines the packet files into spectra and light curves again.
  [Post-processing with exspec](exspec.md) describes when you need it.
- `unittests` runs the unit tests. `make unittests` builds it (see [Testing](testing.md)).

[artistools](https://github.com/artis-mcrt/artistools) is a companion Python package that plots and analyses the
output files.

## Citing ARTIS

ADS holds a list of [papers that use ARTIS](https://ui.adsabs.harvard.edu/user/libraries/g5NyA9gKT5KdDFLY6SixWg) and a
list of [papers that use ARTIS with non-LTE enabled](https://ui.adsabs.harvard.edu/user/libraries/CX8fnPInSu2q1rAE4wWQrg).

If you use ARTIS, please cite it using the [DOI from Zenodo](https://zenodo.org/records/18670358).

These papers describe an early version of the code: [Sim (2007), MNRAS, 375, 154-162,
doi:10.1111/j.1365-2966.2006.11271.x](https://ui.adsabs.harvard.edu/abs/2007MNRAS.375..154S/abstract) and
[Kromer & Sim (2009), MNRAS, 398, 1809-1826,
doi:10.1111/j.1365-2966.2009.15256.x](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1809K/abstract). These papers describe
specific features:

- Polarisation and virtual packets: [Bulla, Sim & Kromer (2015), MNRAS, 450, 967-981,
  doi:10.1093/mnras/stv657](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450..967B/abstract)
- Non-LTE level populations, multibin radiation field model, trajectory-based photoionisation estimators, and
  non-thermal ionisation/excitation: [Shingles et al. (2020), MNRAS, 492, 2029-2043,
  doi:10.1093/mnras/stz3412](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.2029S/abstract)
- Alpha, beta, and fission decay, and time-dependent particle thermalisation for kilonovae:
  [Shingles et al. (2023), ApJL, 954, L41,
  doi:10.3847/2041-8213/acf29a](https://ui.adsabs.harvard.edu/abs/2023ApJ...954L..41S/abstract)
- Expansion opacities and parameterised scattering/thermalisation ratio (instead of default line-by-line opacity and
  macroatom): Shingles et al. (in prep)

## Source code availability and license

The ARTIS source code is available under a
[BSD 3-Clause license](https://github.com/artis-mcrt/artis/blob/main/LICENSE), which requires attribution and
preservation of copyright notices on any substantial copies. If you find the ARTIS code useful in any way, we request
that you cite us as described above and star the repository to help show impact in funding proposals.

## Support

No one is employed specifically to maintain ARTIS. The developers improve the code when their own research needs it.
Thus ARTIS comes as it is, and the developers do not promise to correct a defect or to answer a question.

They accept a [pull request](https://github.com/artis-mcrt/artis/pulls) sometimes. You can open an
[issue](https://github.com/artis-mcrt/artis/issues), but it can stay without an answer. The log file of rank 0 holds
the git version and the compile-time options of the run. Give both in an issue.

:::{toctree}
:caption: Getting started
:maxdepth: 1
:hidden:

quickstart
installation
presets
data
example
:::

:::{toctree}
:caption: Running a simulation
:maxdepth: 1
:hidden:

simulation
methods
planning
running
clusters
outputfiles
exspec
troubleshooting
:::

:::{toctree}
:caption: Reference
:maxdepth: 1
:hidden:

options
inputfiles
makeoptions
scripts
glossary
:::

:::{toctree}
:caption: Development
:maxdepth: 1
:hidden:

development
testing
:::
