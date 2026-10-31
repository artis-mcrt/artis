# ARTIS

[![DOI](https://zenodo.org/badge/11279591.svg)](https://zenodo.org/badge/latestdoi/11279591)
[![CI](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml/badge.svg)](https://github.com/artis-mcrt/artis/actions/workflows/ci.yml)
[![Documentation](https://readthedocs.org/projects/artis-mcrt/badge/?version=latest)](https://artis-mcrt.readthedocs.io/)
[![License](https://img.shields.io/github/license/artis-mcrt/artis)](https://github.com/artis-mcrt/artis/blob/main/LICENSE)

<img src="data/artislogo.png" alt="ARTIS logo" width="280" />

ARTIS is a 3D radiative transfer code that uses Monte Carlo methods with indivisible energy packets ([Lucy 2002, A&A, 384, 725-735, doi:10.1051/0004-6361:20011756](https://ui.adsabs.harvard.edu/abs/2002A%26A...384..725L/abstract)) for ejecta in homologous (ballistic) expansion such as supernovae and kilonovae. The code is designed for high performance on modern HPC clusters, with a focus on physics fidelity and multi-dimensional geometry.

## Documentation

The full documentation is at [artis-mcrt.readthedocs.io](https://artis-mcrt.readthedocs.io/). It comes from the [docs](docs) folder of this repository.

- [Quickstart](https://artis-mcrt.readthedocs.io/en/latest/quickstart.html): a first simulation and a plot on a laptop.
- [Installation](https://artis-mcrt.readthedocs.io/en/latest/installation.html): the compiler, the MPI library, and the build.
- [Option presets](https://artis-mcrt.readthedocs.io/en/latest/presets.html) and [Compile-time options](https://artis-mcrt.readthedocs.io/en/latest/options.html): the physics of the run.
- [Atomic and model data](https://artis-mcrt.readthedocs.io/en/latest/data.html): the data that a simulation reads.
- [Runs on a cluster](https://artis-mcrt.readthedocs.io/en/latest/clusters.html): the Slurm job scripts and the wall-time mechanism.
- [Development setup](https://artis-mcrt.readthedocs.io/en/latest/development.html) and [Testing](https://artis-mcrt.readthedocs.io/en/latest/testing.html): the tools, the linters, and the tests.

## Quickstart

```sh
git clone https://github.com/artis-mcrt/artis.git
cd artis
tests/runtest.sh kilonova_1d
```

This builds the code and runs a small 1D kilonova model in two jobs, and then post-processes the packets with `exspec`. The [Quickstart](https://artis-mcrt.readthedocs.io/en/latest/quickstart.html) page describes the steps and makes a model folder of your own. The build needs gcc 14 or a later version, clang, nvc++, or hipcc, and an MPI library that supplies `mpicxx`. The command above also needs libzstd, `zstd`, `md5sum`, `curl`, `tar`, `python3`, `perl`, and `rsync`. [Installation](https://artis-mcrt.readthedocs.io/en/latest/installation.html) gives the package of each system.

## Key features

- **Multi-dimensional geometry**: 1D spherical, 2D cylindrical, and 3D Cartesian ejecta, with polarised radiative transfer through real and virtual packets.
- **Physics fidelity**: line-by-line Sobolev opacities, the Lucy macroatom, non-LTE level populations with a multibin radiation field model, a Spencer-Fano solver for the non-thermal electrons, and alpha, beta, and fission decays.
- **Full-phase coverage**: one framework for the photospheric and the nebular phase of a transient.
- **High-performance computing**: modern C++ with MPI, shared-memory windows inside a node, optional OpenMP or C++ standard parallelism, and GPU builds with nvc++ or hipcc.

[Key features](https://artis-mcrt.readthedocs.io/en/latest/#key-features) on the documentation site describes each one, with the references.

## Citing ARTIS

If you use ARTIS, please cite it using the [DOI from Zenodo](https://zenodo.org/records/18670358). [Citing ARTIS](https://artis-mcrt.readthedocs.io/en/latest/#citing-artis) lists the paper of each feature, and we maintain a list of [papers that use ARTIS](https://ui.adsabs.harvard.edu/user/libraries/g5NyA9gKT5KdDFLY6SixWg).

An early version of the code is described in [Sim (2007), MNRAS, 375, 154-162, doi:10.1111/j.1365-2966.2006.11271.x](https://ui.adsabs.harvard.edu/abs/2007MNRAS.375..154S/abstract) and [Kromer & Sim (2009), MNRAS, 398, 1809-1826, doi:10.1111/j.1365-2966.2009.15256.x](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1809K/abstract).

## Source code availability and license

The ARTIS source code is available under a [BSD 3-Clause license](https://github.com/artis-mcrt/artis/blob/main/LICENSE), which requires attribution and preservation of copyright notices on any substantial copies. If you find the ARTIS code useful in any way, we request that you cite us as described above and star the repository to help show impact in funding proposals.

## Support

No one is employed specifically to maintain ARTIS. The developers improve the code when their own research needs it. Thus ARTIS comes as it is, and the developers do not promise to correct a defect or to answer a question. They accept a pull request sometimes.
