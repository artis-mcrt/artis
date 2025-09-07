# ARTIS Copilot Instructions

## Repository Overview

**ARTIS** is a 3D radiative transfer code for Type Ia supernovae using the Monte Carlo method. It's written in modern C++23 and scales to thousands of CPU cores using MPI with shared memory windows. The codebase contains ~26,000 lines of C++ across ~40 source files.

**Project Type**: Scientific simulation software (computational astrophysics)  
**Languages**: C++23, Python (utilities), Shell scripts  
**Target Runtime**: Linux/Unix systems with MPI  
**Dependencies**: MPI (mpicxx), GSL (GNU Scientific Library), clang/gcc compilers

## Build and Validation Instructions

### Prerequisites
Always install these dependencies before building:
```bash
sudo apt install -y libgsl-dev openmpi-bin libopenmpi-dev
```

### Required Configuration
Before building, **always** copy one of the artisoptions*.h files to artisoptions.h:
```bash
cp artisoptions_classic.h artisoptions.h  # For classic mode
# OR choose appropriate options file for your simulation type
```

### Build Commands

**Primary build system**: Use Makefile (NOT CMake). CMake build may fail due to C++23 features. Build times: ~2-5 minutes.

**Basic build** (development/debugging):
```bash
make OPTIMIZE=OFF TESTMODE=OFF -j$(nproc) sn3d exspec
```

**If build fails with dependency errors**, create build directory first:
```bash
mkdir -p build/x86_64_optimizeoff  # or appropriate arch directory
make OPTIMIZE=OFF TESTMODE=OFF -j$(nproc) sn3d exspec
```

**Optimized build** (production):
```bash
make -j$(nproc) sn3d exspec
```

**Common build options**:
- `TESTMODE=ON/OFF` - Enables sanitizers and debug assertions
- `OPTIMIZE=OFF` - Disables optimizations for faster builds  
- `OPENMP=ON` - Enable OpenMP parallelization
- `GPU=ON` - Enable GPU support (experimental)
- `REPRODUCIBLE=ON` - Ensure reproducible results
- `STATICGSL=ON` - Link GSL statically

**Output**: Creates symlinks `sn3d` and `exspec` pointing to executables in `build/$(arch)/`

### Validation Steps

**Clean build test**:
```bash
make clean && make OPTIMIZE=OFF -j$(nproc) sn3d exspec
```

**Linting** (clang-tidy requires compilation database):
```bash
# Note: make check currently fails without compilation database
# Pre-commit hooks handle formatting automatically
```

**Pre-commit hooks** (always run before commits):
```bash
pip install pre-commit
pre-commit install
pre-commit run --all-files  # Manual run
```

### Test Infrastructure

**CI runs 12 test configurations** in `.github/workflows/ci.yml`:
- classicmode (1D/3D grids)
- kilonova (various grid types) 
- nebular (NLTE physics)
- Different physics options

**Local test setup**:
```bash
cd tests/
source setup_classicmode_1d_3dgrid.sh  # Downloads test data
# Build with: make REPRODUCIBLE=ON TESTMODE=OFF MAX_NODE_SIZE=2 FASTMATH=OFF
# Run with: mpirun -np 4 --oversubscribe ./sn3d
```

**Test timing**: Small tests ~5-10 minutes, full CI suite ~60-120 minutes

## Project Layout and Architecture

### Root Directory Structure
```
artisoptions*.h     - Configuration files for different simulation modes
sn3d.cc/sn3d.h      - Main Monte Carlo radiative transfer code  
exspec.cc/exspec.h  - Spectrum extraction utility
*.cc/*.h            - Core simulation modules (~40 files)
Makefile            - Primary build system (complex, 300+ lines)
CMakeLists.txt      - Secondary build system (simpler)
data/               - Physical constants and atomic data
tests/              - Test configurations and setup scripts
scripts/            - Utility scripts for job management
.github/workflows/  - CI configuration (ci.yml, ci-checks.yml)
```

### Key Source Files by Function
**Core simulation**: `sn3d.cc`, `grid.cc`, `packet.cc`, `radfield.cc`  
**Physics modules**: `ltepop.cc`, `nltepop.cc`, `nonthermal.cc`, `decay.cc`  
**Packet types**: `rpkt.cc` (r-packets), `kpkt.cc` (k-packets), `vpkt.cc` (virtual)  
**I/O and utilities**: `input.cc`, `spectrum_lightcurve.cc`, `stats.cc`

### Configuration System
**artisoptions.h is required for all builds**. Different predefined options:
- `artisoptions_classic.h` - Classic supernova models
- `artisoptions_nltenebular.h` - NLTE nebular phase
- `artisoptions_kilonova_lte.h` - Kilonova simulations

Key compile-time options include:
- `MPKTS` - Number of energy packets per MPI process
- `GRID_TYPE` - Spatial grid type (Cartesian/cylindrical/spherical) 
- `NLTEITER` - Max NLTE iterations
- Various physics toggles

### CI/CD Pipeline
**GitHub Actions workflows**:
- `ci.yml` - Main test suite (12 test configurations)
- `ci-checks.yml` - Quick checks
- `cislowtestmode.yml` - Extended testing with sanitizers

**Pre-commit hooks**:
- clang-format (code formatting)
- Standard checks (trailing whitespace, JSON/YAML validation)
- Basic make build test

### Common Build Issues and Workarounds

**Issue**: `mpicxx: command not found`  
**Fix**: Install MPI development packages first

**Issue**: GSL library errors during linking  
**Fix**: Install libgsl-dev, check with `pkg-config --libs gsl`

**Issue**: C++23 not supported by compiler  
**Fix**: Makefile auto-detects and falls back to C++20/C++23 as needed

**Issue**: CMake build fails with C++23 feature errors  
**Fix**: Use Makefile instead - CMake support is incomplete

**Issue**: OpenMPI cast-function-type warnings  
**Expected**: These are harmless warnings from OpenMPI headers, not code issues

**Issue**: Build fails with "opening dependency file" errors  
**Fix**: Create build directory manually: `mkdir -p build/x86_64_optimizeoff`

**Issue**: run-clang-tidy fails with "could not find compilation database"  
**Workaround**: Use pre-commit hooks for formatting instead

**Issue**: Tests failing due to missing atomic data  
**Fix**: Test setup scripts automatically download required data files

### Performance Notes
- Build uses architecture-specific optimizations (`-march=native`)
- Link-time optimization enabled by default in optimized builds
- Fast-math optimizations enabled unless `FASTMATH=OFF`
- Build artifacts cached in `build/$(arch)_$(options)/` directories

### File Locations Summary
- **Source code**: Repository root (`*.cc`, `*.h`)
- **Configuration**: `artisoptions*.h` (copy to `artisoptions.h`)  
- **Build outputs**: `build/` directory, symlinks in root
- **Test data**: Downloaded to `tests/*/` directories
- **Documentation**: `README.md`, `artisoptions_doc.md`
- **Scripts**: `scripts/` (job management, post-processing)

### Trust These Instructions
The information above has been verified by building the code and examining the CI workflows. Only search for additional information if these instructions are incomplete or if you encounter errors not covered here.