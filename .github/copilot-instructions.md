# ARTIS Radiative Transfer Code

ARTIS is a 3D radiative transfer code for Type Ia supernovae using Monte Carlo methods with indivisible energy packets. The code is written in modern C++23 and scales to thousands of CPU cores across multiple nodes using MPI with shared memory windows.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

### Prerequisites and Dependencies
Install required dependencies:
```bash
sudo apt-get update
sudo apt install -y openmpi-bin libopenmpi-dev
```

Verify installations:
```bash
mpicxx --version  # Should show Open MPI version
python3 --version  # Should show Python 3.x
clang-format --version  # Should show clang-format (for code formatting)
```

### Building the Code
**NEVER CANCEL builds - they complete in 30-60 seconds with proper timeout settings.**

**Required setup**: Copy an artisoptions configuration file:
```bash
cp artisoptions_classic.h artisoptions.h
```

**Primary build system (Makefile - RECOMMENDED):**
```bash
make -j$(nproc) sn3d exspec
# Build time: ~30 seconds. NEVER CANCEL. Set timeout to 90+ seconds.
```

**CI/Reproducible build (for testing):**
```bash
make REPRODUCIBLE=ON TESTMODE=OFF MAX_NODE_SIZE=2 FASTMATH=OFF -j$(nproc) sn3d exspec
# Build time: ~30 seconds. NEVER CANCEL. Set timeout to 90+ seconds.
```

**Verification**: After building, check executables exist:
```bash
ls -la sn3d exspec  # Should show symlinks to build/x86_64/ executables
```

### Running Tests
**NEVER CANCEL test runs - they can take 120+ minutes in CI environments.**

Setup a test case:
```bash
cd tests/
bash setup_classicmode_1d_3dgrid.sh
# Setup time: ~30 seconds (downloads atomic data if needed)
```

Run simulation:
```bash
cd classicmode_1d_3dgrid_testrun/
cp input-newrun.txt input.txt
mpirun -np 4 --oversubscribe --mca mpi_yield_when_idle 1 ./sn3d
# Runtime: Variable, can be 30+ minutes. NEVER CANCEL. Set timeout to 150+ minutes.
```

Run spectrum extraction:
```bash
mpirun -np 1 ./exspec
# Runtime: Usually < 10 minutes. Set timeout to 30+ minutes.
```

## Validation Scenarios
**CRITICAL**: Always run complete validation scenarios after making changes:

1. **Build validation**: Verify both executables build successfully
2. **Basic functionality**: Run at least one test case end-to-end
3. **Code formatting**: Run `clang-format` on modified .cc/.h files
4. **Scientific validation**: Check that simulation produces expected output files:
   - `estimators*.out` - Main estimator output
   - `*.out` files - Various output data
   - `spec*.out` - Spectral output from exspec

### Code Quality and Pre-commit
Run formatting and checks:
```bash
clang-format -i *.cc *.h  # Format modified source files
```

Pre-commit hooks (if available):
```bash
pip install pre-commit
pre-commit install
pre-commit run --all-files
```

## Key Codebase Information

### Directory Structure
```
├── *.cc, *.h           # Main C++ source code
├── artisoptions*.h     # Configuration templates for different simulation modes
├── data/               # Atomic and nuclear decay data files
├── tests/              # Test cases with input files and setup scripts
├── scripts/            # Utility scripts (mostly for HPC environments)
├── .github/workflows/  # CI configuration
├── Makefile           # Primary build system
└── README.md          # Basic project information
```

### Important Files
- **Main executables**: `sn3d.cc` (simulation), `exspec.cc` (spectrum extraction)
- **Configuration**: `artisoptions.h` (must be copied from a template)
- **Build**: `Makefile` (primary)
- **Core modules**: `grid.*, packet.*, radfield.*, nltepop.*`

### Configuration Templates
Choose appropriate artisoptions template for your simulation:
- `artisoptions_classic.h` - Classic mode simulations
- `artisoptions_nltenebular.h` - NLTE nebular phase
- `artisoptions_kilonova_lte.h` - Kilonova simulations in LTE

### Test Cases Available
CI tests 12 different scenarios - examples:
- `classicmode_1d_3dgrid` - Simple 1D classic mode (good for testing)
- `kilonova_1d_1dgrid` - 1D kilonova simulation
- `nebular_1d_3dgrid` - NLTE nebular phase simulation

## Troubleshooting

### Common Build Issues
1. **Missing artisoptions.h**: Copy from a template (see Building section)
2. **Missing dependencies**: Install MPI and GSL (see Prerequisites)
3. **CMake C++23 errors**: Use Makefile instead of CMake
4. **Cast function type warnings**: Normal with some MPI versions, can be ignored

### Common Runtime Issues
1. **"artis.pid exists"**: Delete `artis.pid` file if stale
2. **Missing input files**: Use test setup scripts to get proper inputs
3. **MPI errors**: Ensure proper MPI flags (`--oversubscribe --mca mpi_yield_when_idle 1`)

### Performance Notes
- Build time: ~30 seconds on modern hardware
- Test simulation time: Highly variable (minutes to hours)
- Memory usage: Scales with simulation size and MPI processes
- Always use MPI even for single-node runs

## Key Principles
- **NEVER CANCEL long-running builds or simulations** - set appropriate timeouts
- **Always validate changes** with at least one complete test scenario
- **Test with MPI** even for development (single process: `mpirun -np 1`)
- **Check artisoptions.h** is properly configured for your simulation mode
- **Scientific output validation** is critical - check MD5 sums when available

## Working with the CI
The GitHub Actions CI (`.github/workflows/ci.yml`) runs comprehensive tests:
- 12 different simulation scenarios
- Reproducible build settings
- MD5 checksum validation
- 120-minute timeout (reflects real runtime expectations)

When making changes, ensure your modifications don't break any existing test scenarios.
