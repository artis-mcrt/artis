# it's recommended that you add the following to your startup script:
# export MAKEFLAGS="--check-symlink-times --jobs=$(nproc --all)"
.DEFAULT_GOAL := all

# artisoptions.h is gitignored and must be supplied (normally as a symlink to a preset). Check up front
# to give one actionable message rather than a compiler error per translation unit. 'make clean' is exempt
ifneq ($(strip $(filter-out clean,$(if $(MAKECMDGOALS),$(MAKECMDGOALS),all))),)
  ifeq ($(wildcard artisoptions.h),)
    $(error artisoptions.h not found. Select a compile-time option preset, e.g. 'ln -s artisoptions_classic.h artisoptions.h'. Symlinking (rather than copying) keeps it in sync when the preset changes. The options are documented in artisoptions_doc.md)
  endif
  # without --check-symlink-times ('L' in the first dashless word of MAKEFLAGS), make compares timestamps
  # through the symlink, so a preset switch via ln -sf may not trigger a rebuild
  ifneq ($(shell test -L artisoptions.h && echo is_symlink),)
    ifeq (,$(findstring L,$(filter-out -%,$(firstword $(MAKEFLAGS)))))
      $(warning artisoptions.h is a symlink, but make was started without --check-symlink-times: switching presets with ln -sf may not trigger a rebuild. Recommended: export MAKEFLAGS="--check-symlink-times --jobs=$$(nproc --all)")
    endif
  endif
endif

$(info mpicxx version: $(shell mpicxx --showme:version 2> /dev/null))

ifeq ($(TESTMODE),ON)
else ifeq ($(TESTMODE),OFF)
else ifeq ($(TESTMODE),)
else
  $(error bad value for TESTMODE option. Should be ON or OFF)
endif

CXX := mpicxx
COMPILER_VERSION := $(shell $(CXX) --version)
COMPILER_VERSION_NUMBER := $(shell $(CXX) -dumpfullversion -dumpversion | head -n1)
$(info $(COMPILER_VERSION))
CXX_STD := c++26

ifeq ($(shell uname -m),arm64)
	ARCH_FLAGS = -mcpu=native -mtune=native
else
	ARCH_FLAGS = -march=native
endif

COMPILER_NAME := unknown
CPU_ARCH := unknown
ifneq '' '$(findstring HIP version,$(COMPILER_VERSION))'
	COMPILER_NAME := hipcc
	CXX_STD := c++23
	CXXFLAGS += -Wno-macro-redefined -Wno-unused-command-line-argument
else ifneq '' '$(findstring clang,$(COMPILER_VERSION))'
	COMPILER_NAME := clang
	CXXFLAGS += -Wunsafe-buffer-usage -Wno-unsafe-buffer-usage-in-libc-call -fsafe-buffer-usage-suggestions -Wno-unneeded-internal-declaration
	LDFLAGS += -Wno-unused-command-line-argument

	ifeq '' '$(findstring Apple,$(COMPILER_VERSION))'
		ifeq ($(if $(shell command -v lld),'true','false'), 'true')
			LDFLAGS += -fuse-ld=lld
		endif
	endif

define GET_ARCH_CMD
	mpicxx $(ARCH_FLAGS) -### -c -x c++ /dev/null 2>&1 \
	| tr ' ' '\n' \
	| awk '/-target-cpu/ {getline; gsub(/"/,""); print; exit}'
endef
	CPU_ARCH := $(shell $(GET_ARCH_CMD))

else ifneq (,$(or $(findstring g++,$(COMPILER_VERSION)),$(findstring gcc,$(COMPILER_VERSION))))
	COMPILER_NAME := gcc
	MIN_GCC_VERSION := 14
	ifeq ($(shell expr $(COMPILER_VERSION_NUMBER) \< $(MIN_GCC_VERSION)),1)
$(warning WARNING: Detected GCC version $(COMPILER_VERSION_NUMBER) but minimum supported version is $(MIN_GCC_VERSION))
	endif
	CXXFLAGS += -Wno-psabi -Wno-interference-size
# 	CXXFLAGS += -Wsuggest-attribute=pure -Wsuggest-attribute=const
	CPU_ARCH := $(shell $(CXX) $(ARCH_FLAGS) -Q --help=target  | grep -- '-mtune=  ' | cut -f3)

else ifneq '' '$(findstring nvc++,$(COMPILER_VERSION))'
	COMPILER_NAME := nvhpc
	CXX_STD := c++23
	# to use the pixi installed libstdc++
	CXXFLAGS += -Minfo=accel
# 	CXXFLAGS += --gcc-toolchain=$(PWD)/.pixi/envs/default -Wl,-rpath,$(PWD)/.pixi/envs/default/lib

	ifneq ($(GCCTOOLCHAIN),)
		CXXFLAGS += --gcc-toolchain=$(GCCTOOLCHAIN)
	else ifneq (,$(shell hostname -A | grep .cosma.))
		CXXFLAGS += --gcc-toolchain=/cosma/local/gcc/14.1.0/
	else ifneq (,$(shell hostname -A | grep heavymetalgb10.mp.qub.ac.uk))
		CXXFLAGS += --gcc-toolchain=$(shell which g++-14)
	endif
	CPU_ARCH := $(shell g++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)
else
$(warning Unknown compiler)
	COMPILER_NAME := unknown
endif

# A host can have more than one gcc installation, and the newest one can have no libstdc++ headers.
# The clang family then gives -Wgcc-install-dir-libstdcxx, because a future release will select a
# different installation. The ROCm container is such a host. The build cannot correct the host, so
# the warning stays visible, but -Werror must not make it an error. An older compiler does not know
# the option, and gives -Wunknown-warning-option, so first ask the compiler for that message.
ifneq (,$(filter $(COMPILER_NAME),hipcc clang))
	ifeq (,$(shell $(CXX) -Wno-error=gcc-install-dir-libstdcxx -fsyntax-only -x c++ /dev/null 2>&1 | grep -m1 'unknown warning option'))
		CXXFLAGS += -Wno-error=gcc-install-dir-libstdcxx
	endif
endif

$(info detected compiler is $(COMPILER_NAME) $(COMPILER_VERSION_NUMBER))
$(info detected CPU is $(CPU_ARCH))

# Use a custom build directory for each combination of compiler, CPU architecture, and options to avoid conflicts and ensure that the correct binaries are used
BUILD_DIR = build/$(COMPILER_NAME)-$(COMPILER_VERSION_NUMBER)_$(CPU_ARCH)

CXXFLAGS += -std=$(CXX_STD) $(ARCH_FLAGS) -Wall -Wextra -Wpedantic -Wredundant-decls -Wno-unused-parameter -Wsign-compare -Wshadow -isystem third_party

# generate and use .d header dependency files, so that header edits trigger recompilation of the
# objects that include them (every compiler, including nvc++, supports these GCC-style options)
CXXFLAGS += -MD -MP

ifneq ($(COMPILER_NAME),nvhpc)
	CXXFLAGS += -Wunused-macros -Werror -Wextra-semi -Wno-unknown-pragmas -Wno-error=cast-function-type -Wno-unused-function
endif

ifeq ($(REPRODUCIBLE),ON)
	CXXFLAGS += -DREPRODUCIBLE=true -ffp-contract=off -DEIGEN_DONT_VECTORIZE
	BUILD_DIR := $(BUILD_DIR)_reproducible
	FASTMATH := OFF
else ifeq ($(REPRODUCIBLE),OFF)
else ifeq ($(REPRODUCIBLE),)
else
  $(error bad value for REPRODUCIBLE option. Should be ON or OFF)
endif

# CXXFLAGS += -DUSE_SIMPSON_INTEGRATOR

ifeq ($(GPU),ON)
	CXXFLAGS += -DGPU_ON -DUSE_SIMPSON_INTEGRATOR -U_GLIBCXX_ASSERTIONS
	BUILD_DIR := $(BUILD_DIR)_gpu
else ifeq ($(GPU),OFF)
else ifeq ($(GPU),)
else
    $(error bad value for GPU option. Should be ON or OFF)
endif

ifeq ($(OPENMP),ON)
	ifeq ($(STDPAR),ON)
		$(error cannot combine OPENMP and STDPAR)
	endif
	BUILD_DIR := $(BUILD_DIR)_openmp

	ifeq ($(COMPILER_NAME),nvhpc)
		ifeq ($(GPU),ON)
			CXXFLAGS += -mp=gpu
		else
			CXXFLAGS += -mp
		endif
	else ifeq ($(COMPILER_NAME),clang)
		CXXFLAGS += -Xpreprocessor -fopenmp
		LDFLAGS += -lomp
	else ifeq ($(COMPILER_NAME),gcc)
		CXXFLAGS += -fopenmp
	endif

else ifeq ($(OPENMP),OFF)
else ifeq ($(OPENMP),)
else
    $(error bad value for OPENMP option. Should be ON or OFF)
endif

ifeq ($(STDPAR),ON)
	CXXFLAGS += -DSTDPAR_ON=true
	BUILD_DIR := $(BUILD_DIR)_stdpar

	ifeq ($(COMPILER_NAME),nvhpc)
		ifeq ($(GPU),ON)
			CXXFLAGS += -stdpar=gpu
		else
			CXXFLAGS += -stdpar=multicore
		endif
	else ifeq ($(COMPILER_NAME),hipcc)
		CXXFLAGS += -fexperimental-library
		ifeq ($(GPU),ON)
			# MI300
			CXXFLAGS += --offload-arch=gfx942 -fgpu-rdc --hipstdpar
		endif
	else ifeq ($(COMPILER_NAME),clang)
		CXXFLAGS += -fexperimental-library
		# LDFLAGS += -ltbb
		# LDFLAGS += -Xlinker -debug_snapshot

	else ifeq ($(COMPILER_NAME),gcc)
		LDFLAGS += -ltbb
	endif

else ifeq ($(STDPAR),OFF)
else ifeq ($(STDPAR),)
else
  $(error bad value for STDPAR option. Should be ON or OFF)
endif

ifeq ($(COMPILER_NAME),nvhpc)
	ifeq ($(GPU),ON)
			CXXFLAGS += -gpu=mem:unified -gpu=ccnative
			CXXFLAGS += -Minfo=stdpar,accel
# 			CXXFLAGS += -gpu=debug -g
# 			CXXFLAGS += -gpu=maxregcount:64
$(info to debug, run mpirun -np 1 compute-sanitizer --tool memcheck ./artis/sn3d)
	endif
endif

ifeq ($(shell uname -s),Darwin)
# 	macOS
	CXXFLAGS += -fno-omit-frame-pointer -g
	# gcc
	# CXXFLAGS += -fopt-info-vec-missed
	# clang
	# CXXFLAGS += -Rpass=loop-vectorize
	# CXXFLAGS += -Rpass-missed=loop-vectorize
	# CXXFLAGS += -Rpass-analysis=loop-vectorize
endif

ifneq ($(MAX_NODE_SIZE),)
	CXXFLAGS += -DMAX_NODE_SIZE=$(MAX_NODE_SIZE)
	# must be in the build directory name like every other option that changes CXXFLAGS
	BUILD_DIR := $(BUILD_DIR)_maxnodesize$(MAX_NODE_SIZE)
endif

# The flags of TESTMODE=ON. compile_commands.json always uses them, also when the build does
# not, so that clangd and clang-tidy examine the body of each assert_testmodeonly() macro.
TESTMODE_CXXFLAGS := -DTESTMODE=true

TESTMODE_CXXFLAGS += -fno-omit-frame-pointer -g

TESTMODE_CXXFLAGS += -D_GLIBCXX_ASSERTIONS
# TESTMODE_CXXFLAGS += -D_GLIBCXX_DEBUG
# TESTMODE_CXXFLAGS += -D_GLIBCXX_DEBUG_BACKTRACE

# TESTMODE_CXXFLAGS += -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST
TESTMODE_CXXFLAGS += -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_EXTENSIVE -DEIGEN_MAX_ALIGN_BYTES=0

TESTMODE_CXXFLAGS += -fsanitize=undefined,address

ifeq ($(TESTMODE),ON)
	CXXFLAGS += $(TESTMODE_CXXFLAGS)

	BUILD_DIR := $(BUILD_DIR)_testmode
endif

ifeq ($(OPTIMIZE),OFF)
	BUILD_DIR := $(BUILD_DIR)_optimizeoff
	CXXFLAGS += -O0 -DEIGEN_FAST_MATH=0
else
	CXXFLAGS += -O3

	ifeq ($(COMPILER_NAME),clang)
		CXXFLAGS += -flto=thin
	else ifeq ($(COMPILER_NAME),gcc)
		CXXFLAGS += -flto=auto -pipe
	endif

	ifeq ($(FASTMATH),OFF)
		BUILD_DIR := $(BUILD_DIR)_nofastmath
		CXXFLAGS += -DEIGEN_FAST_MATH=0
	else
		ifeq ($(COMPILER_NAME),nvhpc)
			CXXFLAGS += -fast
		else
			CXXFLAGS += -ffast-math -funsafe-math-optimizations -fno-finite-math-only
		endif
	endif
endif

# profile-guided optimisation (gcc and clang):
#   1. build with 'make PGO=GENERATE', then run a representative simulation to collect profile data
#      (written into the build directory)
#   2. clang only: merge the raw profiles with
#      llvm-profdata merge -output=<build dir>/default.profdata <build dir>/*.profraw
#   3. rebuild with 'make PGO=USE' (switching phase automatically discards the instrumented objects)
ifneq ($(PGO),)
  ifeq ($(filter $(PGO),GENERATE USE),)
    $(error bad value for PGO option. Should be GENERATE or USE)
  endif
  ifeq ($(filter $(COMPILER_NAME),gcc clang),)
    $(error PGO is only supported for gcc and clang)
  endif

  BUILD_DIR := $(BUILD_DIR)_pgo

  # the profile output path is embedded in the binary and resolved against the runtime working
  # directory, so it must be absolute for profiles to land in the build directory
  PGO_PROFDIR := $(abspath $(BUILD_DIR))

  ifeq ($(PGO),GENERATE)
      CXXFLAGS += -fprofile-generate="$(PGO_PROFDIR)" -fprofile-update=prefer-atomic
  else
    ifeq ($(COMPILER_NAME),clang)
      $(shell llvm-profdata merge -output="$(PGO_PROFDIR)/default.profdata" "$(PGO_PROFDIR)/"*.profraw)
      CXXFLAGS += -fprofile-use="$(PGO_PROFDIR)/default.profdata"
#       CXXFLAGS += -Wno-profile-instr-unprofiled -Wno-profile-instr-out-of-date
    else
      CXXFLAGS += -fprofile-use="$(PGO_PROFDIR)" -fprofile-correction
    endif
  endif

  # instrumented and optimised objects cannot be mixed, so start clean when the phase changes
  ifneq ($(shell cat $(BUILD_DIR)/pgo_phase 2>/dev/null),$(PGO))
    $(shell mkdir -p $(BUILD_DIR))
    $(shell rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.d $(BUILD_DIR)/sn3d $(BUILD_DIR)/sn3dwhole $(BUILD_DIR)/exspec)
    $(shell echo $(PGO) > $(BUILD_DIR)/pgo_phase)
  endif
endif

.ONESHELL:
define version_cc
extern const char* const GIT_VERSION = \"$(shell git describe --dirty --always --tags)\";
extern const char* const GIT_BRANCH = \"$(shell git symbolic-ref --short HEAD 2>/dev/null || git rev-parse --short HEAD )\";
extern const char* const GIT_STATUS = \"$(shell git status --short)\";
endef

# the git metadata changes at almost every git operation. version.h declares the strings and stays
# the same, and this generated file holds the definitions. Only this small file compiles again,
# and not sn3d.cc and exspec.cc
$(shell echo "$(version_cc)" > .version_tmp.cc)
$(shell test -f version.cc || touch version.cc)

ifneq ($(shell cat version.cc),$(shell cat .version_tmp.cc))
  $(info updating version.cc)
  $(shell mv .version_tmp.cc version.cc)
else
  $(shell rm .version_tmp.cc)
endif

# the sources of all three programs. sn3d.cc, exspec.cc, and unittests.cc have main() defined, so
# they are not here. Add a new source file to this list. A wildcard is not correct here, because
# make caches the content of the directory and does not find a version.cc that this run made
common_files := \
	chargetransfer.cc \
	decay.cc \
	gammapkt.cc \
	grid.cc \
	input.cc \
	kpkt.cc \
	ltepop.cc \
	macroatom.cc \
	mpi_logging.cc \
	nltepop.cc \
	nonthermal.cc \
	packet.cc \
	radfield.cc \
	ratecoeff.cc \
	rpkt.cc \
	spectrum_lightcurve.cc \
	stats.cc \
	thermalbalance.cc \
	update_grid.cc \
	update_packets.cc \
	version.cc \
	vpkt.cc

sn3d_files = $(common_files) sn3d.cc
sn3d_objects = $(addprefix $(BUILD_DIR)/,$(sn3d_files:.cc=.o))
sn3d_dep = $(sn3d_objects:%.o=%.d)

exspec_files = $(common_files) exspec.cc
exspec_objects = $(addprefix $(BUILD_DIR)/,$(exspec_files:.cc=.o))
exspec_dep = $(exspec_objects:%.o=%.d)

unittests_files = $(common_files) unittests.cc
unittests_objects = $(addprefix $(BUILD_DIR)/,$(unittests_files:.cc=.o))
unittests_dep = $(unittests_objects:%.o=%.d)

$(shell mkdir -p $(BUILD_DIR))
$(info build directory: $(BUILD_DIR))

all: sn3d exspec

# artisoptions.h is explicit so the most consequential header dependency does not rely on .d generation
$(BUILD_DIR)/%.o: %.cc artisoptions.h Makefile
	$(CXX) $(CXXFLAGS) -c $< -o $@

# runs the same check as CI, over sn3d, exspec, and unittests
check: compile_commands.json
	@test -f compile_commands.json || { echo 'error: no compile_commands.json. python3 makes it'; exit 1; }
	run-clang-tidy

# clangd and clang-tidy read compile_commands.json. Each entry must name the compiler that
# $(CXX) calls, because clang-tidy finds its resource directory from that name. See
# "Linting and formatting" in AGENTS.md. Open MPI and MPICH both print the full command
# with -show. Every build reads these variables, so ':=' keeps the probe to one call.
CDB_PYTHON := $(shell command -v python3)
CDB_SHOW := $(shell $(CXX) -show 2>/dev/null)
CDB_CXX := $(shell command -v $(firstword $(CDB_SHOW)))
CDB_INCFLAGS := $(patsubst -I%,-isystem%,$(filter -I%,$(CDB_SHOW)))
CDB_SRC = $(sort $(sn3d_files) $(exspec_files) $(unittests_files))
# the database always uses TESTMODE=ON, so add those flags when the build does not have them
CDB_CXXFLAGS = $(CXXFLAGS) $(if $(filter ON,$(TESTMODE)),,$(TESTMODE_CXXFLAGS))

# Each build writes the database. The content comes from the make options and not from a file
# that make can examine, so the target is .PHONY. The recipe replaces the file only when the
# content changes.
# clangd therefore reads the files again only after a real change.
# make clean keeps the database, because clangd reads it and compiles nothing.
# The build needs only a compiler and MPI, so an absent python3 gives a message and no error.
# The clang tools need python3, and the check target stops if the database is absent.
.PHONY: compile_commands.json
compile_commands.json:
ifeq ($(CDB_PYTHON),)
	@echo '$@: python3 is absent, so make writes no compilation database'
else ifeq ($(CDB_CXX),)
	@echo '$@: "$(CXX) -show" does not name the compiler, so make writes no compilation database'
else
	@CDB_CXX='$(CDB_CXX)' CDB_FLAGS='$(CDB_INCFLAGS) $(CDB_CXXFLAGS)' CDB_SRC='$(CDB_SRC)' \
		$(CDB_PYTHON) -c 'import json,os,shlex,sys; args=[os.environ["CDB_CXX"],*shlex.split(os.environ["CDB_FLAGS"])]; json.dump([{"directory":os.getcwd(),"arguments":[*args,"-c",src],"file":src} for src in os.environ["CDB_SRC"].split()],sys.stdout,indent=1)' > $@.tmp
	@if cmp -s $@.tmp $@; then rm -f $@.tmp; else mv $@.tmp $@; echo '$@: $(words $(CDB_SRC)) entries for $(CDB_CXX)'; fi
endif

$(BUILD_DIR)/sn3d: $(sn3d_objects)
	$(CXX) $(CXXFLAGS) $(sn3d_objects) $(LDFLAGS) -o $(BUILD_DIR)/sn3d
-include $(sn3d_dep)

sn3d: $(BUILD_DIR)/sn3d compile_commands.json
	ln -sf $(BUILD_DIR)/sn3d sn3d

$(BUILD_DIR)/sn3dwhole: $(sn3d_files) version.h artisoptions.h Makefile
	$(CXX) $(CXXFLAGS) $(sn3d_files) $(LDFLAGS) -o $(BUILD_DIR)/sn3dwhole
-include $(sn3d_dep)

sn3dwhole: $(BUILD_DIR)/sn3dwhole compile_commands.json
	ln -sf $(BUILD_DIR)/sn3dwhole sn3d

$(BUILD_DIR)/exspec: $(exspec_objects)
	$(CXX) $(CXXFLAGS) $(exspec_objects) $(LDFLAGS) -o $(BUILD_DIR)/exspec
-include $(exspec_dep)

exspec: $(BUILD_DIR)/exspec compile_commands.json
	ln -sf $(BUILD_DIR)/exspec exspec

$(BUILD_DIR)/unittests: $(unittests_objects)
	$(CXX) $(CXXFLAGS) $(unittests_objects) $(LDFLAGS) -o $(BUILD_DIR)/unittests
-include $(unittests_dep)

unittests: $(BUILD_DIR)/unittests compile_commands.json
	ln -sf $(BUILD_DIR)/unittests unittests

.PHONY: clean sn3d sn3dwhole exspec unittests check

clean:
	rm -rf sn3d exspec unittests build *.o *.d
