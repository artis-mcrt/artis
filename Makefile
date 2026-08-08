# it's recommended that you add the following to your startup script:
# export MAKEFLAGS="--check-symlink-times --jobs=$(nproc --all)"
.DEFAULT_GOAL := all

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

$(info detected compiler is $(COMPILER_NAME) $(COMPILER_VERSION_NUMBER))
$(info detected CPU is $(CPU_ARCH))

# Use a custom build directory for each combination of compiler, CPU architecture, and options to avoid conflicts and ensure that the correct binaries are used
BUILD_DIR = build/$(COMPILER_NAME)-$(COMPILER_VERSION_NUMBER)_$(CPU_ARCH)

CXXFLAGS += -std=$(CXX_STD) $(ARCH_FLAGS) -Wall -Wextra -Wpedantic -Wredundant-decls -Wno-unused-parameter -Wsign-compare -Wshadow -isystem third_party

ifneq ($(COMPILER_NAME),nvhpc)
	CXXFLAGS += -Wunused-macros -Werror -Wextra-semi -Wno-unknown-pragmas -Wno-error=cast-function-type -MD -MP -Wno-unused-function
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
endif

ifeq ($(TESTMODE),ON)
	CXXFLAGS += -DTESTMODE=true

	CXXFLAGS += -fno-omit-frame-pointer -g

	CXXFLAGS += -D_GLIBCXX_ASSERTIONS
	# CXXFLAGS += -D_GLIBCXX_DEBUG
	# CXXFLAGS += -D_GLIBCXX_DEBUG_BACKTRACE

	# CXXFLAGS += -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST
	CXXFLAGS += -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_EXTENSIVE -DEIGEN_MAX_ALIGN_BYTES=0

	CXXFLAGS += -fsanitize=undefined,address

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

# sn3d.cc, exspec.cc, and unittests.cc have main() defined
common_files := $(filter-out sn3d.cc exspec.cc unittests.cc, $(wildcard *.cc))

sn3d_files = $(common_files) sn3d.cc
sn3d_objects = $(addprefix $(BUILD_DIR)/,$(sn3d_files:.cc=.o))
sn3d_dep = $(sn3d_objects:%.o=%.d)

exspec_files = $(common_files) exspec.cc
exspec_objects = $(addprefix $(BUILD_DIR)/,$(exspec_files:.cc=.o))
exspec_dep = $(exspec_objects:%.o=%.d)

unittests_files = $(common_files) unittests.cc
unittests_objects = $(addprefix $(BUILD_DIR)/,$(unittests_files:.cc=.o))
unittests_dep = $(unittests_objects:%.o=%.d)

.ONESHELL:
define version_h
constexpr const char* GIT_VERSION = \"$(shell git describe --dirty --always --tags)\";
constexpr const char* GIT_BRANCH = \"$(shell git symbolic-ref --short HEAD 2>/dev/null || git rev-parse --short HEAD )\";
constexpr const char* GIT_STATUS = \"$(shell git status --short)\";
endef

$(shell echo "$(version_h)" > version_tmp.h)
$(shell test -f version.h || touch version.h)

ifneq ($(shell cat version.h),$(shell cat version_tmp.h))
  $(info updating version.h)
  $(shell mv version_tmp.h version.h)
else
  $(shell rm version_tmp.h)
endif

$(shell mkdir -p $(BUILD_DIR))
$(info build directory: $(BUILD_DIR))

all: sn3d exspec

$(BUILD_DIR)/%.o: %.cc Makefile
	$(CXX) $(CXXFLAGS) -c $< -o $@

check: $(sn3d_files)
	run-clang-tidy $(sn3d_files)

$(BUILD_DIR)/sn3d: $(sn3d_objects)
	$(CXX) $(CXXFLAGS) $(sn3d_objects) $(LDFLAGS) -o $(BUILD_DIR)/sn3d
-include $(sn3d_dep)

sn3d: $(BUILD_DIR)/sn3d
	ln -sf $(BUILD_DIR)/sn3d sn3d

$(BUILD_DIR)/sn3dwhole: $(sn3d_files) version.h artisoptions.h Makefile
	$(CXX) $(CXXFLAGS) $(sn3d_files) $(LDFLAGS) -o $(BUILD_DIR)/sn3dwhole
-include $(sn3d_dep)

sn3dwhole: $(BUILD_DIR)/sn3dwhole
	ln -sf $(BUILD_DIR)/sn3dwhole sn3d

$(BUILD_DIR)/exspec: $(exspec_objects)
	$(CXX) $(CXXFLAGS) $(exspec_objects) $(LDFLAGS) -o $(BUILD_DIR)/exspec
-include $(exspec_dep)

exspec: $(BUILD_DIR)/exspec
	ln -sf $(BUILD_DIR)/exspec exspec

$(BUILD_DIR)/unittests: $(unittests_objects)
	$(CXX) $(CXXFLAGS) $(unittests_objects) $(LDFLAGS) -o $(BUILD_DIR)/unittests
-include $(unittests_dep)

unittests: $(BUILD_DIR)/unittests
	ln -sf $(BUILD_DIR)/unittests unittests

.PHONY: clean sn3d sn3dwhole exspec unittests

clean:
	rm -rf sn3d exspec unittests build *.o *.d
