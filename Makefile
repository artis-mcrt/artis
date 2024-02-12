.DEFAULT_GOAL := all

# place in architecture folder, e.g. build/arm64
BUILD_DIR = build/$(shell uname -m)

ifeq ($(MPI),)
	# MPI option not specified. set to true by default
	MPI := ON
endif
ifeq ($(MPI),ON)
	CXX = mpicxx
	CXXFLAGS += -DMPI_ON=true
	BUILD_DIR := $(BUILD_DIR)_mpi
else ifeq ($(MPI),OFF)
else
	$(error bad value for MPI option. Should be ON or OFF)
endif

ifeq ($(TESTMODE),ON)
else ifeq ($(TESTMODE),OFF)
else ifeq ($(TESTMODE),)
else
$(error bad value for testmode option. Should be ON or OFF)
endif

COMPILER_VERSION := $(shell $(CXX) --version)

ifneq '' '$(findstring clang,$(COMPILER_VERSION))'
  COMPILER_NAME := CLANG
else ifneq '' '$(findstring g++,$(COMPILER_VERSION))'
  COMPILER_NAME := GCC
else ifneq '' '$(findstring nvc++,$(COMPILER_VERSION))'
  COMPILER_NAME := NVHPC
else
  $(warning Unknown compiler)
  COMPILER_NAME := unknown
endif

$(info detected compiler is $(COMPILER_NAME))

ifeq ($(COMPILER_NAME),NVHPC)
	CXXFLAGS += -std=c++20
else
	CXXFLAGS += -std=c++23 -ftree-vectorize -flto=auto -Wunknown-pragmas -Wunused-macros
endif


CXXFLAGS += -fstrict-aliasing -ftrivial-auto-var-init=zero


ifeq ($(OPENMP),ON)
  BUILD_DIR := $(BUILD_DIR)_openmp

  ifeq ($(COMPILER_NAME),CLANG)
    CXXFLAGS += -Xpreprocessor -fopenmp
    LDFLAGS += -lomp
  else
    CXXFLAGS += -fopenmp
  endif
else ifeq ($(OPENMP),OFF)
else ifeq ($(OPENMP),)
else
  $(error bad value for openmp option. Should be ON or OFF)
endif

ifeq ($(STDPAR),ON)
  ifeq ($(OPENMP),ON)
    $(error cannot combine OPENMP and STDPAR)
  endif

  CXXFLAGS += -DSTDPAR_ON=true
  BUILD_DIR := $(BUILD_DIR)_stdpar

  ifeq ($(COMPILER_NAME),NVHPC)
	CXXFLAGS += -stdpar=gpu
  else ifeq ($(COMPILER_NAME),GCC)
    LDFLAGS += -ltbb
  else ifeq ($(COMPILER_NAME),CLANG)
    LDFLAGS += -Xlinker -debug_snapshot
  endif
else ifeq ($(STDPAR),OFF)
else ifeq ($(STDPAR),)
else
  $(error bad value for STDPAR option. Should be ON or OFF)
endif


ifeq ($(shell uname -s),Darwin)
# 	macOS

    ifeq ($(COMPILER_NAME),CLANG)
    #   fixes linking on macOS with gcc
	  LDFLAGS += -Wl,-ld_classic
    endif

	ifeq ($(shell uname -m),arm64)
#	 	On Arm, -mcpu combines -march and -mtune
		CXXFLAGS += -mcpu=native
	else
#		On x86, -march implies -mtune
		CXXFLAGS += -march=native
	endif

	CXXFLAGS += -fno-omit-frame-pointer
#	CXXFLAGS += -Rpass=loop-vectorize
#	CXXFLAGS += -Rpass-missed=loop-vectorize
#	CXXFLAGS += -Rpass-analysis=loop-vectorize

	# CXXFLAGS += -fopenmp-simd

	# enable OpenMP for Clang
	# CXXFLAGS += -Xpreprocessor -fopenmp -lomp

	# add -lprofiler for gperftools
	# LDFLAGS += $(LIB)
	# LDFLAGS += -lprofiler
	CXXFLAGS += $(shell pkg-config --cflags ompi)

else
	# sometimes the login nodes have slighty different CPUs
	# to the job nodes. Try to find the lowest common denominator here
	# to enable vector extensions
	# CXXFLAGS += -march=cascadelake
	# CXXFLAGS += -march=skylake-avx512
	# CXXFLAGS += -march=icelake-server

	# to get the current CPU architecture, run this:
	# g++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3
	ifneq (,$(shell hostname -A | grep gsi.de))
		# virgo has some AMD nodes and some Intel.
		# As of Feb 2024, the login nodes are zen3, and we select the same arch for jobs
		CXXFLAGS += -march=native
	else
		# for GitHub actions, checksums must match with different assigned CPUs, so avoid -march=native (use lowest common denominator)
		# update: all zenver3 now?

		CXXFLAGS += -march=native
	endif

endif

# GSL (GNU Scientific Library)
LDFLAGS += $(shell pkg-config --libs gsl)
CXXFLAGS += $(shell pkg-config --cflags gsl)
# GSL option 1: Use pkg-config or gsl-config to find GSL
#
# GSL option 2: Use default search paths to find GSL
# LDFLAGS += -lgsl -lgslcblas -lm

# Use GSL inline functions
CXXFLAGS += -DHAVE_INLINE -DGSL_C99_INLINE

ifeq ($(TESTMODE),ON)
	CXXFLAGS += -DTESTMODE=true -D_LIBCPP_DEBUG=0

	CXXFLAGS += -D_GLIBCXX_ASSERTIONS
	# CXXFLAGS += -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_BACKTRACE=1

	CXXFLAGS +=  -fno-omit-frame-pointer

	# CXXFLAGS += -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_EXTENSIVE
	CXXFLAGS += -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_DEBUG

	CXXFLAGS += -fsanitize=address,undefined

	BUILD_DIR := $(BUILD_DIR)_testmode
else
	# skip GSL range checking for better performance
	CXXFLAGS += -DTESTMODE=false -DGSL_RANGE_CHECK_OFF
endif

ifeq ($(OPTIMIZE),OFF)
	BUILD_DIR := $(BUILD_DIR)_optimizeoff
	CXXFLAGS += -O0
else
	# ifeq ($(TESTMODE),ON)
	# 	CXXFLAGS += -Og
	# else
		ifeq ($(FASTMATH),OFF)
			CXXFLAGS += -O3
			BUILD_DIR := $(BUILD_DIR)_nofastmath
		else
			CXXFLAGS += -Ofast

			ifeq ($(COMPILER_NAME),NVHPC)
				CXXFLAGS += -fast
			else
				CXXFLAGS += -ffast-math -funsafe-math-optimizations -fno-finite-math-only
			endif
		endif
	# endif
endif

CXXFLAGS += -Werror -Werror=undef -Winline -Wall -Wpedantic -Wredundant-decls -Wundef -Wno-unused-parameter -Wno-unused-function -Wno-inline -Wsign-compare


### use pg when you want to use gprof profiler
#CXXFLAGS = -g -pg -Wall -I$(INCLUDE)

# sn3d.cc and exspec.cc have main() defined
common_files := $(filter-out sn3d.cc exspec.cc, $(wildcard *.cc))

sn3d_files = sn3d.cc $(common_files)
sn3d_objects = $(addprefix $(BUILD_DIR)/,$(sn3d_files:.cc=.o))
sn3d_dep = $(sn3d_objects:%.o=%.d)

exspec_files = exspec.cc $(common_files)
exspec_objects = $(addprefix $(BUILD_DIR)/,$(exspec_files:.cc=.o))
exspec_dep = $(exspec_objects:%.o=%.d)

all: sn3d exspec

$(BUILD_DIR)/%.o: %.cc artisoptions.h Makefile
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -MD -MP -c $< -o $@

$(BUILD_DIR)/sn3d.o $(BUILD_DIR)/exspec.o: version.h artisoptions.h Makefile

check: $(sn3d_files)
	run-clang-tidy $(sn3d_files)

sn3d: $(sn3d_objects) artisoptions.h Makefile
	$(CXX) $(CXXFLAGS) $(sn3d_objects) $(LDFLAGS) -o sn3d
-include $(sn3d_dep)

sn3dwhole: version.h artisoptions.h Makefile
	$(CXX) $(CXXFLAGS) -g $(sn3d_files) $(LDFLAGS) -o sn3d

exspec: $(exspec_objects) artisoptions.h Makefile
	$(CXX) $(CXXFLAGS) $(exspec_objects) $(LDFLAGS) -o exspec
-include $(exspec_dep)

.PHONY: clean version.h TESTMODE TESTMODEON

version.h:
	@echo "constexpr const char* GIT_VERSION = \"$(shell git describe --dirty --always --tags)\";" > version.h
	@echo "constexpr const char* GIT_HASH = \"$(shell git rev-parse HEAD)\";" >> version.h
# requires git > 2.22
# @echo "constexpr const char* GIT_BRANCH = \"$(shell git branch --show)\";" >> version.h
	@echo "constexpr const char* GIT_BRANCH = \"$(shell git symbolic-ref --short HEAD 2>/dev/null || git rev-parse --short HEAD )\";" >> version.h
	@echo "constexpr const char* GIT_STATUS = \"$(shell git status --short)\";" >> version.h

clean:
	rm -rf sn3d exspec build *.o *.d
