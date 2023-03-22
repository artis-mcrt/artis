.DEFAULT_GOAL := all

# place in architecture folder, e.g. build/arm64
BUILD_DIR = build/$(shell uname -m)

CXXFLAGS += -std=c++20 -fstrict-aliasing -ftree-vectorize -g -flto=auto -Werror -Werror=undef
# CXXFLAGS += -Wunreachable-code

ifeq ($(shell uname -s),Darwin)
# 	macOS

	ifeq ($(shell uname -m),arm64)
#	 	On Arm, -mcpu combines -march and -mtune
		CXXFLAGS += -mcpu=native
	else
#		On x86, -march implies -mtune
		CXXFLAGS += -march=native
	endif

#	CXXFLAGS += -Rpass=loop-vectorize
#	CXXFLAGS += -Rpass-missed=loop-vectorize
#	CXXFLAGS += -Rpass-analysis=loop-vectorize

	# CXXFLAGS += -fopenmp-simd

	# enable OpenMP for Clang
	# CXXFLAGS += -Xpreprocessor -fopenmp -lomp

	# add -lprofiler for gperftools
	# LDFLAGS += $(LIB)
	# LDFLAGS += -lprofiler

else ifeq ($(USER),localadmin_ccollins)
	# CXX = c++
	LDFLAGS= -lgsl -lgslcblas -lm -I/home/localadmin_ccollins/gsl/include
	INCLUDE = /home/localadmin_ccollins/gsl/include
	LIB = /home/localadmin_ccollins/gsl/lib
	CXXFLAGS += -g -I$(INCLUDE)
	LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm
	CXXFLAGS += -std=c++17 -Wstrict-aliasing -fstrict-aliasing #-fopenmp=libomp

else
	# sometimes the login nodes have slighty different CPUs
	# to the job nodes. Try to find the lowest common denominator here
	# to enable vector extensions
	# CXXFLAGS += -march=haswell
	# CXXFLAGS += -march=skylake-avx512

	# to get the current CPU architecture, run this:
	# g++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3
	ifneq (,$(findstring juwels,$(HOSTNAME)))
		CXXFLAGS += -march=skylake-avx512
	else ifneq (,$(findstring lxbk,$(HOSTNAME)))
		# virgo has some AMD nodes (znver1 arch) and some Intel
		CXXFLAGS += -march=haswell
	endif

endif

# GSL (GNU Scientific Library)
# GSL option 1: Use pkg-config to find GSL
LDFLAGS += $(shell pkg-config --libs gsl)
CXXFLAGS += $(shell pkg-config --cflags gsl)
#
# GSL option 2: Use default search paths to find GSL
# LDFLAGS += -lgsl -lgslcblas -lm

# Use GSL inline functions
CXXFLAGS += -DHAVE_INLINE -DGSL_C99_INLINE

ifeq ($(TESTMODE),ON)
	CXXFLAGS += -DTESTMODE=true -O3 -DLIBCXX_ENABLE_DEBUG_MODE
	CXXFLAGS += -fsanitize=address -fno-omit-frame-pointer -fno-common
	BUILD_DIR := $(BUILD_DIR)_testmode
else
	# skip array range checking for better performance and use optimizations
	CXXFLAGS += -DTESTMODE=false -DGSL_RANGE_CHECK_OFF -O3
endif

CXXFLAGS += -Winline -Wall -Wpedantic -Wredundant-decls -Wundef -Wno-unused-parameter -Wno-unused-function -Wstrict-aliasing -Wno-inline

ifeq ($(MPI),ON)
else ifeq ($(MPI),OFF)
else ifeq ($(MPI),)
	# MPI option not specified. set to true by default
	MPI := ON
else
$(error bad value for MPI option. Should be ON or OFF)
endif

ifeq ($(TESTMODE),ON)
else ifeq ($(TESTMODE),OFF)
else ifeq ($(TESTMODE),)
else
$(error bad value for testmode option. Should be ON or OFF)
endif

ifeq ($(MPI),ON)
	CXX = mpicxx
	CXXFLAGS += -DMPI_ON=true
	BUILD_DIR := $(BUILD_DIR)_mpi
endif

ifeq ($(OPENMP),ON)
	CXXFLAGS += -Xpreprocessor
	CXXFLAGS += -fopenmp
	LDFLAGS += -lomp
	BUILD_DIR := $(BUILD_DIR)_openmp
endif

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

sn3d: $(sn3d_objects)
	$(CXX) $(CXXFLAGS) $(sn3d_objects) $(LDFLAGS) -o sn3d
#	$(LINK.cpp) $(filter %.o,$^) -o $@
-include $(sn3d_dep)

sn3dwhole: version.h
	$(CXX) $(CXXFLAGS) $(sn3d_files) $(LDFLAGS) -o sn3d

$(BUILD_DIR)/%.o: %.cc artisoptions.h Makefile
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -MD -MP -c $< -o $@

$(BUILD_DIR)/sn3d.o $(BUILD_DIR)/exspec.o: version.h artisoptions.h Makefile

exspec: $(exspec_objects)
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
	rm -rf sn3d exspec build version.h *.o *.d
