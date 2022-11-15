GIT_VERSION := $(shell git describe --dirty --always --tags)
GIT_HASH := $(shell git rev-parse HEAD)
GIT_BRANCH := $(shell git branch | sed -n '/\* /s///p')
.DEFAULT_GOAL := all
SYSNAME := $(shell uname -s)

BUILD_DIR = build/$(shell uname -m)

CXXFLAGS += -std=c++20 -fstrict-aliasing -ftree-vectorize -g -flto=auto

ifeq ($(SYSNAME),Darwin)
	# macOS

	# march=native also works on Apple Silicon with Clang >=15
	CXXFLAGS += -march=native
#	CXXFLAGS += -Rpass=loop-vectorize
#	CXXFLAGS += -Rpass-missed=loop-vectorize
#	CXXFLAGS += -Rpass-analysis=loop-vectorize

	# ifeq ($(shell uname -m),arm64)
	# 	CXXFLAGS += -mcpu=apple-m1
	# else
	# 	CXXFLAGS += -march=native
	# endif

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
	CXXFLAGS += -std=c++17 -march=native -Wstrict-aliasing -fstrict-aliasing #-fopenmp=libomp

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
	CXXFLAGS += -DTESTMODE=true -O3
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
	CXXFLAGS += -DMPI_ON
	BUILD_DIR := $(BUILD_DIR)_mpi
endif

ifeq ($(OPENMP),ON)
	CXXFLAGS += -Xpreprocessor
	CXXFLAGS += -fopenmp
	LDFLAGS += -lomp
	BUILD_DIR := $(BUILD_DIR)_openmp
endif

sn3dcuda sn3dcudawhole: LDFLAGS += -lcudart
sn3dcuda sn3dcudawhole: CXXFLAGS += -DCUDA_ENABLED=true

# Gadi
ifneq (,$(findstring gadi,$(HOSTNAME)))
	# Tesla V100
	CUDA_NVCC_FLAGS += -arch=sm_70 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_70,code=compute_70
	CXX = mpic++
	# CXX = icpc
	# CXXFLAGS += -qopenmp
	INCLUDE += -I/home/120/ljs120/cuda_samples/common/inc
endif

# CXXFLAGS += -std=c++17
# CXXFLAGS += -fPIC -shared
# CUDA_NVCC_FLAGS += -Xcompiler -fPIC -shared -rdc=true
CUDA_NVCC_FLAGS += -ccbin=$(CXX) -std=c++17 -O3 -use_fast_math -Xcompiler "$(CXXFLAGS)" -rdc=true --expt-relaxed-constexpr
# CUDA_NVCC_FLAGS += -G -g

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

sn3d: version.h artisoptions.h $(sn3d_objects) Makefile
	$(CXX) $(CXXFLAGS) $(sn3d_objects) $(LDFLAGS) -o sn3d
#	$(LINK.cpp) $(filter %.o,$^) -o $@
-include $(sn3d_dep)

sn3dwhole: version.h
	$(CXX) $(CXXFLAGS) $(sn3d_files) $(LDFLAGS) -o sn3d

sn3dcudawhole: version.h
	nvcc -x cu $(CUDA_NVCC_FLAGS) $(INCLUDE) $(LDFLAGS) $(sn3d_files) -o sn3dcuda

sn3dcuda: version.h $(sn3d_objects)
	nvcc --gpu-architecture=sm_70 --device-link $(sn3d_objects) --output-file gpucode.o
	$(CXX) $(CXXFLAGS) gpucode.o $(INCLUDE) -lcudadevrt $(LDFLAGS) $(sn3d_objects) -o sn3dcuda

# %.o: %.cc
# 	nvcc -x cu $(CUDA_NVCC_FLAGS) $(INCLUDE) --device-c $< -c

$(BUILD_DIR)/%.o: %.cc Makefile artisoptions.h
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -MD -MP -c $< -o $@

exspec: version.h artisoptions.h $(exspec_objects) Makefile
	$(CXX) $(CXXFLAGS) $(exspec_objects) $(LDFLAGS) -o exspec
-include $(exspec_dep)

.PHONY: clean version.h TESTMODE TESTMODEON

version.h:
	@echo "#define GIT_VERSION \"$(GIT_VERSION)\"" > version.h
	@echo "#define GIT_HASH \"$(GIT_HASH)\"" >> version.h
	@echo "#define GIT_BRANCH \"$(GIT_BRANCH)\"" >> version.h

clean:
	rm -rf sn3d exspec build version.h *.o *.d

