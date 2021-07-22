GIT_VERSION := $(shell git describe --dirty --always --tags)
GIT_HASH := $(shell git rev-parse HEAD)
GIT_BRANCH := $(shell git branch | sed -n '/\* /s///p')
.DEFAULT_GOAL := all
SYSNAME := $(shell uname -s)

ifeq ($(SYSNAME),Darwin)
	# macOS

	# CXX = c++
	CXXFLAGS += -std=c++17 -march=native -fstrict-aliasing -ftree-vectorize -flto

	CXXFLAGS += -Winline -Wall -Wextra -Wredundant-decls -Wundef -Wno-unused-parameter -Wno-unused-function -Wstrict-aliasing

	# CXXFLAGS += -fopenmp-simd
	# CXXFLAGS += -fvectorize

	# enable OpenMP for Clang
	# CXXFLAGS += -Xpreprocessor -fopenmp -lomp

	# in GCC, -Wmisleading-indentation will be useful
	# also -fopenmp after -I$(INCLUDE)
	# maybe  -fopt-info-vec-missed
	#  -fwhole-program
	# add -lprofiler for gperftools
	LDFLAGS += $(LIB)
	# LDFLAGS += -lprofiler

else ifneq (,$(findstring kelvin,$(HOSTNAME)))
	# QUB Kelvin cluster
	# needs
	#  mpi/openmpi/1.8.5/gcc-4.4.7
	#  compilers/gcc/system(default)
	#  libs/gsl/1.16/gcc-4.4.7

	CXX = mpicxx
	CXXFLAGS += -std=c++17 -mcmodel=medium #-fopenmp=libomp
sn3d: CXXFLAGS += -DMPI_ON

else ifneq (, $(shell which mpicxx))
	# any other system which has mpicxx available (Juwels, Cambridge, Gadi, etc)

	CXX = mpicxx
	CXXFLAGS += -std=c++17 -march=native #-fopenmp=libomp

sn3d sn3dcuda: CXXFLAGS += -DMPI_ON

else ifeq ($(USER),localadmin_ccollins)
	# CXX = c++
	LDFLAGS= -lgsl -lgslcblas -lm -I/home/localadmin_ccollins/gsl/include
	INCLUDE = /home/localadmin_ccollins/gsl/include
	LIB = /home/localadmin_ccollins/gsl/lib
	CXXFLAGS += -g -I$(INCLUDE)
	LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm
	CXXFLAGS += -std=c++17 -march=native -Wstrict-aliasing -fstrict-aliasing #-fopenmp=libomp

else
	# CXX = c++
	# CXX = icpc
	CXXFLAGS += -std=c++17 -march=native -Wstrict-aliasing -fstrict-aliasing #-fopenmp=libomp
endif


# ** GSL (GNU Scientific Library) **

# option 1: Use pkg-config to find GSL and use dynamic linking
LDFLAGS += $(shell pkg-config --libs gsl)
CXXFLAGS += $(shell pkg-config --cflags gsl)
#
# option 2: Use compiler default search paths to find GSL and use dynamic linking
# LDFLAGS += -lgsl -lgslcblas -lm
#
# option 3: Specify the path to libgsl.a and libgslclas.a and use static linking (GSL needed to compile but not to run)
# CXXFLAGS += /usr/local/Cellar/gsl/2.6/lib/libgsl.a
# CXXFLAGS += /usr/local/Cellar/gsl/2.6/lib/libgslcblas.a

# Use GSL inline functions
CXXFLAGS += -DHAVE_INLINE -DGSL_C99_INLINE

ifeq ($(TESTMODE),ON)
	CXXFLAGS += -DTESTMODE=true -O3 -g -fsanitize=address -fno-omit-frame-pointer
else
	# skip array range checking for better performance and use optimizations
	CXXFLAGS += -DTESTMODE=false -DGSL_RANGE_CHECK_OFF -O3 -flto
endif

sn3dmpi: CXX = mpicxx
sn3dmpi: CXXFLAGS += -DMPI_ON
sn3dmpi: sn3d

sn3dopenmp: CXXFLAGS += -Xpreprocessor
sn3dopenmp: CXXFLAGS += -fopenmp
sn3dopenmp: LDFLAGS += -lomp
sn3dopenmp: sn3d

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
CUDA_NVCC_FLAGS += -ccbin=$(CXX) -std=c++17 -use_fast_math -Xcompiler "$(CXXFLAGS)" -rdc=true --expt-relaxed-constexpr
# CUDA_NVCC_FLAGS += -G -g

### use pg when you want to use gprof profiler
#CXXFLAGS = -g -pg -Wall -I$(INCLUDE)

sn3d_files = sn3d.cc atomic.cc boundary.cc decay.cc emissivities.cc gamma.cc globals.cc grey_emissivities.cc grid.cc gsl_managed.cc input.cc kpkt.cc light_curve.cc ltepop.cc macroatom.cc md5.cc nltepop.cc nonthermal.cc packet_init.cc photo_electric.cc polarization.cc radfield.cc ratecoeff.cc rpkt.cc stats.cc thermalbalance.cc update_grid.cc update_packets.cc vectors.cc vpkt.cc

sn3d_objects = sn3d.o atomic.o boundary.o decay.o emissivities.o gamma.o globals.o grey_emissivities.o grid.o gsl_managed.o input.o kpkt.o light_curve.o ltepop.o macroatom.o md5.o nltepop.o nonthermal.o packet_init.o photo_electric.o polarization.o radfield.o ratecoeff.o rpkt.o stats.o thermalbalance.o update_grid.o update_packets.o vectors.o vpkt.o

exspec_files = exspec.cc atomic.cc boundary.cc decay.cc emissivities.cc gamma.cc globals.cc grey_emissivities.cc grid.cc gsl_managed.cc input.cc kpkt.cc light_curve.cc ltepop.cc macroatom.cc md5.cc nltepop.cc nonthermal.cc packet_init.cc photo_electric.cc polarization.cc radfield.cc ratecoeff.cc rpkt.cc spectrum.cc stats.cc thermalbalance.cc update_grid.cc update_packets.cc vectors.cc vpkt.cc

all: sn3d exspec

sn3d: version $(sn3d_objects)
	$(CXX) $(CXXFLAGS) $(sn3d_objects) $(LDFLAGS) -o sn3d
#	$(CXX) $(CXXFLAGS) $(sn3d_files) $(LDFLAGS) -o sn3d

sn3dcudawhole: version
	nvcc -x cu $(CUDA_NVCC_FLAGS) $(INCLUDE) $(LDFLAGS) $(sn3d_files) -o sn3d

sn3dcuda: version $(sn3d_objects)
	nvcc --gpu-architecture=sm_70 --device-link $(sn3d_objects) --output-file gpucode.o
	$(CXX) $(CXXFLAGS) gpucode.o $(INCLUDE) -lcudadevrt $(LDFLAGS) $(sn3d_objects) -o sn3d

# %.o: %.cc
# 	nvcc -x cu $(CUDA_NVCC_FLAGS) $(INCLUDE) --device-c $< -c

exspec: version
	$(CXX) $(CXXFLAGS) -DDO_EXSPEC $(exspec_files) $(LDFLAGS) -o exspec

.PHONY: clean version

version:
	@echo "#define GIT_VERSION \"$(GIT_VERSION)\"" > version.h
	@echo "#define GIT_HASH \"$(GIT_HASH)\"" >> version.h
	@echo "#define GIT_BRANCH \"$(GIT_BRANCH)\"" >> version.h
	@echo "#define COMPILETIME \"`date`\"" >> version.h

clean:
	rm -f sn3d exspec *.o version.h
