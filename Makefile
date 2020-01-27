WALLTIMEHOURS := 24
# WALLTIMEHOURS IS DEPRECATED: use command line argument e.g. "./sn3d -w 24" in the future
GIT_VERSION := $(shell git describe --dirty --always --tags)
GIT_HASH := $(shell git rev-parse HEAD)
GIT_BRANCH := $(shell git branch | sed -n '/\* /s///p')

SYSNAME := $(shell uname -s)

ifeq ($(SYSNAME),Darwin)
	# macOS

	#if using homebrew gsl, gperftools, etc, you might need to add:
	# export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/include/
	# export LIBRARY_PATH=$LIBRARY_PATH:/usr/local/lib/
	# to your .zshrc/.bashrc startup script

	CC = clang
	# CC = gcc-9
	# CC = icc
	# CC = mpicc
	CFLAGS = -std=c17 -O3 -fstrict-aliasing -ftree-vectorize -flto

	CFLAGS += -Winline -Wall -Wextra -Wredundant-decls -Wundef -Wstrict-prototypes -Wmissing-prototypes -Wno-unused-parameter -Wno-unused-function -Wstrict-aliasing

	# CFLAGS += -fopenmp-simd
	# CFLAGS += -fvectorize

	# enable OpenMP for Clang
	# CFLAGS += -Xpreprocessor -fopenmp -lomp

	# in GCC6, -Wmisleading-indentation will be useful
	# also -fopenmp after -I$(INCLUDE)
	# maybe  -fopt-info-vec-missed
	#  -fwhole-program
	# add -lprofiler for gperftools
	LDFLAGS = $(LIB) -lgsl -lgslcblas -lprofiler

else ifneq (,$(findstring kelvin,$(HOSTNAME)))
  # QUB Kelvin cluster
	# needs
	#  mpi/openmpi/1.8.5/gcc-4.4.7
	#  compilers/gcc/system(default)
	#  libs/gsl/1.16/gcc-4.4.7

	CC = mpicc
	CFLAGS = -DWALLTIMELIMITSECONDS=\($(WALLTIMEHOURS)\*3600\) -mcmodel=medium -O3 -std=c11 -I$(GSLINCLUDE) #-fopenmp=libomp
	LDFLAGS= -lgsl -lgslcblas -lm -L$(GSLLIB)

	sn3d: CFLAGS += -DMPI_ON

else ifneq (, $(shell which mpicc))
  # any other system which has mpicc available (Juwels, Cambridge, Gadi, etc)

  CC = mpicc
  CFLAGS = -DWALLTIMELIMITSECONDS=\($(WALLTIMEHOURS)\*3600\) -mcmodel=medium -march=native -Wstrict-aliasing -O3 -fstrict-aliasing -std=c11 #-fopenmp=libomp
  LDFLAGS= -lgsl -lgslcblas -lm

  ifeq (,$(findstring raijin,$(HOSTNAME)))
    LDFLAGS += -lgslcblas
	endif

  sn3d: CFLAGS += -DMPI_ON
else
	  CFLAGS = -mcmodel=medium -march=native -Wstrict-aliasing -O3 -fstrict-aliasing -std=c11 #-fopenmp=libomp
	  LDFLAGS= -lgsl -lgslcblas -lm
endif


# Use GSL inline functions and skip array range checking for performance
CFLAGS += -DHAVE_INLINE -DGSL_C99_INLINE -DGSL_RANGE_CHECK_OFF

exspec exgamma: CFLAGS += -DDO_EXSPEC
exgamma: CFLAGS += -DDO_EXGAMMA

sn3dmpi: CC = mpicc
sn3dmpi: CFLAGS += -DMPI_ON
sn3dmpi: sn3d

sn3dopenmp: CFLAGS += -Xpreprocessor
sn3dopenmp: CFLAGS += -fopenmp
sn3dopenmp: LDFLAGS += -lomp
sn3dopenmp: sn3d


### use pg when you want to use gprof the profiler
#CFLAGS = -g -pg -Wall -I$(INCLUDE)
sn3d_files = sn3d.c atomic.c boundary.c emissivities.c gamma.c globals.c grey_emissivities.c grid_init.c input.c kpkt.c ltepop.c macroatom.c nltepop.c nonthermal.c decay.c packet_init.c photo_electric.c polarization.c radfield.c ratecoeff.c rpkt.c thermalbalance.c update_grid.c update_packets.c vectors.c vpkt.c md5.c

sn3d_objects = sn3d.o atomic.o boundary.o emissivities.o gamma.o globals.o grey_emissivities.o grid_init.o input.o kpkt.o ltepop.o macroatom.o nltepop.o nonthermal.o decay.o packet_init.o photo_electric.o polarization.o radfield.o ratecoeff.o rpkt.o thermalbalance.o update_grid.o update_packets.o vectors.o vpkt.o md5.o

all: sn3d exspec

sn3d: clean version
	$(CC) $(CFLAGS) $(sn3d_files) $(LDFLAGS) -o sn3d

sn3ddebug: clean version $(sn3d_objects)
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) $(sn3d_objects) -o sn3d

exspec_files = exspec.c grid_init.c globals.c input.c vectors.c packet_init.c update_grid.c update_packets.c gamma.c boundary.c macroatom.c decay.c rpkt.c kpkt.c photo_electric.c emissivities.c grey_emissivities.c ltepop.c atomic.c ratecoeff.c thermalbalance.c light_curve.c spectrum.c polarization.c nltepop.c radfield.c nonthermal.c vpkt.c md5.c

exspec: clean version
	$(CC) $(CFLAGS) $(exspec_files) $(LDFLAGS) -o exspec

exgamma: clean version
	$(CC) $(CFLAGS) $(exspec_files) $(LDFLAGS) -o exgamma


.PHONY: clean version

version:
	@echo "#define GIT_VERSION \"$(GIT_VERSION)\"" > version.h
	@echo "#define GIT_HASH \"$(GIT_HASH)\"" >> version.h
	@echo "#define GIT_BRANCH \"$(GIT_BRANCH)\"" >> version.h
	@echo "#define COMPILETIME \"`date`\"" >> version.h

clean:
	rm -f *.o version.h
