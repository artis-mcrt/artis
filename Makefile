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
	CFLAGS = -std=c17 -O3 -fstrict-aliasing -ftree-vectorize -flto -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF

	CFLAGS += -Winline -Wall -Wextra -Wredundant-decls -Wundef -Wstrict-prototypes -Wmissing-prototypes -Wno-unused-parameter -Wno-unused-function -Wstrict-aliasing

	# CFLAGS += -fopenmp-simd
	# CFLAGS += -fvectorize

	# enable OpenMP (for Clang)
	# CFLAGS += -Xpreprocessor -fopenmp -lomp

	# in GCC6, -Wmisleading-indentation will be useful
	# also -fopenmp after -I$(INCLUDE)
	# maybe  -fopt-info-vec-missed
	#  -fwhole-program
	# add -lprofiler for gperftools
#-Wl
	LDFLAGS = $(LIB) -lgsl -lgslcblas -lprofiler
	# sn3d: CFLAGS += -fopenmp

else ifneq (,$(findstring kelvin,$(HOSTNAME)))
	# needs
	#  mpi/openmpi/1.8.5/gcc-4.4.7
	#  compilers/gcc/system(default)
	#  libs/gsl/1.16/gcc-4.4.7

	CC = mpicc
	CFLAGS = -DWALLTIMELIMITSECONDS=\($(WALLTIMEHOURS)\*3600\) -mcmodel=medium -O3 -std=c11 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -I$(GSLINCLUDE) #-fopenmp=libomp
	LDFLAGS= -lgsl -lgslcblas -lm -L$(GSLLIB)

	sn3d: CFLAGS += -DMPI_ON

else ifneq (, $(shell which mpicc))
	# recommended for NCI Raijin cluster:
	# module load intel-cc
	# module load intel-mpi
	# module load gsl

  CC = mpicc
  CFLAGS = -DWALLTIMELIMITSECONDS=\($(WALLTIMEHOURS)\*3600\) -mcmodel=medium -march=native -Wstrict-aliasing -O3 -fstrict-aliasing -std=c11 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF #-fopenmp=libomp
  LDFLAGS= -lgsl -lgslcblas -lm

  ifeq (,$(findstring raijin,$(HOSTNAME)))
    LDFLAGS += -lgslcblas
	endif

  sn3d: CFLAGS += -DMPI_ON
else
	  CFLAGS = -mcmodel=medium -march=native -Wstrict-aliasing -O3 -fstrict-aliasing -std=c11 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF #-fopenmp=libomp
	  LDFLAGS= -lgsl -lgslcblas -lm
endif



### Settings for MPA machines
ifeq ($(DOMAIN),MPA-Garching.MPG.DE)
  #CC = gcc
  CC    = /afs/rzg/@sys/bin/icc
  INCLUDE=/afs/mpa/project/artis/code/lib/gsl32/include/
  LIB=/afs/mpa/project/artis/code/lib/gsl32
  ifeq ($(MACHTYPE),x86_64)
    INCLUDE=/afs/mpa/project/artis/code/lib/gsl64/include/
    LIB=/afs/mpa/project/artis/code/lib/gsl64
  endif
  #CFLAGS =  -Wall -g -I$(INCLUDE)
  #CFLAGS =  -O3 -g -I$(INCLUDE)
  CFLAGS = -O3 -g -pg -I$(INCLUDE)
  #CFLAGS = -openmp  -g -I$(INCLUDE)
  #CFLAGS = -openmp -O3 -g -I$(INCLUDE)
  LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm
  exspec: override CFLAGS =  -g -O3 -I$(INCLUDE) -DDO_EXSPEC
  exgamma: override CFLAGS =  -g -O3 -I$(INCLUDE) -DDO_EXSPEC

  ### desktop MPI
  #CC    =  /opt/mpich2-1.0.3/bin/mpicc #/usr/common/pdsoft/appl/mpich-1.2.6/bin/mpicc #/usr/common/pdsoft/bin/mpicc
  #INCLUDE=/afs/mpa/common/pdsoft/include/
  #LIB=/afs/mpa/common/pdsoft/lib/
  #CFLAGS = -O3 -g -I$(INCLUDE)  -DMPI_ON
  #LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm
  #exspec: override CFLAGS =  -g -O3 -I$(INCLUDE) -DMPI_ON -DDO_EXSPEC
endif



### Settings for IBM Regatta
ifeq ($(OSTYPE),aix)
  #
  #ifeq ($(HOSTTYPE),rs6000)
  CC = mpcc_r
  #INCLUDE=/u/mkromer/lib/regatta64/include
  #LIB=/u/mkromer/lib/regatta64/lib
#  INCLUDE=/u/mkromer/lib/power6_xlc_32_o3/include
#  LIB=/u/mkromer/lib/power6_xlc_32_o3/lib
  #INCLUDE=/u/mkromer/lib/power6_xlcr_64_o3/include
  #LIB=/u/mkromer/lib/power6_xlcr_64_o3/lib
  INCLUDE=/u/mkromer/lib/test_64_r/include
  LIB=/u/mkromer/lib/test_64_r/lib
  #CFLAGS = -O4 -I$(INCLUDE) -qsmp -q64   -qcpluscmt -DMPI_ON
  CFLAGS = -O3 -g -I$(INCLUDE) -q64 -qstrict -qcpluscmt -DMPI_ON #-DPOWER6
  CFLAGS = -O3 -g -I$(INCLUDE) -qsmp -q64 -qstrict -qcpluscmt -DMPI_ON #-DPOWER6
#  CFLAGS = -O3 -g -I$(INCLUDE) -qstrict -qcpluscmt -bmaxdata:0x80000000 -DMPI_ON #-DPOWER6
  #CFLAGS = -O3 -g -I$(INCLUDE) -qsmp -q64 -qstrict -qcpluscmt
  #CFLAGS = -O3 -g -I$(INCLUDE) -q64 -qstrict -qcpluscmt -bmaxdata:0x160000000 -bmaxstack:0x160000000
  LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm
  exspec: override CFLAGS = -O3 -g -q64 -I$(INCLUDE) -qstrict -qcpluscmt -DMPI_ON -DDO_EXSPEC
  exgamma: override CFLAGS = -O3 -g -q64 -I$(INCLUDE) -qstrict -qcpluscmt -DMPI_ON -DDO_EXSPEC
endif



### Settings for the RZG BlueGene
ifeq ($(HOSTNAME),genius1.rzg.mpg.de)
  CC     = mpixlc_r
  INCLUDE= /u/mkromer/lib/BG/include
  LIB    = /u/mkromer/lib/BG/lib
  CFLAGS = -g -O3 -I$(INCLUDE) -qstrict -qcpluscmt -DMPI_ON
  CFLAGS = -g -O3 -I$(INCLUDE) -qsmp -qstrict -qcpluscmt -DMPI_ON
  ##try qsmp=omp the higher optimisation levels up to O5 and no -qstrict, arch
  LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm
endif



### Settings for the Juelich BlueGene
ifeq ($(findstring jugene,$(HOSTNAME)), jugene)
  #this requires a
  #  module load gsl
  #check available module with module avail
  CC     = mpixlc_r
  CFLAGS = -g -O3 -I$(GSL_INCLUDE) -qarch=450 -qtune=450 -qsmp -qstrict -qcpluscmt -DMPI_ON
  ##try qsmp=omp the higher optimisation levels up to O5 and no -qstrict, arch
  LDFLAGS= -L$(GSL_LIB) -lgsl -lgslcblas -lm
endif



### Settings for JUROPA
ifeq ($(WORK),/lustre/jwork1)
  #this requires a
  #  module load gsl
  #check available module with module avail
  CC     = mpicc
  CFLAGS = -O3 -I$(GSL_ROOT) -mcmodel medium -shared-intel -openmp -DMPI_ON
  CFLAGS = -O3 -I$(GSL_ROOT) -openmp -DMPI_ON
  LDFLAGS= -L$(GSL_ROOT) -lgsl -lgslcblas -lm
  exspec exgamma: override CFLAGS = -O3 -I$(GSL_ROOT) -mcmodel medium -shared-intel -DDO_EXSPEC
endif


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
sn3d_files = sn3d.c atomic.c boundary.c emissivities.c gamma.c grey_emissivities.c grid_init.c input.c kpkt.c ltepop.c macroatom.c nltepop.c nonthermal.c decay.c packet_init.c photo_electric.c polarization.c radfield.c ratecoeff.c rpkt.c thermalbalance.c update_grid.c update_packets.c vectors.c vpkt.c md5.c

sn3d_objects = sn3d.o atomic.o boundary.o emissivities.o gamma.o grey_emissivities.o grid_init.o input.o kpkt.o ltepop.o macroatom.o nltepop.o nonthermal.o decay.o packet_init.o photo_electric.o polarization.o radfield.o ratecoeff.o rpkt.o thermalbalance.o update_grid.o update_packets.o vectors.o vpkt.o md5.o

all: sn3d exspec

sn3d: clean version
	$(CC) $(CFLAGS) $(sn3d_files) $(LDFLAGS) -o sn3d

sn3ddebug: clean version $(sn3d_objects)
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) $(sn3d_objects) -o sn3d

exspec_files = exspec.c grid_init.c input.c vectors.c packet_init.c update_grid.c update_packets.c gamma.c boundary.c macroatom.c decay.c rpkt.c kpkt.c photo_electric.c emissivities.c grey_emissivities.c ltepop.c atomic.c ratecoeff.c thermalbalance.c light_curve.c spectrum.c polarization.c nltepop.c radfield.c nonthermal.c vpkt.c md5.c

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
