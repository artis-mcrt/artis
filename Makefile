SHELL = /bin/sh
GIT_HASH = $(shell git rev-parse HEAD)
GIT_BRANCH = $(shell git branch | sed -n '/\* /s///p')

### Settings for the miner
ifeq ($(SYSTYPE),Haswell)
  CC = mpicc
  CFLAGS = -O3 -DMPI_ON
  LDFLAGS= -lm -L /hits/sw/tap/libs/gsl-2.1/lib -lgsl -lgslcblas
  exspec: override CFLAGS = -O3 -DDO_EXSPEC
  exspec_dd: override CFLAGS = -O3 -DDO_EXSPEC
  exgamma: override CFLAGS = -O3 -DDO_EXSPEC
endif


### Settings for the miner
ifeq ($(OSTYPE),linux)
  CC = cc
  INCLUDE = /home/ssim/gsl/include
  LIB = /home/ssim/gsl/lib
  CFLAGS = -O3 -g -I$(INCLUDE)
  LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm
endif


### Settings for Coala
ifeq ($(OSTYPE),linux)
  CC = /pkg/linux/SS12/sunstudio12/bin/cc
#  CC = gcc
  INCLUDE = /home/ssim/gsl/include
  INCLUDE2 = /usr/local/openmpi/include
  LIB2 = /usr/local/openmpi/lib
  LIB = /home/ssim/gsl/lib
#  CFLAGS = -O3 -I$(INCLUDE) -I$(INCLUDE2) -fast -xtarget=nehalem -xipo=2 -xvector=simd -DMPI_ON
  CFLAGS = -O3 -I$(INCLUDE) -I$(INCLUDE2)  -pthread -DMPI_ON
#-fast -xtarget=nehalem -xipo=2 -xvector=simd -DMPI_ON
  LDFLAGS= -L$(LIB) -L$(LIB2) -R$(LIB2) -lgsl -lgslcblas -lm -pthread -L/usr/local/openmpi/lib -lmpi_cxx -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl
#  LDFLAGS= -L$(LIB) -L$(LIB2) -L$(LIB3) -lgsl -lgslcblas -lm -pthread
endif


### Settings for mime
ifeq ($(HOST),mime)
  CC = mpicc
  CFLAGS = -O3 -m64 -DMPI_ON -DTIMED_RESTARTS
  LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm -m64  
#-pthread -lmpi_cxx -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl
endif



### Settings for vayu
ifeq ($(USER),mbu552)

#needs 
#module load intel-cc/11.1.046
#module load openmpi/1.4.3
#module load gsl/1.15

  CC = mpicc
  CFLAGS = -O3 -DMPI_ON -DTIMED_RESTARTS
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
  exspec_dd: override CFLAGS =  -g -O3 -I$(INCLUDE) -DDO_EXSPEC
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
  CFLAGS = -O3 -g -I$(INCLUDE) -q64 -qsmp=omp -qthreaded -qstrict -qcpluscmt -DMPI_ON #-DPOWER6
#  CFLAGS = -O3 -g -I$(INCLUDE) -qstrict -qcpluscmt -bmaxdata:0x80000000 -DMPI_ON #-DPOWER6
  #CFLAGS = -O3 -g -I$(INCLUDE) -qsmp -q64 -qstrict -qcpluscmt 
  #CFLAGS = -O3 -g -I$(INCLUDE) -q64 -qstrict -qcpluscmt -bmaxdata:0x160000000 -bmaxstack:0x160000000
  LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm -qthreaded
  exspec: override CFLAGS = -O3 -g -q64 -I$(INCLUDE) -qstrict -qcpluscmt  -DDO_EXSPEC
  exspec_dd: override CFLAGS = -O3 -g -q64 -I$(INCLUDE) -qstrict -qcpluscmt  -DDO_EXSPEC
  exgamma: override CFLAGS = -O3 -g -q64 -I$(INCLUDE) -qstrict -qcpluscmt  -DDO_EXSPEC
endif



### Settings for the OPA cluster
ifeq ($(DOMAIN),opt.rzg.mpg.de)
  CC    = mpiicc
  #this requires a 
  #  module load intel
  #  module load impi
  #  module load gsl
  #check available module with module avail
  #Read gsl now from system wide installation
  #INCLUDE=/afs/ipp-garching.mpg.de/home/m/mkromer/lib/opa/include/
  #LIB=/afs/ipp-garching.mpg.de/home/m/mkromer/lib/opa/lib
  #INCLUDE=/afs/ipp-garching.mpg.de/home/s/ssim/gsl-opa/include/
  #LIB=/afs/ipp-garching.mpg.de/home/s/ssim/gsl-opa/lib
  #CFLAGS = -O2 -openmp -I$(INCLUDE)
  #CFLAGS = -m64 -O2 -mcmodel medium -shared-intel -I$(INCLUDE) -DMPI_ON
  CFLAGS = -m64 -O2 -mcmodel medium -shared-intel $(GSL_CFLAGS) -DMPI_ON
  #"-mcmodel medium -shared-intel" are needed to hold > 2GB static data
  #in memory http://software.intel.com/en-us/forums/showthread.php?t=43717#18089
  #LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm
  LDFLAGS= $(GSL_LDFLAGS) -lgsl -lgslcblas -lm
  exspec: override CFLAGS = -m64 -O2 -mcmodel medium -shared-intel $(GSL_CFLAGS) -DDO_EXSPEC
  exspec_dd: override CFLAGS = -m64 -O2 -mcmodel medium -shared-intel $(GSL_CFLAGS) -DDO_EXSPEC
  exgamma: override CFLAGS = -m64 -O2 -mcmodel medium -shared-intel $(GSL_CFLAGS) -DDO_EXSPEC
endif



### Settings for the RZG BlueGene
ifeq ($(HOSTNAME),genius1.rzg.mpg.de)
  CC     = mpixlc_r
  INCLUDE= /u/mkromer/lib/BG/include
  LIB    = /u/mkromer/lib/BG/lib
  #CFLAGS = -g -O3 -I$(GSL_INCLUDE) -qarch=450 -qtune=450 -qsmp=omp -qthreaded -qstrict -qcpluscmt -DMPI_ON
  CFLAGS = -O5 -I$(GSL_INCLUDE) -qarch=450d -qtune=450 -qsmp=omp -qthreaded -qcpluscmt -DMPI_ON
  ##try qsmp=omp the higher optimisation levels up to O5 and no -qstrict, arch
  LDFLAGS= -L$(LIB) -lgsl -lgslcblas -lm -qthreaded
endif



### Settings for the Juelich BlueGene/P
ifeq ($(findstring jugene,$(HOSTNAME)), jugene)
  #this requires a
  #  module load gsl
  #check available module with module avail
  CC     = mpixlc_r
  #CFLAGS = -g -O3 -I$(GSL_INCLUDE) -qarch=450 -qtune=450 -qsmp=omp -qthreaded -qstrict -qcpluscmt -DMPI_ON
  CFLAGS = -O5 -I$(GSL_INCLUDE) -qarch=450d -qtune=450 -qsmp=omp -qthreaded -qcpluscmt -DMPI_ON
  LDFLAGS= -L$(GSL_LIB) -lgsl -lgslcblas -lm -qthreaded
endif


### Settings for the Juelich BlueGene/Q                                                                           
ifeq ($(findstring juqueen,$(HOSTNAME)), juqueen)
 #this requires a                                                                                                 
 #  module load gsl                                                                                               
 #check available module with module avail                                                                        
 #CC     = mpixlc_r
 CC     = mpixlc  
#CFLAGS = -g -O3 -I$(GSL_INCLUDE) -qarch=qp -qtune=qp -qsmp=omp -qthreaded -qstrict -qcpluscmt -DMPI_ON          
 CFLAGS = -O3 -I$(GSL_INCLUDE) -qarch=qp -qtune=qp -qinline -qsmp=omp -qthreaded -qcpluscmt -DMPI_ON
# CFLAGS = -O3 -I$(GSL_INCLUDE) -qarch=qp -qtune=qp -qinline -qcpluscmt -DMPI_ON 
#CFLAGS =  -O4 -I$(GSL_INCLUDE) -qarch=qp -qtune=qp -qnoipa -qinline -qsmp=omp -qthreaded -qcpluscmt -DMPI_ON    
 LDFLAGS= -L$(GSL_LIB) -lgsl -lgslcblas -lm -qthreaded

  exspec: override CFLAGS =  -O3  -DDO_EXSPEC -I$(GSL_INCLUDE) -qarch=qp -qtune=qp -qinline -qcpluscmt
  exspec_dd: override CFLAGS =  -O3  -DDO_EXSPEC
  exgamma: override CFLAGS =  -O3  -DDO_EXSPEC

endif


### Settings for JUROPA
ifeq ($(WORK),/lustre/jwork/hmu14/hmu146)
  #this requires a
  #  module load gsl
  #check available module with module avail
  CC     = mpicc
  CFLAGS = -O3 -I$(GSL_ROOT)/include -openmp -DMPI_ON
  LDFLAGS= -L$(GSL_ROOT)/lib -lgsl -lgslcblas -lm
  exspec: override CFLAGS = -O3 -I$(GSL_ROOT)/include -DDO_EXSPEC
  exspec_dd: override CFLAGS = -O3 -I$(GSL_ROOT)/include -DDO_EXSPEC
  exgamma: override CFLAGS = -O3 -I$(GSL_ROOT)/include -DDO_EXSPEC
endif


### Settings for MB MACBOOK
ifeq ($(USER),mattia)

  CC = mpicc
  CFLAGS = -O3 -DMPI_ON
  LDFLAGS= -lgsl -lgslcblas -lm

  exspec: override CFLAGS =  -O3  -DDO_EXSPEC
  exspec_dd: override CFLAGS =  -O3  -DDO_EXSPEC
  exgamma: override CFLAGS =  -O3  -DDO_EXSPEC

endif


### Settings for STARBASE                                                                                          
ifeq ($(USER),mb)

#needs                                                                                                           
#module load intel_comp/c4                                                                                       
#module load platform_mpi/8.2.1                                                                                
#module load gsl    
         
  INCLUDE=/home/mb/gsl/include/
  LIB=/home/mb/gsl/lib/
                                                                                       
  CC = mpicc
  CFLAGS = -O3 -DMPI_ON -I$(INCLUDE)
  LDFLAGS= -lgsl -lgslcblas -lm -L$(LIB)

  exspec: override CFLAGS =  -O3  -DDO_EXSPEC -I$(INCLUDE)
  exspec_dd: override CFLAGS =  -O3  -DDO_EXSPEC
  exgamma: override CFLAGS =  -O3  -DDO_EXSPEC

endif

### Settings for cosma                                                                                            
ifeq ($(USER),dc-bull1)

#needs                                                                                                           
#module load intel_comp/c4/2013.0.028                                                                                       
#module load platform_mpi/8.2.1                                                                                
#module load gsl/1.15                                                                                                    
  CC = mpicc
  CFLAGS = -O3 -DMPI_ON
  LDFLAGS= -lgsl -lgslcblas -lm

  exspec: override CFLAGS =  -O3  -DDO_EXSPEC
  exspec_dd: override CFLAGS =  -O3  -DDO_EXSPEC
  exgamma: override CFLAGS =  -O3  -DDO_EXSPEC

endif



### use pg when you want to use gprof the profiler
#CFLAGS = -g -pg -Wall -I$(INCLUDE)


# Get the compiletime and the git hash of the current version
# and write it to version.h to so that the code knows which
# revision was used for a particular run



sn3d_objects = sn3d.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o packet_prop.o compton.o macroatom.o rpkt.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o grey_emissivities.o syn_lc.o  ltepop.o atomic.o ratecoeff.o thermalbalance.o polarization.o vpkt.o energy_init.o

sn3d: version $(sn3d_objects) 
	$(CC) $(CFLAGS) $(sn3d_objects) $(LDFLAGS) -o sn3d.exe

exspec_objects = exspec.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o packet_prop.o compton.o macroatom.o rpkt.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o grey_emissivities.o syn_lc.o  ltepop.o atomic.o ratecoeff.o thermalbalance.o light_curve.o gamma_light_curve.o spectrum.o polarization.o specpol.o vpkt.o energy_init.o

exspec: version $(exspec_objects) 
	$(CC) $(CFLAGS) $(exspec_objects) $(LDFLAGS) -o exspec.exe

exspec_dd_objects = exspec_dd.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o packet_prop.o compton.o macroatom.o rpkt.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o grey_emissivities.o syn_lc.o  ltepop.o atomic.o ratecoeff.o thermalbalance.o light_curve.o gamma_light_curve.o spectrum.o polarization.o specpol.o vpkt.o

exspec_dd: version $(exspec_dd_objects) 
	$(CC) $(CFLAGS) $(exspec_dd_objects) $(LDFLAGS) -o exspec_dd.exe

exgamma_objects = exgamma.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o packet_prop.o compton.o macroatom.o rpkt.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o grey_emissivities.o syn_lc.o  ltepop.o atomic.o ratecoeff.o thermalbalance.o light_curve.o gamma_light_curve.o spectrum.o polarization.o specpol.o vpkt.o energy_init.o

exgamma: version $(exgamma_objects) 
	$(CC) $(CFLAGS) $(exgamma_objects) $(LDFLAGS) -o exgamma.exe

sn3dsyn_objects = sn3dsyn.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o spectrum.o packet_prop.o compton.o rpkt.o light_curve.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o gamma_light_curve.o grey_emissivities.o syn_lc.o light_curve_res.o polarization.o specpol.o vpkt.o

sn3dsyn: version $(sn3dsyn_objects) 
	$(CC) $(CFLAGS) $(sn3dsyn_objects) $(LDFLAGS) -o sn3dsyn.exe

sn3dlcsyn_objects = sn3dlcsyn.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o spectrum.o packet_prop.o compton.o rpkt.o light_curve.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o gamma_light_curve.o grey_emissivities.o syn_lc.o light_curve_res.o polarization.o specpol.o vpkt.o

sn3dlcsyn: version $(sn3dlcsyn_objects) 
	$(CC) $(CFLAGS) $(sn3dlcsyn_objects) $(LDFLAGS) -o sn3dlcsyn.exe

version:
	@echo "#define GIT_HASH \"$(GIT_BRANCH)$(GIT_HASH)\"" > version.h
	@echo "#define COMPILETIME \"`date`\"" >> version.h




clean:
	rm -f *o

veryclean:
	rm -f *o *exe *~ 
