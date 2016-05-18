SHELL = /bin/sh

INCLUDE=/usr/local/opt/gsl/include
LIB=/usr/local/opt/gsl/lib


### Settings for MACBOOK
#ifeq ($(USER),mattia)

#needs
#module load intel-cc/11.1.046
#module load openmpi/1.4.3
#module load gsl/1.12                                       

#  CC = mpicc
#  CFLAGS = -O3 -march=native -DMPI_ON -I$(INCLUDE)
#  CC = cc
#  CFLAGS = -O3 -march=native -flto -I$(INCLUDE)
  CC = clang-omp
  CFLAGS = -O3 -march=native -flto -I$(INCLUDE) -fopenmp

  LDFLAGS=  -L$(LIB) -lgsl -lgslcblas -lm

  exspec: override CFLAGS =  -O3  -DDO_EXSPEC -I$(INCLUDE)
  exspec_dd: override CFLAGS =  -O3  -DDO_EXSPEC -I$(INCLUDE)
  exgamma: override CFLAGS =  -O3  -DDO_EXSPEC -I$(INCLUDE)


#endif

### Settings for STARBASE                                                                                          
ifeq ($(USER),mb)

#needs                                                                                                           
#module load intel_comp/c4                                                                                       
#module load platform_mpi/8.2.1                                                                                
#module load gsl    
                                                                                                
  CC = mpicc
  CFLAGS = -O3 -DMPI_ON -I$(INCLUDE)
  LDFLAGS= -lgsl -lgslcblas -lm -L$(LIB)

  exspec: override CFLAGS =  -O3  -DDO_EXSPEC -I$(INCLUDE)
  exspec_dd: override CFLAGS =  -O3  -DDO_EXSPEC
  exgamma: override CFLAGS =  -O3  -DDO_EXSPEC

endif


### Settings for cosma                                                                                            
ifeq ($(HOSTNAME),login)

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


### Settings for leicester                                                                                          
                                                                                                                  
ifeq ($(USER),dc-bull1)

#needs                                                                                                           
#module load intel/compilers/13.0.0
#module load intel/impi/4.1.3
#module load gsl/intel/1.15                                                                                         
                                                                                                                  
  CC = mpicc
  CFLAGS = -O3 -DMPI_ON
  LDFLAGS= -lgsl -lgslcblas -lm

  exspec: override CFLAGS =  -O3  -DDO_EXSPEC
  exspec_dd: override CFLAGS =  -O3  -DDO_EXSPEC
  exgamma: override CFLAGS =  -O3  -DDO_EXSPEC

endif


### Settings for the Juelich BlueGene/Q                                                                           
ifeq ($(findstring juqueen,$(HOSTNAME)), juqueen)
 #this requires a                                                                                                 
 #  module load gsl                                                                                               
 #check available module with module avail                                                                        
 CC     = mpixlc_r
 #CFLAGS = -g -O3 -I$(GSL_INCLUDE) -qarch=qp -qtune=qp -qsmp=omp -qthreaded -qstrict -qcpluscmt -DMPI_ON          
 CFLAGS = -O3 -I$(GSL_INCLUDE) -qarch=qp -qtune=qp -qinline -qsmp=omp -qthreaded -qcpluscmt -DMPI_ON
 #CFLAGS =  -O4 -I$(GSL_INCLUDE) -qarch=qp -qtune=qp -qnoipa -qinline -qsmp=omp -qthreaded -qcpluscmt -DMPI_ON    
 LDFLAGS= -L$(GSL_LIB) -lgsl -lgslcblas -lm -qthreaded
endif



### use pg when you want to use gprof the profiler
#CFLAGS = -g -pg -Wall -I$(INCLUDE)


# Get the compiletime and the git hash of the current version
# and write it to version.h to so that the code knows which
# revision was used for a particular run



sn3d_objects = sn3d.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o packet_prop.o compton.o macroatom.o rpkt.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o grey_emissivities.o syn_lc.o  ltepop.o atomic.o ratecoeff.o thermalbalance.o polarization.o vpkt.o

sn3d: version $(sn3d_objects) 
	$(CC) $(CFLAGS) $(sn3d_objects) $(LDFLAGS) -o sn3d

exspec_objects = exspec.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o packet_prop.o compton.o macroatom.o rpkt.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o grey_emissivities.o syn_lc.o  ltepop.o atomic.o ratecoeff.o thermalbalance.o light_curve.o gamma_light_curve.o spectrum.o polarization.o specpol.o vpkt.o

exspec: version $(exspec_objects) 
	$(CC) $(CFLAGS) $(exspec_objects) $(LDFLAGS) -o exspec

exspec_dd_objects = exspec_dd.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o packet_prop.o compton.o macroatom.o rpkt.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o grey_emissivities.o syn_lc.o  ltepop.o atomic.o ratecoeff.o thermalbalance.o light_curve.o gamma_light_curve.o spectrum.o polarization.o specpol.o vpkt.o

exspec_dd: version $(exspec_dd_objects) 
	$(CC) $(CFLAGS) $(exspec_dd_objects) $(LDFLAGS) -o exspec_dd

exgamma_objects = exgamma.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o packet_prop.o compton.o macroatom.o rpkt.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o grey_emissivities.o syn_lc.o  ltepop.o atomic.o ratecoeff.o thermalbalance.o light_curve.o gamma_light_curve.o spectrum.o polarization.o specpol.o vpkt.o

exgamma: version $(exgamma_objects) 
	$(CC) $(CFLAGS) $(exgamma_objects) $(LDFLAGS) -o exgamma

sn3dsyn_objects = sn3dsyn.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o spectrum.o packet_prop.o compton.o rpkt.o light_curve.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o gamma_light_curve.o grey_emissivities.o syn_lc.o light_curve_res.o polarization.o specpol.o vpkt.o

sn3dsyn: version $(sn3dsyn_objects) 
	$(CC) $(CFLAGS) $(sn3dsyn_objects) $(LDFLAGS) -o sn3dsyn

sn3dlcsyn_objects = sn3dlcsyn.o grid_init.o input.o vectors.o packet_init.o time_init.o update_grid.o update_packets.o gamma.o boundary.o move.o spectrum.o packet_prop.o compton.o rpkt.o light_curve.o kpkt.o photo_electric.o linelist.o syn_gamma.o ray_prop.o update_gamma_rays.o emissivities.o gamma_light_curve.o grey_emissivities.o syn_lc.o light_curve_res.o polarization.o specpol.o vpkt.o

sn3dlcsyn: version $(sn3dlcsyn_objects) 
	$(CC) $(CFLAGS) $(sn3dlcsyn_objects) $(LDFLAGS) -o sn3dlcsyn

version:
	@echo "#define GIT_HASH \"`cat .git/refs/heads/polarization`\"" > version.h
	@echo "#define COMPILETIME \"`date`\"" >> version.h




clean:
	rm -f *o

veryclean:
	rm -f *o *exe *~ 







