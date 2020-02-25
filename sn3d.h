#ifndef SN3D_H
#define SN3D_H

#ifndef __CUDA_ARCH__
  #define printout(...) fprintf (output_file, __VA_ARGS__)

  #ifdef DEBUG_ON
    #ifdef assert
      #undef assert
    #endif
    #define assert(e) if (!(e)) { printout("%s:%u: failed assertion `%s' in function %s\n", __FILE__, __LINE__, #e, __PRETTY_FUNCTION__); abort(); }
  #else
    #define	assert(e)	((void)0)
  #endif

  #ifdef _OPENMP
    #define safeadd(var, val) #pragma omp atomic\\
    var += val
  #else
    #define safeadd(var, val) var += val
  #endif

#else

  #define printout(...) printf (__VA_ARGS__)

  #define safeadd(var, val) atomicAdd(&var, val)

#endif

#define safeincrement(var) safeadd(var, 1)

#include "cuda.h"

#include "globals.h"

#include <unistd.h>
#include <stdarg.h>  /// MK: needed for printout()
#include <stdbool.h>
#ifndef __CUDA_ARCH__
#include <gsl/gsl_integration.h>
#endif
#include <iostream>

#define DEBUG_ON
// #define DO_TITER
// #define FORCE_LTE

#include "types.h"

#if (DETAILED_BF_ESTIMATORS_ON && !NO_LUT_PHOTOION)
  #error Must use NO_LUT_PHOTOION with DETAILED_BF_ESTIMATORS_ON
#endif

#if !defined DO_EXSPEC && !defined MPI_ON
  // #define MPI_ON //only needed for debugging MPI, the makefile will switch this on
#endif

#ifdef MPI_ON
  #include "mpi.h"
#endif

//#define _OPENMP
#ifdef _OPENMP
  #include "omp.h"
#endif

#ifndef _OPENMP
typedef int omp_int_t;
__host__ __device__ static inline omp_int_t omp_get_thread_num(void) {
#ifdef __CUDA_ARCH__
  return threadIdx.x + blockDim.x * blockIdx.x;
#else
  return 0;
#endif
}
// static inline omp_int_t omp_get_num_threads(void) { return 1; }
__host__ __device__ static inline omp_int_t omp_get_num_threads(void) { return MTHREADS; }
#endif

#define GRID_UNIFORM 1 // Simple cuboidal cells.
#define GRID_SPHERICAL1D 2 // radial shells

#define COOLING_UNDEFINED       -99

#define COOLINGCUT              0.99 //1.01
//#define TAKE_N_BFCONTINUA       900 //40 //900 //20 //900 //40 //20

#define RPKT_EVENTTYPE_BB 550
#define RPKT_EVENTTYPE_CONT 551

#define PACKET_SAME -929 //MUST be negative

#define MAX_RSCAT 50000
#define MIN_XS 1e-40



// number of ion stats counters that should be divided by the ion populations
#define nstatcounters_ratecoeff 18

void increment_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstatscounters ion_counter_type, const double increment);

double get_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstatscounters ion_counter_type);

void set_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstatscounters ion_counter_type, const double newvalue);

// extern int tid;
#ifndef __CUDA_ARCH__
extern gsl_rng *rng;  // pointer for random number generator
extern gsl_integration_workspace *gslworkspace;
#else
extern __device__ void *rng;
extern __device__ void *gslworkspace;
#endif

extern __managed__ int myGpuId;
extern __managed__ bool use_cellhist;
extern __managed__ bool neutral_flag;
extern FILE *output_file;

#ifdef _OPENMP
  #pragma omp threadprivate(tid, myGpuId, use_cellhist, neutral_flag, rng, gslworkspace, output_file)
#endif



// inline int printout(const char *format, ...)
// {
//    va_list args;
//    va_start(args, format);
//    const int ret_status = vfprintf(output_file, format, args);
//    // fprintf(output_file, "vfprintf return code %d\n", va_arg(args, char*) == NULL);
//    va_end(args);
//
//    return ret_status;
// }


inline void gsl_error_handler_printout(const char *reason, const char *file, int line, int gsl_errno)
{
  if (gsl_errno != 18) // roundoff error
  {
    printout("WARNING: gsl (%s:%d): %s (Error code %d)\n", file, line, reason, gsl_errno);
    // abort();
  }
}


inline FILE *fopen_required(const char *filename, const char *mode)
{
  FILE *file = fopen(filename, mode);
  if (file == NULL)
  {
    printout("ERROR: Could not open file '%s' for mode '%s'.\n", filename, mode);
    abort();
  }

  return file;
}

#if CUDA_ENABLED
void* makemanaged(void* ptr, size_t curSize);
#endif

#endif // SN3D_H
