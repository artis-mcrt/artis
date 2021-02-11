#ifndef SN3D_H
#define SN3D_H

#ifndef __CUDA_ARCH__
  // host code

  #define printout(...) fprintf (output_file, __VA_ARGS__)

  #ifdef DEBUG_ON
    #ifdef assert
      #undef assert
    #endif
    #define assert(e) if (!(e)) { printout("%s:%u: failed assertion `%s' in function %s\n", __FILE__, __LINE__, #e, __PRETTY_FUNCTION__); abort(); }

    #if defined TESTMODE && TESTMODE
      #define assert_testmodeonly(e) if (!(e)) { printout("%s:%u: failed testmode assertion `%s' in function %s\n", __FILE__, __LINE__, #e, __PRETTY_FUNCTION__); abort(); }
    #else
      #define	assert_testmodeonly(e)	((void)0)
    #endif

  #else
    #ifdef assert
      #undef assert
    #endif
    #define	assert(e) ((void)0)
    #define	assert_testmodeonly(e) ((void)0)
  #endif

  #ifdef _OPENMP
  #ifndef __CUDACC__
    #define safeadd(var, val) _Pragma("omp atomic update") \
    var += val
  #else
    #define safeadd(var, val) var += val
  #endif
  #else
    #define safeadd(var, val) var += val
  #endif

#else
  // device code

  #define printout(...) printf (__VA_ARGS__)

  #if defined TESTMODE && TESTMODE
    #define assert_testmodeonly(e) assert(e)
  #else
    #define	assert_testmodeonly(e) ((void)0)
  #endif


  #define safeadd(var, val) atomicAdd(&var, val)

#endif

#define safeincrement(var) safeadd(var, 1)

#include "cuda.h"

#include "globals.h"

#include <stdarg.h>  /// MK: needed for printout()
#include <gsl/gsl_integration.h>

#define DEBUG_ON
// #define DO_TITER
// #define FORCE_LTE

#include "types.h"
#include "vectors.h"

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

#define GRID_UNIFORM 1 // Simple cuboidal cells.
#define GRID_SPHERICAL1D 2 // radial shells

#define COOLING_UNDEFINED       -99

#define RPKT_EVENTTYPE_BB 550
#define RPKT_EVENTTYPE_CONT 551

extern int tid;
extern __managed__ bool use_cellhist;
extern __managed__ bool neutral_flag;
#ifndef __CUDA_ARCH__
extern gsl_rng *rng;  // pointer for random number generator
#else
extern __device__ void *rng;
#endif
extern gsl_integration_workspace *gslworkspace;
extern FILE *output_file;
extern __managed__ int myGpuId;

#ifdef _OPENMP
  #pragma omp threadprivate(tid, myGpuId, use_cellhist, neutral_flag, rng, gslworkspace, output_file)
#endif


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


inline int get_timestep(const double time)
{
  assert(time >= globals::tmin);
  assert(time < globals::tmax);
  for (int nts = 0; nts < globals::ntstep; nts++)
  {
    const double tsend = (nts < (globals::ntstep - 1)) ? globals::time_step[nts + 1].start : globals::tmax;
    if (time >= globals::time_step[nts].start && time < tsend)
    {
      return nts;
    }
  }
  assert(false); // could not find matching timestep

  return -1;
}


inline double get_arrive_time(const PKT *pkt_ptr)
/// We know that a packet escaped at "escape_time". However, we have
/// to allow for travel time. Use the formula in Leon's paper. The extra
/// distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
{
  return pkt_ptr->escape_time - (dot(pkt_ptr->pos, pkt_ptr->dir) / globals::CLIGHT_PROP);
}


inline double get_arrive_time_cmf(const PKT *pkt_ptr)
{
  return pkt_ptr->escape_time * sqrt(1. - (globals::vmax * globals::vmax / CLIGHTSQUARED));
}


__host__ __device__
inline int get_max_threads(void)
{
#ifdef __CUDA_ARCH__
  return MCUDATHREADS;
#elif defined _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}


__host__ __device__
inline int get_num_threads(void)
{
#ifdef __CUDA_ARCH__
  return blockDim.x * blockDim.y * blockDim.z;
#elif defined _OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}


__host__ __device__
inline int get_thread_num(void)
{
#ifdef __CUDA_ARCH__
  return threadIdx.x + blockDim.x * blockIdx.x;
#elif defined _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}


#endif // SN3D_H
