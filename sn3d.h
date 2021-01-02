#ifndef SN3D_H
#define SN3D_H

#include "globals.h"

#include <stdarg.h>  /// MK: needed for printout()
#include <cstdbool>
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

#define COOLINGCUT              0.99 //1.01
//#define TAKE_N_BFCONTINUA       900 //40 //900 //20 //900 //40 //20

#define RPKT_EVENTTYPE_BB 550
#define RPKT_EVENTTYPE_CONT 551

#define MAX_RSCAT 50000
#define MIN_XS 1e-40

extern int tid;
extern bool use_cellhist;
extern bool neutral_flag;
extern gsl_rng *rng;  // pointer for random number generator
extern gsl_integration_workspace *gslworkspace;
extern FILE *output_file;

#ifdef _OPENMP
  #pragma omp threadprivate(tid, use_cellhist, neutral_flag, rng, gslworkspace, output_file)
#endif


#define printout(...) fprintf (output_file, __VA_ARGS__)

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

#ifdef DEBUG_ON
  #ifdef assert
    #undef assert
  #endif
  #define assert(e) if (!(e)) { printout("%s:%u: failed assertion `%s' in function %s\n", __FILE__, __LINE__, #e, __PRETTY_FUNCTION__); abort(); }

  #if defined TESTMODE && TESTMODE
    #define assert_testmodeonly(e) if (!(e)) { printout("%s:%u: failed assertion `%s' in function %s\n", __FILE__, __LINE__, #e, __PRETTY_FUNCTION__); abort(); }
  #else
    #define	assert_testmodeonly(e)	((void)0)
  #endif

#else
  #define	assert(e)	((void)0)
  #define	assert_testmodeonly(e)	((void)0)
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
  else
    return file;
}


inline int get_timestep(const double time)
{
  assert(time >= globals::tmin);
  assert(time < globals::tmax);
  for (int nts = 0; nts < globals::ntstep; nts++)
  {
    const double tsend = globals::time_step[nts].start + globals::time_step[nts].width;
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


inline int get_nthreads(void)
{
  #ifdef _OPENMP
  return omp_get_num_threads();
  #else
  return 1;
  #endif
}


#endif // SN3D_H
