#ifndef SN3D_H
#define SN3D_H

#include <getopt.h>
#include <gsl/gsl_integration.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cassert>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>

#include "artisoptions.h"
#include "globals.h"

// #define _OPENMP
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef MPI_ON
#include <mpi.h>
#endif

extern FILE *output_file;
extern int tid;
extern bool use_cellhist;

extern std::mt19937 stdrng;

extern gsl_integration_workspace *gslworkspace;

#ifdef _OPENMP
#pragma omp threadprivate(tid, use_cellhist, stdrng, gslworkspace, output_file)
#endif

#define __artis_assert(e)                                                                                              \
  {                                                                                                                    \
    const bool pass = static_cast<bool>(e);                                                                            \
    if (!pass) {                                                                                                       \
      if (output_file != nullptr) {                                                                                    \
        (void)fprintf(output_file, "[rank %d] %s:%d: failed assertion `%s' in function %s\n", globals::rank_global,    \
                      __FILE__, __LINE__, #e, __PRETTY_FUNCTION__);                                                    \
      }                                                                                                                \
      (void)fprintf(stderr, "[rank %d] %s:%d: failed assertion `%s' in function %s\n", globals::rank_global, __FILE__, \
                    __LINE__, #e, __PRETTY_FUNCTION__);                                                                \
      abort();                                                                                                         \
    }                                                                                                                  \
    assert(pass);                                                                                                      \
  }

#define assert_always(e) __artis_assert(e)

#ifndef TESTMODE
#define TESTMODE false
#endif

#if defined TESTMODE && TESTMODE
#define assert_testmodeonly(e) __artis_assert(e)
#else
#define assert_testmodeonly(e) \
  if (!(e)) {                  \
    __builtin_unreachable();   \
  }
#endif

// #define printout(...) fprintf(output_file, __VA_ARGS__)

extern int tid;

template <typename... Args>
static int printout(const char *format, Args... args) {
  if (globals::startofline[tid]) {
    const time_t now_time = time(nullptr);
    char s[32] = "";
    strftime(s, 32, "%FT%TZ", gmtime(&now_time));
    fprintf(output_file, "%s ", s);
  }
  globals::startofline[tid] = (format[strlen(format) - 1] == '\n');
  return fprintf(output_file, format, args...);
}

static int printout(const char *format) {
  if (globals::startofline[tid]) {
    const time_t now_time = time(nullptr);
    char s[32] = "";
    strftime(s, 32, "%FT%TZ", gmtime(&now_time));
    fprintf(output_file, "%s ", s);
  }
  globals::startofline[tid] = (format[strlen(format) - 1] == '\n');
  return fprintf(output_file, "%s", format);
}

static inline int get_bflutindex(const int tempindex, const int element, const int ion, const int level,
                                 const int phixstargetindex) {
  const int contindex = -1 - globals::elements[element].ions[ion].levels[level].cont_index + phixstargetindex;

  const int bflutindex = tempindex * globals::nbfcontinua + contindex;
  assert_testmodeonly(bflutindex <= TABLESIZE * globals::nbfcontinua);
  return bflutindex;
}

#ifdef _OPENMP
#define safeadd(var, val) _Pragma("omp atomic update") var += val
#else
#define safeadd(var, val) var += val
#endif

#define safeincrement(var) safeadd(var, 1)

// #define DO_TITER

static inline void gsl_error_handler_printout(const char *reason, const char *file, int line, int gsl_errno) {
  if (gsl_errno != 18)  // roundoff error
  {
    printout("WARNING: gsl (%s:%d): %s (Error code %d)\n", file, line, reason, gsl_errno);
    // abort();
  }
}

static FILE *fopen_required(const char *filename, const char *mode) {
  assert_always(filename != nullptr);
  std::string datafolderfilename("data/");

  std::filesystem::path pid_file_path(datafolderfilename + filename);
  if (mode[0] == 'r' && std::filesystem::exists(datafolderfilename)) {
    return fopen_required(datafolderfilename.c_str(), mode);
  }
  FILE *file = std::fopen(filename, mode);
  if (file == nullptr) {
    printout("ERROR: Could not open file '%s' for mode '%s'.\n", filename, mode);
    abort();
  }

  return file;
}

static int get_timestep(const double time) {
  assert_always(time >= globals::tmin);
  assert_always(time < globals::tmax);
  for (int nts = 0; nts < globals::ntstep; nts++) {
    const double tsend = (nts < (globals::ntstep - 1)) ? globals::time_step[nts + 1].start : globals::tmax;
    if (time >= globals::time_step[nts].start && time < tsend) {
      return nts;
    }
  }
  assert_always(false);  // could not find matching timestep

  return -1;
}

inline int get_max_threads(void) {
#if defined _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

inline int get_num_threads(void) {
#if defined _OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}

inline int get_thread_num(void) {
#if defined _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

inline float rng_uniform(void) {
  float zrand;
  do {
    zrand = std::generate_canonical<float, std::numeric_limits<float>::digits>(stdrng);
  } while (zrand == 1.);
  return zrand;
}

inline float rng_uniform_pos(void) {
  float zrand = 0.;
  do {
    zrand = rng_uniform();
  } while (zrand <= 0.);
  return zrand;
}

inline void rng_init(const uint_fast64_t zseed) {
  printout("rng is a std::mt19937 generator\n");
  stdrng.seed(zseed);
}

inline bool is_pid_running(pid_t pid) {
  while (waitpid(-1, 0, WNOHANG) > 0) {
    // Wait for defunct....
  }

  if (0 == kill(pid, 0)) return true;  // Process exists

  return false;
}

inline void check_already_running(void) {
  pid_t artispid = getpid();

  if (std::filesystem::exists("artis.pid")) {
    std::ifstream pidfile("artis.pid", std::ifstream::in);
    pid_t artispid_in;
    pidfile >> artispid_in;
    pidfile.close();
    if (is_pid_running(artispid_in)) {
      fprintf(
          stderr,
          "\nERROR: artis or exspec is already running in this folder with existing pid %d. Refusing to start. (delete "
          "artis.pid if you are sure this is incorrect)\n",
          artispid_in);
      abort();
    }
  }

  std::ofstream pidfile("artis.pid", std::ofstream::out | std::ofstream::trunc);
  pidfile << artispid;
  pidfile.close();
}

#endif  // SN3D_H
