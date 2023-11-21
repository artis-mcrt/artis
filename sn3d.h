#pragma once
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
#include "atomic.h"
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
extern bool use_cellcache;

extern std::mt19937 stdrng;

extern gsl_integration_workspace *gslworkspace;

#ifdef _OPENMP
#pragma omp threadprivate(tid, use_cellcache, stdrng, gslworkspace, output_file)
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
      std::abort();                                                                                                    \
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

template <typename... Args>
static auto printout(const char *format, Args... args) -> int {
  if (globals::startofline[tid]) {
    const time_t now_time = time(nullptr);
    char s[32] = "";
    strftime(s, 32, "%FT%TZ", gmtime(&now_time));  // NOLINT[concurrency-mt-unsafe]
    fprintf(output_file, "%s ", s);
  }
  globals::startofline[tid] = (format[strlen(format) - 1] == '\n');
  return fprintf(output_file, format, args...);
}

static auto printout(const char *format) -> int {
  if (globals::startofline[tid]) {
    const time_t now_time = time(nullptr);
    char s[32] = "";
    strftime(s, 32, "%FT%TZ", gmtime(&now_time));  // NOLINT[concurrency-mt-unsafe]
    fprintf(output_file, "%s ", s);
  }
  globals::startofline[tid] = (format[strlen(format) - 1] == '\n');
  return fprintf(output_file, "%s", format);
}

[[nodiscard]] static inline auto get_bflutindex(const int tempindex, const int element, const int ion, const int level,
                                                const int phixstargetindex) -> int {
  const int contindex = -1 - globals::elements[element].ions[ion].levels[level].cont_index + phixstargetindex;

  const int bflutindex = tempindex * globals::nbfcontinua + contindex;
  assert_testmodeonly(bflutindex <= TABLESIZE * globals::nbfcontinua);
  return bflutindex;
}

#ifdef _OPENMP
#define safeadd(var, val) _Pragma("omp atomic update") var += val
#else
#define safeadd(var, val) var = (var) + val
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

static auto fopen_required(const std::string &filename, const char *mode) -> FILE * {
  // look in the data folder first
  const std::string datafolderfilename = "data/" + filename;
  if (mode[0] == 'r' && std::filesystem::exists(datafolderfilename)) {
    return fopen_required(datafolderfilename, mode);
  }

  FILE *file = std::fopen(filename.c_str(), mode);
  if (file == nullptr) {
    printout("ERROR: Could not open file '%s' for mode '%s'.\n", filename.c_str(), mode);
    std::abort();
  }

  return file;
}

static auto fstream_required(const std::string &filename, std::ios_base::openmode mode) -> std::fstream {
  const std::string datafolderfilename = "data/" + filename;
  if (mode == std::ios::in && std::filesystem::exists(datafolderfilename)) {
    return fstream_required(datafolderfilename, mode);
  }
  auto file = std::fstream(filename, mode);
  if (!file.is_open()) {
    printout("ERROR: Could not open file '%s'\n", filename.c_str());
    std::abort();
  }
  return file;
}

[[nodiscard]] static auto get_timestep(const double time) -> int {
  assert_always(time >= globals::tmin);
  assert_always(time < globals::tmax);
  for (int nts = 0; nts < globals::ntimesteps; nts++) {
    const double tsend = (nts < (globals::ntimesteps - 1)) ? globals::timesteps[nts + 1].start : globals::tmax;
    if (time >= globals::timesteps[nts].start && time < tsend) {
      return nts;
    }
  }
  assert_always(false);  // could not find matching timestep

  return -1;
}

[[nodiscard]] inline auto get_max_threads() -> int {
#if defined _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

[[nodiscard]] inline auto get_num_threads() -> int {
#if defined _OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}

[[nodiscard]] inline auto get_thread_num() -> int {
#if defined _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

inline auto rng_uniform() -> float {
  const auto zrand = std::generate_canonical<float, std::numeric_limits<float>::digits>(stdrng);
  if (zrand == 1.) {
    return rng_uniform();
  }
  return zrand;
}

inline auto rng_uniform_pos() -> float {
  const auto zrand = std::generate_canonical<float, std::numeric_limits<float>::digits>(stdrng);
  if (zrand <= 0.) {
    return rng_uniform_pos();
  }
  return zrand;
}

inline void rng_init(const uint_fast64_t zseed) {
  printout("rng is a std::mt19937 generator\n");
  stdrng.seed(zseed);
}

[[nodiscard]] inline auto is_pid_running(pid_t pid) -> bool {
  while (waitpid(-1, nullptr, WNOHANG) > 0) {
    // Wait for defunct....
  }

  return (kill(pid, 0) == 0);
}

inline void check_already_running() {
  const pid_t artispid = getpid();

  if (std::filesystem::exists("artis.pid")) {
    auto pidfile = std::fstream("artis.pid", std::ios::in);
    pid_t artispid_in = 0;
    pidfile >> artispid_in;
    pidfile.close();
    if (is_pid_running(artispid_in)) {
      fprintf(stderr,
              "\nERROR: artis or exspec is already running in this folder with existing pid %d. Refusing to start. "
              "(delete "
              "artis.pid if you are sure this is incorrect)\n",
              artispid_in);
      std::abort();
    }
  }

  auto pidfile = std::fstream("artis.pid", std::ofstream::out | std::ofstream::trunc);
  pidfile << artispid;
  pidfile.close();
}

[[nodiscard]] inline auto get_ionestimindex(const int mgi, const int element, const int ion) -> int {
  return mgi * get_includedions() + get_uniqueionindex(element, ion);
}

#endif  // SN3D_H
