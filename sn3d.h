#pragma once
#include <string_view>
#ifndef SN3D_H
#define SN3D_H

#include <getopt.h>
#include <gsl/gsl_integration.h>
#include <sys/wait.h>
#include <unistd.h>

#ifdef STDPAR_ON
#include <execution>

#ifndef __cpp_lib_execution
// homebrew llvm doesn't support execution policy yet, so brew install onedpl tbb
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#endif

#define EXEC_PAR_UNSEQ std::execution::par_unseq,
#define EXEC_PAR std::execution::par,
#else
#define EXEC_PAR_UNSEQ
#define EXEC_PAR
#endif

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

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef MPI_ON
#include <mpi.h>
#endif

extern std::ofstream output_file;
inline constexpr int cellcacheslotid = 0;
extern bool use_cellcache;

extern long long int rngseed;

extern gsl_integration_workspace *gslworkspace;

#define __artis_assert(e)                                                                                              \
  {                                                                                                                    \
    const bool pass = static_cast<bool>(e);                                                                            \
    if (!pass) {                                                                                                       \
      if (output_file) {                                                                                               \
        output_file << "\n[rank " << globals::rank_global << "] " << __FILE__ << ":" << __LINE__                       \
                    << ": failed assertion `" << #e << "` in function " << __PRETTY_FUNCTION__ << "\n";                \
        output_file.flush();                                                                                           \
      }                                                                                                                \
      std::cerr << "\n[rank " << globals::rank_global << "] " << __FILE__ << ":" << __LINE__ << ": failed assertion `" \
                << #e << "` in function " << __PRETTY_FUNCTION__ << "\n";                                              \
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

// make these thread_local if we want separate log files for STDPAR threads
inline char outputlinebuf[1024] = "";
inline bool outputstartofline = true;
inline struct tm timebuf {};

#ifdef _OPENMP
#pragma omp threadprivate(cellcacheslotid, rngseed, gslworkspace, output_file, outputlinebuf, outputstartofline, \
                              timebuf)
#endif

static inline void print_line_start() {
  if (outputstartofline) {
    const time_t now_time = time(nullptr);
    strftime(outputlinebuf, 32, "%FT%TZ", gmtime_r(&now_time, &timebuf));
    output_file << outputlinebuf << " ";
  }
}

#define printout(...)                                                       \
  {                                                                         \
    print_line_start();                                                     \
    snprintf(outputlinebuf, 1024, __VA_ARGS__);                             \
    outputstartofline = (outputlinebuf[strlen(outputlinebuf) - 1] == '\n'); \
    output_file << outputlinebuf;                                           \
    output_file.flush();                                                    \
  }

#include "globals.h"

[[nodiscard]] static inline auto get_bflutindex(const int tempindex, const int element, const int ion, const int level,
                                                const int phixstargetindex) -> int {
  const int contindex = -1 - globals::elements[element].ions[ion].levels[level].cont_index + phixstargetindex;

  const int bflutindex = tempindex * globals::nbfcontinua + contindex;
  assert_testmodeonly(bflutindex >= 0);
  assert_testmodeonly(bflutindex <= TABLESIZE * globals::nbfcontinua);
  return bflutindex;
}

template <typename T>
inline void safeadd(T &var, T val) {
#ifdef _OPENMP
#pragma omp atomic update
  var += val;
#else
#ifdef STDPAR_ON
#ifdef __cpp_lib_atomic_ref
  static_assert(std::atomic<T>::is_always_lock_free);
  std::atomic_ref<T>(var).fetch_add(val, std::memory_order_relaxed);
#else
  // this works on clang but not gcc for doubles.
  __atomic_fetch_add(&var, val, __ATOMIC_RELAXED);
#endif
#else
  var += val;
#endif
#endif
}

#define safeincrement(var) safeadd((var), 1)

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

[[nodiscard]] inline auto get_thread_num() -> int {
#if defined _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

inline auto rng_uniform() -> float {
  thread_local std::mt19937 stdrng(rngseed);
  while (true) {
    const auto zrand = std::generate_canonical<float, std::numeric_limits<float>::digits>(stdrng);
    if (zrand != 1.) {
      return zrand;
    }
  }
}

inline auto rng_uniform_pos() -> float {
  while (true) {
    const auto zrand = rng_uniform();
    if (zrand > 0) {
      return zrand;
    }
  }
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

constexpr auto get_range_chunk(int size, int nchunks, int nchunk) -> std::tuple<int, int> {
  const int minchunksize = size / nchunks;  // integer division, minimum non-empty cells per process
  const int n_remainder = size % nchunks;
  const auto nstart =
      (minchunksize + 1) * std::min(n_remainder, nchunk) + minchunksize * std::max(0, nchunk - n_remainder);
  const auto nsize = (nchunk < n_remainder) ? minchunksize + 1 : minchunksize;
  return {nstart, nsize};
}

#endif  // SN3D_H
