#ifndef SN3D_H
#define SN3D_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <getopt.h>
#include <gsl/gsl_integration.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#ifndef GPU_ON
#include <random>
#endif
#include <string>
#include <tuple>

#ifdef STDPAR_ON
#include <execution>

// #ifndef __cpp_lib_execution
// // homebrew llvm doesn't support execution policy yet, so brew install onedpl tbb
// #include <oneapi/dpl/algorithm>
// #include <oneapi/dpl/execution>
// #endif

#define EXEC_PAR_UNSEQ std::execution::par_unseq,
#define EXEC_PAR std::execution::par,
#else
#define EXEC_PAR_UNSEQ
#define EXEC_PAR
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef MPI_ON
#include <mpi.h>
#endif

#ifdef __NVCOMPILER_CUDA_ARCH__
#define THREADLOCALONHOST
#else
#define THREADLOCALONHOST thread_local static
#endif

#include "constants.h"

constexpr int cellcacheslotid = 0;
inline bool use_cellcache = false;

#ifndef GPU_ON
extern std::mt19937 stdrng;
#endif

extern std::ofstream output_file;

inline char outputlinebuf[1024] = "";
inline bool outputstartofline = true;
inline tm timebuf{};

// if not set, force Simpson integrator on GPU mode (since gsl doesn't work there!)
#ifndef USE_SIMPSON_INTEGRATOR
#define USE_SIMPSON_INTEGRATOR false
#endif

inline thread_local auto gslworkspace =
    std::unique_ptr<gsl_integration_workspace, void (*)(gsl_integration_workspace *)>{
        USE_SIMPSON_INTEGRATOR ? nullptr : gsl_integration_workspace_alloc(GSLWSIZE),
        USE_SIMPSON_INTEGRATOR ? [](gsl_integration_workspace *const w) {} : gsl_integration_workspace_free};

#ifdef _OPENMP

#ifndef GPU_ON
#pragma omp threadprivate(stdrng, output_file, outputlinebuf, outputstartofline, timebuf)
#endif

#endif

#ifdef __NVCOMPILER_CUDA_ARCH__
#define printout(...) printf(__VA_ARGS__)

#define __artis_assert(e)                         \
  {                                               \
    const bool assertpass = static_cast<bool>(e); \
    assert(assertpass);                           \
  }

#else
inline void print_line_start() {
  if (outputstartofline) {
    const time_t now_time = time(nullptr);
    strftime(outputlinebuf, 32, "%FT%TZ", gmtime_r(&now_time, &timebuf));
    output_file << outputlinebuf << " ";
  }
}

__attribute__((__format__(__printf__, 1, 2))) inline auto printout(const char *format, ...) -> void {
  print_line_start();
  va_list args{};
  va_start(args, format);
  vsnprintf(outputlinebuf, sizeof(outputlinebuf), format, args);
  va_end(args);

  outputstartofline = (outputlinebuf[strlen(outputlinebuf) - 1] == '\n');
  output_file << outputlinebuf;
  output_file.flush();
}

#define __artis_assert(e)                                                                                              \
  {                                                                                                                    \
    const bool assertpass = static_cast<bool>(e);                                                                      \
    if (!assertpass) {                                                                                                 \
      if (output_file) {                                                                                               \
        output_file << "\n[rank " << globals::my_rank << "] " << __FILE__ << ":" << __LINE__ << ": failed assertion `" \
                    << #e << "` in function " << __PRETTY_FUNCTION__ << "\n";                                          \
        output_file.flush();                                                                                           \
      }                                                                                                                \
      std::cerr << "\n[rank " << globals::my_rank << "] " << __FILE__ << ":" << __LINE__ << ": failed assertion `"     \
                << #e << "` in function " << __PRETTY_FUNCTION__ << "\n";                                              \
      std::abort();                                                                                                    \
    }                                                                                                                  \
    assert(assertpass);                                                                                                \
  }
#endif

#define assert_always(e) __artis_assert(e)

#ifndef TESTMODE
#define TESTMODE false
#endif

#if defined TESTMODE && TESTMODE
#define assert_testmodeonly(e) __artis_assert(e)
#else
#define assert_testmodeonly(e) (void)0
#endif

#if defined REPRODUCIBLE && REPRODUCIBLE
#define SORT_OR_STABLE_SORT stable_sort
#else
#define SORT_OR_STABLE_SORT sort
#endif

template <typename T>
inline void atomicadd(T &var, const T &val) {
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

inline void gsl_error_handler_printout(const char *reason, const char *file, int line, int gsl_errno) {
  if (gsl_errno != 18)  // roundoff error
  {
    printout("WARNING: gsl (%s:%d): %s (Error code %d)\n", file, line, reason, gsl_errno);
    // abort();
  }
}

[[nodiscard]] inline auto fopen_required(const std::string &filename, const char *mode) -> FILE * {
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

[[nodiscard]] inline auto fstream_required(const std::string &filename, std::ios_base::openmode mode) -> std::fstream {
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
#include "globals.h"

[[nodiscard]] inline auto get_bflutindex(const int tempindex, const int element, const int ion, const int level,
                                         const int phixstargetindex) -> int {
  const int contindex = globals::elements[element].ions[ion].levels[level].cont_index + phixstargetindex;

  const int bflutindex = (tempindex * globals::nbfcontinua) + contindex;
  assert_testmodeonly(bflutindex >= 0);
  assert_testmodeonly(bflutindex <= TABLESIZE * globals::nbfcontinua);
  return bflutindex;
}

[[nodiscard]] inline auto get_timestep(const double time) -> int {
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
  while (true) {
#ifndef GPU_ON
    const auto zrand = std::generate_canonical<float, std::numeric_limits<float>::digits>(stdrng);
#else
    const auto zrand = 0.5;
#endif
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
  if (globals::my_rank == 0) {
    const pid_t artispid = getpid();

    if (std::filesystem::exists("artis.pid")) {
      auto pidfile = std::fstream("artis.pid", std::ios::in);
      pid_t artispid_in = 0;
      pidfile >> artispid_in;
      pidfile.close();
      if (is_pid_running(artispid_in)) {
        fprintf(stderr,
                "\nERROR: artis or exspec is already running in this folder with existing pid %d. Refusing to start. "
                "(delete artis.pid if you are sure this is incorrect)\n",
                artispid_in);
        std::abort();
      }
    }

    auto pidfile = std::fstream("artis.pid", std::ofstream::out | std::ofstream::trunc);
    pidfile << artispid;
    pidfile.close();
  }

// make sure rank 0 checked for a pid file before we proceed
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

constexpr auto get_range_chunk(const ptrdiff_t size, const ptrdiff_t nchunks, const ptrdiff_t nchunk)
    -> std::tuple<ptrdiff_t, ptrdiff_t> {
  assert_always(size >= 0);
  assert_always(nchunks >= 0);
  assert_always(nchunk >= 0);
  const auto minchunksize = size / nchunks;  // integer division, minimum non-empty cells per process
  const auto n_remainder = size % nchunks;
  const auto nstart = ((minchunksize + 1) * std::min(n_remainder, nchunk)) +
                      (minchunksize * std::max(static_cast<ptrdiff_t>(0), nchunk - n_remainder));
  const auto nsize = (nchunk < n_remainder) ? minchunksize + 1 : minchunksize;
  assert_testmodeonly(nstart >= 0);
  assert_testmodeonly(nsize >= 0);
  assert_testmodeonly((nstart + nsize) <= size);
  return std::tuple{nstart, nsize};
}

#endif  // SN3D_H
