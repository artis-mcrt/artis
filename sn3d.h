#ifndef SN3D_H
#define SN3D_H

#include <gsl/gsl_integration.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <ranges>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#ifdef STACKTRACE_ON
#include <stacktrace>
#define STACKTRACEIFSUPPORTED << std::stacktrace::current()
#else
#define STACKTRACEIFSUPPORTED
#endif

#ifdef STDPAR_ON
#include <execution>
#include <thread>

#define EXEC_PAR_UNSEQ std::execution::par_unseq,
#define EXEC_PAR std::execution::par,
#else
#define EXEC_PAR_UNSEQ
#define EXEC_PAR
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#ifdef __NVCOMPILER_CUDA_ARCH__
#define THREADLOCALONHOST
#else
#define THREADLOCALONHOST thread_local static
#endif

#include "constants.h"

constexpr int cellcacheslotid = 0;
inline bool use_cellcache = false;

extern std::fstream output_file;

inline std::array<char, 1024> outputlinebuf = {};
inline std::string outputlinestr = {};
inline bool outputstartofline = true;
inline tm timebuf{};

#ifndef USE_SIMPSON_INTEGRATOR
#define USE_SIMPSON_INTEGRATOR false
#endif

inline thread_local auto gslworkspace =
    std::unique_ptr<gsl_integration_workspace, void (*)(gsl_integration_workspace*)>{
        USE_SIMPSON_INTEGRATOR ? nullptr : gsl_integration_workspace_alloc(GSLWSIZE),
        [](gsl_integration_workspace* const w) -> void {
          if (!USE_SIMPSON_INTEGRATOR) {
            gsl_integration_workspace_free(w);
          }
        }};

#ifdef _OPENMP

#ifndef GPU_ON
#pragma omp threadprivate(output_file, outputlinebuf, outputstartofline, timebuf)
#endif

#endif

#ifdef __NVCOMPILER_CUDA_ARCH__
#define printout(...) printf(__VA_ARGS__)

template <class... Args>
inline auto logprintfmt(const std::format_string<Args...> fmt, Args&&... args) -> void {
  printf("%s", std::format(fmt, std::forward<Args>(args)...).c_str());
}

template <class... Args>
inline auto logprintlnfmt(const std::format_string<Args...> fmt, Args&&... args) -> void {
  printf("%s\n", std::format(fmt, std::forward<Args>(args)...).c_str());
}

#define __artis_assert(e)                         \
  {                                               \
    const bool assertpass = static_cast<bool>(e); \
    assert(assertpass);                           \
  }

#else
inline void print_line_start() {
  if (outputstartofline) {
    const time_t now_time = time(nullptr);
    strftime(outputlinebuf.data(), 32, "%FT%TZ", gmtime_r(&now_time, &timebuf));
    output_file << outputlinebuf.data() << ' ';
  }
}

__attribute__((__format__(__printf__, 1, 2))) inline auto printout(const char* format, ...) -> void {
  print_line_start();
  va_list args{};
  va_start(args, format);
  vsnprintf(outputlinebuf.data(), outputlinebuf.size(), format, args);
  va_end(args);

  const auto linebuflen = strlen(outputlinebuf.data());
  outputstartofline = (linebuflen == 0 || (outputlinebuf[linebuflen - 1] == '\n'));
  output_file << outputlinebuf.data();
  output_file.flush();
}

template <class... Args>
inline auto logprintfmt(const std::format_string<Args...> fmt, Args&&... args) -> void {
  print_line_start();
  outputlinestr = std::format(fmt, std::forward<Args>(args)...);
  outputstartofline = (outputlinestr.back() == '\n');
  output_file << outputlinestr;
  output_file.flush();
}

template <class... Args>
inline auto logprintlnfmt(const std::format_string<Args...> fmt, Args&&... args) -> void {
  print_line_start();
  outputlinestr = std::format(fmt, std::forward<Args>(args)...);
  outputstartofline = true;
  output_file << outputlinestr << '\n';
  output_file.flush();
}

#define __artis_assert(e)                                                                                              \
  {                                                                                                                    \
    const bool assertpass = static_cast<bool>(e);                                                                      \
    if (!assertpass) [[unlikely]] {                                                                                    \
      if (output_file) {                                                                                               \
        output_file << "\n[rank " << globals::my_rank << "] " << __FILE__ << ":" << __LINE__ << ": failed assertion `" \
                    << #e << "` in function " << __PRETTY_FUNCTION__ << '\n';                                          \
        output_file.flush();                                                                                           \
      }                                                                                                                \
      std::cerr << "\n[rank " << globals::my_rank << "] " << __FILE__ << ":" << __LINE__ << ": failed assertion `"     \
                << #e << "` in function " << __PRETTY_FUNCTION__ << '\n' STACKTRACEIFSUPPORTED;                        \
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

#ifdef _OPENMP
#define atomicadd(var, val)                  \
  {                                          \
    _Pragma("omp atomic update") var += val; \
  }

#else
#ifdef STDPAR_ON

#include <atomic>
template <typename T, typename U>
constexpr void atomicadd(T& var, U&& val) {
  std::atomic_ref<T>(var).fetch_add(std::forward<U>(val), std::memory_order_relaxed);
}
#else
#define atomicadd(var, val) var += (val);
#endif
#endif

inline void gsl_error_handler_printout(const char* reason, const char* file, int line, int gsl_errno) {
  if (gsl_errno != 18)  // roundoff error
  {
    printout("WARNING: gsl (%s:%d): %s (Error code %d)\n", file, line, reason, gsl_errno);
  }
}

[[nodiscard]] inline auto fopen_required(const std::string& filename, std::span<const char> mode) -> FILE* {
  // look in the data folder first
  const std::string datafolderfilename = "data/" + filename;
  if (mode[0] == 'r' && std::filesystem::exists(datafolderfilename)) {
    return fopen_required(datafolderfilename, mode);
  }

  auto* file = std::fopen(filename.c_str(), mode.data());
  if (file == nullptr) {
    logprintlnfmt("ERROR: Could not open file '{}' for mode '{}'.", filename, mode.data());
    std::abort();
  }

  return file;
}

[[nodiscard]] inline auto fopen_required_uniqueptr(const std::string& filename, std::span<const char> mode) {
  return std::unique_ptr<FILE, int (*)(FILE*)>(fopen_required(filename, mode),
                                               [](FILE* fp) -> int { return std::fclose(fp); });
}

[[nodiscard]] inline auto fstream_required(const std::string& filename, std::ios_base::openmode mode) -> std::fstream {
  if (filename.empty()) {
    logprintlnfmt("ERROR: Cannot open file with empty filename.");
    std::abort();
  }
  const std::string datafolderfilename = "data/" + filename;
  if (mode == std::ios::in && std::filesystem::exists(datafolderfilename)) {
    return fstream_required(datafolderfilename, mode);
  }
  auto file = std::fstream(filename, mode);
  if (!file.is_open()) {
    logprintlnfmt("ERROR: Could not open file '{}'", filename);
    std::abort();
  }
  return file;
}

#include "globals.h"

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
#ifdef STDPAR_ON
  return static_cast<int>(std::thread::hardware_concurrency());
#else
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
#endif
}

[[nodiscard]] inline auto get_thread_num() -> int {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
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
      std::string line;
      std::getline(pidfile, line);
      std::istringstream(line) >> artispid_in;
      std::getline(pidfile, line);
      std::string working_directory;
      std::istringstream(line) >> working_directory;
      pidfile.close();
      if (is_pid_running(artispid_in) && std::filesystem::current_path().generic_string() == working_directory) {
        fprintf(stderr,
                "\nERROR: artis or exspec is already running in this folder with existing pid %d. Refusing to start. "
                "(delete artis.pid if you are sure this is incorrect)\n",
                artispid_in);
        std::abort();
      }
    }

    auto pidfile = std::fstream("artis.pid", std::ofstream::out | std::ofstream::trunc);
    pidfile << artispid << '\n';
    pidfile << std::filesystem::current_path().generic_string() << '\n';
  }

  // make sure rank 0 checked for a pid file before we proceed
  MPI_Barrier(MPI_COMM_WORLD);
}

constexpr auto get_range_chunk(const ptrdiff_t size, const ptrdiff_t nchunks, const ptrdiff_t nchunk)
    -> std::tuple<ptrdiff_t, ptrdiff_t> {
  assert_always(size >= 0);
  assert_always(nchunks >= 0);
  assert_always(nchunk >= 0);
  const auto minchunksize = size / nchunks;  // integer division, minimum non-empty cells per process
  const auto n_remainder = size % nchunks;
  const auto nstart =
      ((minchunksize + 1) * std::min(n_remainder, nchunk)) + (minchunksize * std::max(0Z, nchunk - n_remainder));
  const auto nsize = (nchunk < n_remainder) ? minchunksize + 1 : minchunksize;
  assert_testmodeonly(nstart >= 0);
  assert_testmodeonly(nsize >= 0);
  assert_testmodeonly((nstart + nsize) <= size);
  return std::tuple{nstart, nsize};
}

static_assert(get_range_chunk(10, 3, 0) == std::tuple{0, 4});
static_assert(get_range_chunk(10, 3, 1) == std::tuple{4, 3});
static_assert(get_range_chunk(10, 3, 2) == std::tuple{7, 3});

template <typename T>
[[nodiscard]] auto MPI_shared_malloc_span_keepwin(const ptrdiff_t num_allranks) -> std::tuple<std::span<T>, MPI_Win> {
  if (num_allranks == 0) {
    return {std::span<T>{}, MPI_WIN_NULL};
  }
  assert_always(num_allranks >= 0);

  const auto [_, num_thisnoderank] = get_range_chunk(num_allranks, globals::node_nprocs, globals::rank_in_node);
  assert_always(num_thisnoderank >= 0);

  auto size = static_cast<MPI_Aint>(num_thisnoderank * sizeof(T));
  int disp_unit = sizeof(T);
  MPI_Win mpiwin{MPI_WIN_NULL};
  T* ptr{};

  assert_always(MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &ptr, &mpiwin) ==
                MPI_SUCCESS);
  assert_always(MPI_Win_shared_query(mpiwin, 0, &size, &disp_unit, &ptr) == MPI_SUCCESS);
  MPI_Barrier(globals::mpi_comm_node);
  assert_always(ptr != nullptr);
#pragma clang unsafe_buffer_usage begin
  return {std::span<T>(ptr, num_allranks), mpiwin};
#pragma clang unsafe_buffer_usage end
}

template <typename T>
[[nodiscard]] auto MPI_shared_malloc_span(const ptrdiff_t num_allranks) -> std::span<T> {
  return std::get<0>(MPI_shared_malloc_span_keepwin<T>(num_allranks));
}

template <typename T>
inline auto GET_MPI_TYPE() -> MPI_Datatype {
  if constexpr (std::is_same_v<T, float>) {
    return MPI_FLOAT;
  } else if constexpr (std::is_same_v<T, double>) {
    return MPI_DOUBLE;
  } else if constexpr (std::is_same_v<T, int>) {
    return MPI_INT;
  } else if constexpr (std::is_same_v<T, std::int64_t>) {
    return MPI_INT64_T;
  } else if constexpr (std::is_same_v<T, bool>) {
    return MPI_C_BOOL;
  } else {
    return MPI_BYTE;  // fallback to byte type for unsupported types
  }
}

// these wrappers add type, bounds, and overflow safety to the MPI calls
template <typename R, typename Op, typename Comm>
  requires std::ranges::random_access_range<R>
inline void MPI_Allreduce_safe(R&& data, Op&& op, Comm&& comm) {
  if (std::forward<R>(data).empty()) {
    return;
  }
  assert_always(!std::forward<R>(data).empty());
  assert_always(std::forward<R>(data).data() != nullptr);
  assert_always(op != MPI_OP_NULL);
  assert_always(comm != MPI_COMM_NULL);

  const auto true_size = std::ssize(std::forward<R>(data));
  const auto int_data_size = static_cast<int>(true_size);
  assert_always(std::cmp_equal(int_data_size, true_size));

  const auto mpi_datatype = GET_MPI_TYPE<std::ranges::range_value_t<R>>();
  assert_always(mpi_datatype != MPI_BYTE);  // we can't reduce MPI_BYTE types

  auto ret = MPI_Allreduce(MPI_IN_PLACE, std::forward<R>(data).data(), int_data_size, mpi_datatype,
                           std::forward<Op>(op), std::forward<Comm>(comm));
  assert_always(ret == MPI_SUCCESS);
}

template <typename T, typename Op, typename Comm>
  requires(!std::ranges::random_access_range<T>)
inline void MPI_Allreduce_safe(T& data, Op&& op, Comm&& comm) {
  MPI_Allreduce_safe(std::span(&data, 1), std::forward<Op>(op), std::forward<Comm>(comm));
}

template <typename R, typename Comm>
  requires std::ranges::random_access_range<R>
inline void MPI_Bcast_safe(R&& data, const int root, Comm&& comm) {
  if (std::forward<R>(data).empty()) {
    return;
  }
  assert_always(std::forward<R>(data).data() != nullptr);
  assert_always(comm != MPI_COMM_NULL);
  using value_t = typename std::ranges::range_value_t<R>;

  const auto mpi_datatype = GET_MPI_TYPE<value_t>();
  // if we're transferring bytes, then we need multiply the array count by the byte size of the type
  const ptrdiff_t sizefactor = mpi_datatype == MPI_BYTE ? sizeof(value_t) : 1;

  const auto true_size = std::ssize(std::forward<R>(data)) * sizefactor;
  const auto int_data_size = static_cast<int>(true_size);
  assert_always(std::cmp_equal(int_data_size, true_size));

  auto ret = MPI_Bcast(std::forward<R>(data).data(), int_data_size, mpi_datatype, root, std::forward<Comm>(comm));
  assert_always(ret == MPI_SUCCESS);
}

template <typename T, typename Comm>
  requires(!std::ranges::random_access_range<T>)
inline void MPI_Bcast_safe(T& data, const int root, Comm&& comm) {
  MPI_Bcast_safe(std::span(&data, 1), root, std::forward<Comm>(comm));
}

template <typename R, typename Op, typename Comm>
  requires std::ranges::random_access_range<R>
inline void MPI_Reduce_safe(R&& data, Op&& op, const int root, Comm&& comm) {
  if (std::forward<R>(data).empty()) {
    return;
  }
  assert_always(std::forward<R>(data).data() != nullptr);
  assert_always(comm != MPI_COMM_NULL);

  const auto true_size = std::ssize(std::forward<R>(data));
  const auto int_data_size = static_cast<int>(true_size);
  assert_always(std::cmp_equal(int_data_size, true_size));

  int my_rank{-1};
  assert_always(MPI_Comm_rank(std::forward<Comm>(comm), &my_rank) == MPI_SUCCESS);
  assert_always(my_rank >= 0);

  const auto mpi_datatype = GET_MPI_TYPE<std::ranges::range_value_t<R>>();
  assert_always(mpi_datatype != MPI_BYTE);  // we can't reduce MPI_BYTE types

  const auto ret = MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : std::forward<R>(data).data(), std::forward<R>(data).data(),
                              int_data_size, mpi_datatype, std::forward<Op>(op), root, std::forward<Comm>(comm));
  assert_always(ret == MPI_SUCCESS);
}

template <typename T>
constexpr void resize_exactly(std::vector<T>& vec, const size_t size) {
  // just resizing can (only with libstdc++?) allocate a larger capacity than needed
  vec.reserve(size);
  vec.resize(size);
}

#endif  // SN3D_H
