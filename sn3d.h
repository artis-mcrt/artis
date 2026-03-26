#ifndef SN3D_H
#define SN3D_H

#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <atomic>
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
#include <limits>
#include <memory>
#include <ranges>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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

#include "constants.h"

constexpr int cellcacheslotid = 0;

extern std::fstream output_file;

inline std::array<char, 1024> outputlinebuf{};
inline bool outputstartofline = true;

#ifdef _OPENMP
#ifndef GPU_ON
#pragma omp threadprivate(output_file, outputlinebuf, outputstartofline)
#endif
#endif

#ifdef __NVCOMPILER_CUDA_ARCH__
#include <string_view>

#define printout(...) printf(__VA_ARGS__)

template <typename... Args>
inline auto printlog(std::string_view fmt, Args&&... args) -> void {
  const auto str = std::vformat(fmt, std::make_format_args(args...));
  printf("%s", str.c_str());
}

template <typename... Args>
inline auto printlnlog(std::string_view fmt, Args&&... args) -> void {
  const auto str = std::vformat(fmt, std::make_format_args(args...));
  printf("%s\n", str.c_str());
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
    THREADLOCALONHOST tm timebuf{};
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

template <typename... Args>
inline auto printlog(const std::format_string<Args...> fmt, Args&&... args) -> void {
  print_line_start();
  THREADLOCALONHOST std::string outputlinestr;
  outputlinestr = std::format(fmt, std::forward<Args>(args)...);
  outputstartofline = (outputlinestr.back() == '\n');
  output_file << outputlinestr;
  output_file.flush();
}

template <typename... Args>
inline auto printlnlog(const std::format_string<Args...> fmt, Args&&... args) -> void {
  print_line_start();
  outputstartofline = true;
  output_file << std::format(fmt, std::forward<Args>(args)...) << '\n';
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

#elifdef STDPAR_ON

template <typename T, typename U>
constexpr void atomicadd(T& var, U&& val) {
  std::atomic_ref<T>(var).fetch_add(std::forward<U>(val), std::memory_order_relaxed);
}

#else

#define atomicadd(var, val) var += (val);

#endif

[[nodiscard]] inline auto fopen_required(const std::string& filename, std::span<const char> mode) -> FILE* {
  if (mode[0] == 'r') {
    // search data folders in order to find file to read
    for (const auto& datadir : datafolders) {
      const std::string datafolderfilename = std::string(datadir) + filename;
      auto* file = std::fopen(datafolderfilename.c_str(), mode.data());
      if (file != nullptr) {
        return file;
      }
    }
  } else {
    auto* file = std::fopen(filename.c_str(), mode.data());
    if (file != nullptr) {
      return file;
    }
  }

  printlnlog("ERROR: Could not open file '{}' for mode '{}'.", filename, mode.data());
  std::abort();
}

[[nodiscard]] inline auto fopen_required_uniqueptr(const std::string& filename, std::span<const char> mode) {
  return std::unique_ptr<FILE, int (*)(FILE*)>(fopen_required(filename, mode),
                                               [](FILE* fp) -> int { return std::fclose(fp); });
}

[[nodiscard]] inline auto fstream_required(const std::string& filename, std::ios_base::openmode mode) -> std::fstream {
  if (filename.empty()) {
    printlnlog("ERROR: Cannot open file with empty filename.");
    std::abort();
  }

  if (mode == std::ios::in) {
    // search data folders in order to find file to read
    for (const auto& datadir : datafolders) {
      auto datafolderfilename = std::string(datadir) + filename;
      auto file = std::fstream(datafolderfilename.c_str(), mode);
      if (file.is_open()) {
        return file;
      }
    }
  } else {
    // don't prepend data folders when writing
    auto file = std::fstream(filename, mode);
    if (file.is_open()) {
      return file;
    }
  }

  printlnlog("ERROR: Could not open file '{}'", filename);
  std::abort();
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
      std::istringstream{line} >> artispid_in;
      std::getline(pidfile, line);
      std::string working_directory;
      std::istringstream{line} >> working_directory;
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
[[nodiscard]] auto MPI_shared_malloc_span_keepwin(const ptrdiff_t num_allranks, const T& initval = {})
    -> std::tuple<std::span<T>, MPI_Win> {
  if (num_allranks == 0) {
    return {std::span<T>{}, MPI_WIN_NULL};
  }
  assert_always(num_allranks >= 0);

  // only rank_in_node 0 on each node allocates memory, but all ranks will get a pointer to it
  const auto num_thisnoderank = (globals::rank_in_node == 0) ? num_allranks : 0;

  auto size = static_cast<MPI_Aint>(num_thisnoderank * sizeof(T));
  int disp_unit = sizeof(T);
  MPI_Win mpiwin{MPI_WIN_NULL};
  T* ptr{};
  MPI_Info info{};
  MPI_Info_create(&info);
  // Request 64-byte alignment (e.g., for AVX-512)
  MPI_Info_set(info, "mpi_minimum_memory_alignment", "64");
  assert_always(MPI_Win_allocate_shared(size, disp_unit, info, globals::mpi_comm_node, static_cast<void*>(&ptr),
                                        &mpiwin) == MPI_SUCCESS);
  MPI_Info_free(&info);
  assert_always(MPI_Win_shared_query(mpiwin, 0, &size, &disp_unit, static_cast<void*>(&ptr)) == MPI_SUCCESS);
  assert_always(ptr != nullptr);
#ifdef __cpp_lib_is_sufficiently_aligned
  assert_always(std::is_sufficiently_aligned<64>(ptr));
#endif
#pragma clang unsafe_buffer_usage begin
  const auto newspan = std::span<T>(ptr, num_allranks);
#pragma clang unsafe_buffer_usage end
  // initialise the shared memory
  if (globals::rank_in_node == 0) {
    std::ranges::fill(newspan, initval);
  }
  MPI_Barrier(globals::mpi_comm_node);
  return {newspan, mpiwin};
}

template <typename T>
[[nodiscard]] auto MPI_shared_malloc_span(const ptrdiff_t num_allranks, const T& initval = {}) -> std::span<T> {
  return std::get<0>(MPI_shared_malloc_span_keepwin<T>(num_allranks, initval));
}

template <typename T>
inline auto GET_MPI_TYPE() -> MPI_Datatype {
  if constexpr (std::is_same_v<T, float>) {
    return MPI_FLOAT;
  } else if constexpr (std::is_same_v<T, double>) {
    return MPI_DOUBLE;
  } else if constexpr (std::is_same_v<T, std::int8_t>) {
    return MPI_INT8_T;
  } else if constexpr (std::is_same_v<T, std::int16_t>) {
    return MPI_INT16_T;
  } else if constexpr (std::is_same_v<T, std::int32_t>) {
    return MPI_INT32_T;
  } else if constexpr (std::is_same_v<T, std::int64_t>) {
    return MPI_INT64_T;
  } else if constexpr (std::is_same_v<T, bool>) {
    return MPI_C_BOOL;
  } else {
    return MPI_BYTE;  // fallback to byte type for unsupported types
  }
}

// MPI operations use a 32-bit int for the count, so we need to chunk large arrays
constexpr std::ptrdiff_t MPI_COUNT_MAX = std::numeric_limits<int>::max();

// these wrappers add type, bounds, and overflow safety to the MPI calls
template <typename R>
  requires std::ranges::random_access_range<R>
inline void MPI_Allreduce_safe(R&& data, MPI_Op op, MPI_Comm comm) {
  auto dataspan = std::span{std::forward<R>(data)};
  if (dataspan.empty()) {
    return;
  }
  assert_always(dataspan.data() != nullptr);
  assert_always(op != MPI_OP_NULL);
  assert_always(comm != MPI_COMM_NULL);

  const auto mpi_datatype = GET_MPI_TYPE<std::ranges::range_value_t<R>>();
  assert_always(mpi_datatype != MPI_BYTE);  // we can't reduce MPI_BYTE types

  const auto nchunks = (std::ssize(dataspan) / MPI_COUNT_MAX) + ((std::ssize(dataspan) % MPI_COUNT_MAX) == 0 ? 0 : 1);
  assert_always(nchunks >= 1);
  std::ptrdiff_t items_processed{0};
  for (auto chunk = 0Z; chunk < nchunks; chunk++) {
    const auto [chunkstart, chunksize] = get_range_chunk(std::ssize(dataspan), nchunks, chunk);
    assert_always(chunksize > 0);
    const auto chunk_span = dataspan.subspan(chunkstart, chunksize);
    const auto int_chunksize = static_cast<int>(chunk_span.size());
    assert_always(std::cmp_equal(int_chunksize, chunksize));
    assert_always(MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(chunk_span.data()), int_chunksize, mpi_datatype, op,
                                comm) == MPI_SUCCESS);
    items_processed += chunksize;
  }
  assert_always(items_processed == std::ssize(dataspan));
}

template <typename T>
  requires(!std::ranges::random_access_range<T>)
inline void MPI_Allreduce_safe(T& data, MPI_Op op, MPI_Comm comm) {
  MPI_Allreduce_safe(std::span{&data, 1}, op, comm);
}

template <typename R>
  requires std::ranges::random_access_range<R>
inline void MPI_Bcast_safe(R&& data, const int root, MPI_Comm comm) {
  auto dataspan = std::span{std::forward<R>(data)};
  if (dataspan.empty()) {
    return;
  }
  assert_always(dataspan.data() != nullptr);
  assert_always(comm != MPI_COMM_NULL);
  using value_t = std::ranges::range_value_t<R>;

  const auto mpi_datatype = GET_MPI_TYPE<value_t>();
  // if we're transferring bytes, then we need multiply the array count by the byte size of the type
  const ptrdiff_t sizefactor = mpi_datatype == MPI_BYTE ? sizeof(value_t) : 1;

  assert_always(MPI_COUNT_MAX > sizefactor);  // otherwise we can't make any progress
  const auto MPI_COUNT_MAX_MPITYPE = MPI_COUNT_MAX / sizefactor;
  const auto datasize_mpitype = std::ssize(dataspan) * sizefactor;
  const auto nchunks =
      (datasize_mpitype / MPI_COUNT_MAX_MPITYPE) + ((datasize_mpitype % MPI_COUNT_MAX_MPITYPE) == 0 ? 0 : 1);
  assert_always(nchunks >= 1);
  std::ptrdiff_t items_processed{0};
  for (auto chunk = 0Z; chunk < nchunks; chunk++) {
    const auto [chunkstart, chunksize] = get_range_chunk(std::ssize(dataspan), nchunks, chunk);
    assert_always(chunksize > 0);
    const auto chunk_span = dataspan.subspan(chunkstart, chunksize);
    const auto chunksize_mpitype = std::ssize(chunk_span) * sizefactor;
    assert_always(chunksize_mpitype <= MPI_COUNT_MAX);
    const auto int_chunksize_mpitype = static_cast<int>(chunksize_mpitype);
    assert_always(std::cmp_equal(int_chunksize_mpitype, chunksize_mpitype));
    assert_always(MPI_Bcast(chunk_span.data(), int_chunksize_mpitype, mpi_datatype, root, comm) == MPI_SUCCESS);
    items_processed += chunksize;
  }
  assert_always(items_processed == std::ssize(dataspan));
}

template <typename T>
  requires(!std::ranges::random_access_range<T>)
inline void MPI_Bcast_safe(T& data, const int root, MPI_Comm comm) {
  MPI_Bcast_safe(std::span{&data, 1}, root, comm);
}

template <typename R>
  requires std::ranges::random_access_range<R>
inline void MPI_Reduce_safe(R&& data, MPI_Op op, const int root, MPI_Comm comm) {
  auto dataspan = std::span{std::forward<R>(data)};
  if (dataspan.empty()) {
    return;
  }
  assert_always(dataspan.data() != nullptr);
  assert_always(op != MPI_OP_NULL);
  assert_always(comm != MPI_COMM_NULL);

  int my_rank{-1};
  assert_always(MPI_Comm_rank(comm, &my_rank) == MPI_SUCCESS);
  assert_always(my_rank >= 0);

  const auto mpi_datatype = GET_MPI_TYPE<std::ranges::range_value_t<R>>();
  assert_always(mpi_datatype != MPI_BYTE);  // we can't reduce MPI_BYTE types

  const auto nchunks = (std::ssize(dataspan) / MPI_COUNT_MAX) + ((std::ssize(dataspan) % MPI_COUNT_MAX) == 0 ? 0 : 1);
  assert_always(nchunks >= 1);
  std::ptrdiff_t items_processed{0};
  for (auto chunk = 0Z; chunk < nchunks; chunk++) {
    const auto [nstart, chunksize] = get_range_chunk(std::ssize(dataspan), nchunks, chunk);
    assert_always(chunksize > 0);
    const auto chunk_span = dataspan.subspan(nstart, chunksize);
    const auto int_chunksize = static_cast<int>(chunk_span.size());
    assert_always(std::cmp_equal(int_chunksize, chunksize));
    assert_always(MPI_Reduce(my_rank == 0 ? MPI_IN_PLACE : chunk_span.data(), chunk_span.data(), int_chunksize,
                             mpi_datatype, op, root, comm) == MPI_SUCCESS);

    items_processed += chunksize;
  }
  assert_always(items_processed == std::ssize(dataspan));
}

template <typename T>
constexpr void resize_exactly(std::vector<T>& vec, const size_t size) {
  // just resizing can (only with libstdc++?) allocate a larger capacity than needed
  vec.reserve(size);
  vec.resize(size);
}

template <double fractional_accuracy>
inline auto ftol(const double a, const double b) -> bool {
  return std::abs(a - b) <= (fractional_accuracy * std::min(std::abs(a), std::abs(b)));
}

class ScopedMutex {
 private:
  int* lock_;

  static void mutex_lock(int& lock) {
    while (std::atomic_ref<int>(lock).exchange(1, std::memory_order_acquire) == 1) {
      std::atomic_ref<int>(lock).wait(1, std::memory_order_relaxed);
      // blocks until lock != 1 (i.e., someone called unlock->notify)
    }
  }

  static void mutex_unlock(int& lock) {
    std::atomic_ref<int>(lock).store(0, std::memory_order_release);
    std::atomic_ref<int>(lock).notify_one();  // wake one sleeping thread
  }

 public:
  explicit ScopedMutex(int& lock) : lock_(&lock) { mutex_lock(*lock_); }
  ~ScopedMutex() { mutex_unlock(*lock_); }

  // disable copying and moving to avoid accidentally sharing locks between threads
  ScopedMutex(const ScopedMutex&) = delete;
  auto operator=(const ScopedMutex&) -> ScopedMutex& = delete;
  ScopedMutex(ScopedMutex&&) = delete;
  auto operator=(ScopedMutex&&) -> ScopedMutex& = delete;
};

#endif  // SN3D_H
