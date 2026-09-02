// MPI setup (global and per-node communicators, rank and node identifiers) plus the logging
// (printlog/printlnlog) and assertion (assert_always/assert_testmodeonly) facilities used throughout the code.

#ifndef MPI_LOGGING_H
#define MPI_LOGGING_H

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <format>
#include <fstream>
#include <ios>
#include <limits>
#include <memory>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <version>

#ifdef _OPENMP
#include <omp.h>
#endif

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "constants.h"

inline void MPI_Barrier_allranks() { MPI_Barrier(MPI_COMM_WORLD); }

namespace globals {

inline MPI_Comm mpi_comm_node{MPI_COMM_NULL};
inline MPI_Comm mpi_comm_internode{MPI_COMM_NULL};

inline int nprocs{-1};
inline int my_rank{-1};

inline int node_nprocs{-1};
inline int rank_in_node{-1};

inline int node_count{-1};
inline int node_id{-1};

// optional folder (sn3d -o option) receiving the per-rank output files; empty means the current directory
inline std::string runoutputfolder;

inline void setup_mpi_vars() {
  MPI_Comm_rank(MPI_COMM_WORLD, &globals::my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &globals::nprocs);

  // make an intra-node communicator (group ranks that can share memory)
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globals::my_rank, MPI_INFO_NULL, &globals::mpi_comm_node);

  // get the local rank within this node
  MPI_Comm_rank(globals::mpi_comm_node, &globals::rank_in_node);

  // get the number of ranks on the node
  MPI_Comm_size(globals::mpi_comm_node, &globals::node_nprocs);

  MPI_Barrier_allranks();

#ifdef MAX_NODE_SIZE
  if (MAX_NODE_SIZE > 0 && globals::node_nprocs > MAX_NODE_SIZE) {
    // limit the number of ranks that can share memory
    MPI_Comm_split(globals::mpi_comm_node, globals::rank_in_node / MAX_NODE_SIZE, globals::my_rank,
                   &globals::mpi_comm_node);

    MPI_Comm_rank(globals::mpi_comm_node, &globals::rank_in_node);
    MPI_Comm_size(globals::mpi_comm_node, &globals::node_nprocs);
  }

  MPI_Barrier_allranks();
#endif

  // make an inter-node communicator (using local rank as the key for group membership)
  MPI_Comm_split(MPI_COMM_WORLD, globals::rank_in_node, globals::my_rank, &globals::mpi_comm_internode);

  // take the node id from the local rank 0 (node master) and broadcast it
  if (globals::rank_in_node == 0) {
    MPI_Comm_rank(globals::mpi_comm_internode, &globals::node_id);
    MPI_Comm_size(globals::mpi_comm_internode, &globals::node_count);
  }

  MPI_Bcast(&globals::node_id, 1, MPI_INT, 0, globals::mpi_comm_node);
  MPI_Bcast(&globals::node_count, 1, MPI_INT, 0, globals::mpi_comm_node);
}
}  // namespace globals

inline void MPI_Barrier_node() { MPI_Barrier(globals::mpi_comm_node); }

extern std::fstream output_file;

#ifdef _OPENMP
#ifndef GPU_ON
#pragma omp threadprivate(output_file)
#endif
#endif

void set_log_file(std::string_view filename) noexcept;

// Write an already-formatted message to output_file, prepending a timestamp at the start of each line. When
// add_newline is set, a trailing newline is appended and the next write starts a new line.
void log_write(std::string_view message, bool add_newline) noexcept;

// Report a failed assertion to output_file (if open) and stderr. Defined out-of-line in mpi_logging.cc so that the
// heavyweight <iostream> dependency does not propagate into every translation unit that includes this header.
[[gnu::cold]] DEVICE_FUNC void report_assert_failure(const char* file, int line, const char* expr,
                                                     const char* func) noexcept;

// Write the message with the rank, the file, the line, and the function to the rank log and to stderr, then stop
// the run. Host only. Call it through the macro fatal_crash(fmt, args...), which also has a device path.
[[noreturn]] [[gnu::cold]] void report_fatal_error_and_abort(const char* file, int line, const char* func,
                                                             const char* message) noexcept;

#define __artis_assert(e)                                                 \
  {                                                                       \
    const bool assertpass = static_cast<bool>(e);                         \
    if (!assertpass) [[unlikely]] {                                       \
      report_assert_failure(__FILE__, __LINE__, #e, __PRETTY_FUNCTION__); \
    }                                                                     \
    assert(assertpass);                                                   \
  }

template <typename... Args>
inline auto printlog(const std::format_string<Args...> fmt, Args&&... args) noexcept -> void {
  MY_IF_DEVICE(const auto str = std::vformat(fmt.get(), std::make_format_args(args...)); printf("%s", str.c_str()););
  MY_IF_HOST(log_write(std::format(fmt, std::forward<Args>(args)...), false););
}

template <typename... Args>
inline auto printlnlog(const std::format_string<Args...> fmt, Args&&... args) noexcept -> void {
  MY_IF_DEVICE(const auto str = std::vformat(fmt.get(), std::make_format_args(args...)); printf("%s\n", str.c_str()););
  MY_IF_HOST(log_write(std::format(fmt, std::forward<Args>(args)...), true););
}

// Device code has no std::format, so these helpers print a std::format string with the device printf. Each {...}
// takes the next argument and the spec inside the braces is ignored. {{ and }} print one brace.
template <typename T>
DEVICE_FUNC inline auto device_printf_arg(const T& value) -> void {
  if constexpr (std::is_same_v<T, bool>) {
    printf("%s", value ? "true" : "false");
  } else if constexpr (std::is_integral_v<T> && std::is_signed_v<T>) {
    printf("%lld", static_cast<long long>(value));
  } else if constexpr (std::is_integral_v<T>) {
    printf("%llu", static_cast<unsigned long long>(value));
  } else if constexpr (std::is_floating_point_v<T>) {
    printf("%g", static_cast<double>(value));
  } else if constexpr (std::is_convertible_v<T, const char*>) {
    printf("%s", static_cast<const char*>(value));
  } else if constexpr (requires { value.c_str(); }) {
    printf("%s", value.c_str());
  } else if constexpr (requires {
                         value.data();
                         value.size();
                       }) {
    printf("%.*s", static_cast<int>(value.size()), value.data());
  } else {
    static_assert(false, "device_printf_arg: no printf conversion for this argument type");
  }
}

// Print the part [start, end) of the format string
DEVICE_FUNC inline auto device_printf_literal(const std::string_view fmt, const size_t start, const size_t end)
    -> void {
  const auto part = fmt.substr(start, end - start);
  printf("%.*s", static_cast<int>(part.size()), part.data());
}

// Print the literal text up to and including the next {...} placeholder. Give back the position after it.
DEVICE_FUNC inline auto device_printf_until_placeholder(const std::string_view fmt, size_t pos) -> size_t {
  size_t literal_start = pos;
  while (pos < fmt.size()) {
    const char c = fmt[pos];
    const char next = (pos + 1 < fmt.size()) ? fmt[pos + 1] : '\0';
    if ((c == '{' && next == '{') || (c == '}' && next == '}')) {
      device_printf_literal(fmt, literal_start, pos + 1);
      pos += 2;
      literal_start = pos;
    } else if (c == '{') {
      device_printf_literal(fmt, literal_start, pos);
      const auto close = fmt.find('}', pos);
      return (close == std::string_view::npos) ? fmt.size() : close + 1;
    } else {
      pos++;
    }
  }
  device_printf_literal(fmt, literal_start, fmt.size());
  return pos;
}

template <typename... Args>
DEVICE_FUNC inline auto device_printf_format(const std::string_view fmt, const Args&... args) -> void {
  size_t pos = 0;
  ((pos = device_printf_until_placeholder(fmt, pos), device_printf_arg(args)), ...);
  device_printf_until_placeholder(fmt, pos);
}

// Stop the run with a message. The host formats the message with std::format and writes it with the rank, the
// file, the line, and the function to the rank log and to stderr. A device prints the same line with printf.
template <typename... Args>
[[noreturn]] DEVICE_FUNC inline auto fatal_crash_at(const char* file, const int line, const char* func,
                                                    const std::format_string<Args...> fmt, Args&&... args) noexcept
    -> void {
  MY_IF_DEVICE(printf("\n[rank %d] [error] %s:%d in %s: ", globals::my_rank, file, line, func);
               device_printf_format(fmt.get(), args...); printf("\n"); assert(false); __builtin_trap(););
  MY_IF_HOST(report_fatal_error_and_abort(file, line, func, std::format(fmt, std::forward<Args>(args)...).c_str()););
}

#define fatal_crash(...) fatal_crash_at(__FILE__, __LINE__, __PRETTY_FUNCTION__, __VA_ARGS__)

#define assert_always(e) __artis_assert(e)

#if defined TESTMODE && TESTMODE
#define assert_testmodeonly(e) __artis_assert(e)
#else
#define assert_testmodeonly(e) (void)0
#endif

// Chunk a range of integers into (approximately) equal contiguous pieces for getting around the MPI 32-bit limit
// on counts.
//
// This won't be necessary after Open MPI 6.0, which supports MPI-4's 64-bit MPI_Count functions (e.g.,
// MPI_Bcast_c instead of MPI_Bcast). For now we need this to be able to use more than ~2 billion items in a single
// array.
constexpr auto get_range_chunk(const ptrdiff_t size, const ptrdiff_t nchunks, const ptrdiff_t nchunk)
    -> std::tuple<ptrdiff_t, ptrdiff_t> {
  assert_always(size >= 0);
  assert_always(nchunks > 0);
  assert_always(nchunk >= 0);
  assert_always(nchunk < nchunks);
  const auto minchunksize = size / nchunks;  // integer division, minimum number of items per chunk
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

constexpr auto get_chunk_count(const ptrdiff_t size, const ptrdiff_t max_chunksize) -> ptrdiff_t {
  assert_always(size >= 0);
  assert_always(max_chunksize > 0);
  return (size / max_chunksize) + ((size % max_chunksize) != 0 ? 1 : 0);
}

static_assert(get_chunk_count(0, 3) == 0);
static_assert(get_chunk_count(1, 3) == 1);
static_assert(get_chunk_count(3, 3) == 1);
static_assert(get_chunk_count(4, 3) == 2);

template <typename T>
  requires(!std::is_const_v<T>)
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
  assert_always(MPI_Info_create(&info) == MPI_SUCCESS);
  // Request alignment (AVX-512 requires 64b, and 128b is Apple Silicon cache line size).
  assert_always(MPI_Info_set(info, "mpi_minimum_memory_alignment", "128") == MPI_SUCCESS);
  assert_always(MPI_Win_allocate_shared(size, disp_unit, info, globals::mpi_comm_node, static_cast<void*>(&ptr),
                                        &mpiwin) == MPI_SUCCESS);
  assert_always(MPI_Info_free(&info) == MPI_SUCCESS);
  assert_always(MPI_Win_shared_query(mpiwin, 0, &size, &disp_unit, static_cast<void*>(&ptr)) == MPI_SUCCESS);
  assert_always(ptr != nullptr);
#ifdef __cpp_lib_is_sufficiently_aligned
  assert_always(std::is_sufficiently_aligned<128>(ptr));
#endif
#pragma clang unsafe_buffer_usage begin
  const auto newspan = std::span<T>(ptr, num_allranks);
#pragma clang unsafe_buffer_usage end
  // initialise the shared memory. Each rank on the node fills its own contiguous slice so that,
  // under the operating system's first-touch page placement policy, the physical pages are
  // distributed across the NUMA domains spanned by the node's ranks instead of all being homed
  // on the domain of rank 0 (relevant for multi-CCD/multi-socket x86; no effect on UMA systems)
  const auto [slicestart, slicecount] = get_range_chunk(num_allranks, globals::node_nprocs, globals::rank_in_node);
  std::ranges::fill(newspan.subspan(slicestart, slicecount), initval);
  MPI_Barrier(globals::mpi_comm_node);
  return {newspan, mpiwin};
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

template <typename T>
  requires std::is_trivially_copyable_v<T>
class MPI_shared_array {
  friend class MPI_shared_array<const T>;  // allow conversion from non-const to const version
  friend class MPI_shared_array<std::remove_const_t<T>>;

 private:
  MPI_Win _win{MPI_WIN_NULL};
  std::span<T> _span{};

 public:
  MPI_shared_array() = default;

  explicit MPI_shared_array(const ptrdiff_t num_allranks, const T& initval = {}) { allocate(num_allranks, initval); }

  // copy constructor is deleted to avoid multiple owners of the same MPI window, but move constructor is allowed
  MPI_shared_array(const MPI_shared_array&) = delete;
  MPI_shared_array(MPI_shared_array&& other) noexcept : _win(other._win), _span(other._span) {
    // prevent the other object from freeing the window in its destructor
    other._win = MPI_WIN_NULL;
    other._span = {};
  }

  template <typename U>
    requires(std::is_same_v<T, const U> && !std::is_const_v<U>)
  // NOLINTNEXTLINE(cppcoreguidelines-rvalue-reference-param-not-moved,*-explicit-constructor,hicpp-explicit-conversions)
  MPI_shared_array(MPI_shared_array<U>&& other) noexcept  // cppcheck-suppress noExplicitConstructor
      : _win(std::exchange(other._win, MPI_WIN_NULL)), _span(std::exchange(other._span, {})) {}

  auto operator=(const MPI_shared_array<T>&) -> MPI_shared_array& = delete;

  auto operator=(MPI_shared_array&& other) noexcept -> MPI_shared_array& {
    if (this->_span.data() == other._span.data()) {
      return *this;
    }
    assert_always(_span.empty() && (_win == MPI_WIN_NULL));
    _span = std::exchange(other._span, {});
    _win = std::exchange(other._win, MPI_WIN_NULL);
    return *this;
  }

  template <typename U>
    requires(std::is_same_v<T, const U> && !std::is_const_v<U>)
  auto operator=(MPI_shared_array<U>&& other_) noexcept -> MPI_shared_array& {
    if (this->_span.data() == other_._span.data()) {
      return *this;
    }
    assert_always(_span.empty() && (_win == MPI_WIN_NULL));
    auto other = std::move(other_);
    _span = static_cast<std::span<T>>(std::exchange(other._span, {}));
    _win = std::exchange(other._win, MPI_WIN_NULL);
    return *this;
  }

  ~MPI_shared_array() { reset(); }

  auto allocate(const ptrdiff_t num_allranks, const T& initval = {}) {
    assert_always(_span.empty() && (_win == MPI_WIN_NULL));  // should not be allocating if we already own a window
    if (globals::node_nprocs > 1) {
      int initialized = 0;
      MPI_Initialized(&initialized);
      assert_always(initialized != 0);  // MPI must be initialized before constructing an MPI_shared_array
      std::tie(_span, _win) = MPI_shared_malloc_span_keepwin<T>(num_allranks, initval);
    } else {
#pragma clang unsafe_buffer_usage begin
      _span = std::span<T>(new T[num_allranks], num_allranks);
#pragma clang unsafe_buffer_usage end
      std::ranges::fill(_span, initval);
    }
  }

  auto reset() {
    if constexpr (TESTMODE) {
      printlnlog("freeing MPI_shared_array of size {}", _span.size());
    }
    if (_win != MPI_WIN_NULL) {
      int finalized = 0;
      MPI_Finalized(&finalized);
      // do not attempt to free MPI windows after MPI_Finalize
      if (finalized == 0) {
        MPI_Win_free(&_win);
      }
      _win = MPI_WIN_NULL;
    } else {
      delete[] _span.data();
    }
    _span = {};
  }

  // Conversion to a const span is allowed on const objects.
  explicit operator std::span<const T>() const { return _span; }

  // mutable span if T is not const
  template <typename U = T>
    requires(!std::is_const_v<U>)
  explicit operator std::span<U>() {
    return _span;
  }
  // Mutable span accessor.
  [[nodiscard]] auto span() -> std::span<T> { return _span; }  // cppcheck-suppress functionConst
  // Read-only span accessor.
  [[nodiscard]] auto span() const -> std::span<const T> { return std::span<const T>{_span}; }
  // Mutable data pointer.
  [[nodiscard]] auto data() -> T* { return _span.data(); }
  // Read-only data pointer.
  [[nodiscard]] auto data() const -> const T* { return _span.data(); }
  // Iterators for mutable access.
  [[nodiscard]] auto begin() { return _span.begin(); }
  [[nodiscard]] auto end() { return _span.end(); }
  // Iterators for read-only access.
  [[nodiscard]] auto begin() const { return std::span<const T>{_span}.begin(); }
  [[nodiscard]] auto end() const { return std::span<const T>{_span}.end(); }
  [[nodiscard]] auto empty() const -> bool { return _span.empty(); }
  // Mutable subspan accessor.
  [[nodiscard]] auto subspan(const size_t offset, const size_t count) -> std::span<T> {
    return _span.subspan(offset, count);
  }
  // Read-only subspan accessor.
  [[nodiscard]] auto subspan(const size_t offset, const size_t count) const -> std::span<const T> {
    return std::span<const T>{_span}.subspan(offset, count);
  }
  [[nodiscard]] auto first(const size_t count) -> std::span<T> { return _span.first(count); }
  [[nodiscard]] auto first(const size_t count) const -> std::span<const T> {
    return std::span<const T>{_span}.first(count);
  }
  [[nodiscard]] auto size() const -> size_t { return _span.size(); }
  // (std::span has no ssize() member, so compute the signed size from size())
  [[nodiscard]] auto ssize() const -> ptrdiff_t { return static_cast<ptrdiff_t>(_span.size()); }

  // define operator[] to allow direct indexing into the span
  auto operator[](const size_t index) noexcept -> T& { return _span[index]; }
  auto operator[](const size_t index) const noexcept -> const T& { return std::span<const T>{_span}[index]; }
};

// MPI operations use a 32-bit int for the count, so we need to chunk large arrays
constexpr std::ptrdiff_t MPI_COUNT_MAX = std::numeric_limits<int>::max();

// these wrappers add type, bounds, and overflow safety to the MPI calls
template <typename R>
  requires std::ranges::random_access_range<R>
inline void MPI_Allreduce_safe(R&& data, MPI_Op op, MPI_Comm comm) {
  const auto dataspan = std::span{std::forward<R>(data)};
  if (dataspan.empty()) {
    return;
  }
  assert_always(dataspan.data() != nullptr);
  assert_always(op != MPI_OP_NULL);
  assert_always(comm != MPI_COMM_NULL);

  const auto mpi_datatype = GET_MPI_TYPE<std::ranges::range_value_t<R>>();
  assert_always(mpi_datatype != MPI_BYTE);  // we can't reduce MPI_BYTE types

  const auto nchunks = get_chunk_count(std::ssize(dataspan), MPI_COUNT_MAX);
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
  const auto dataspan = std::span{std::forward<R>(data)};
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
  // maximum number of items per chunk such that the MPI count (items * sizefactor) stays within MPI_COUNT_MAX
  const auto max_items_per_chunk = MPI_COUNT_MAX / sizefactor;
  const auto nitems = std::ssize(dataspan);
  const auto nchunks = get_chunk_count(nitems, max_items_per_chunk);
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
  const auto dataspan = std::span{std::forward<R>(data)};
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

  const auto nchunks = get_chunk_count(std::ssize(dataspan), MPI_COUNT_MAX);
  assert_always(nchunks >= 1);
  std::ptrdiff_t items_processed{0};
  for (auto chunk = 0Z; chunk < nchunks; chunk++) {
    const auto [nstart, chunksize] = get_range_chunk(std::ssize(dataspan), nchunks, chunk);
    assert_always(chunksize > 0);
    const auto chunk_span = dataspan.subspan(nstart, chunksize);
    const auto int_chunksize = static_cast<int>(chunk_span.size());
    assert_always(std::cmp_equal(int_chunksize, chunksize));
    // MPI_IN_PLACE as the send buffer is only valid at the root rank; all other ranks must pass
    // the data as the send buffer (their receive buffer is ignored)
    assert_always(MPI_Reduce(my_rank == root ? MPI_IN_PLACE : chunk_span.data(), chunk_span.data(), int_chunksize,
                             mpi_datatype, op, root, comm) == MPI_SUCCESS);

    items_processed += chunksize;
  }
  assert_always(items_processed == std::ssize(dataspan));
}

// path for a per-rank output file (rank logs, estimators, nlte/radfield/macroatom files), which the
// sn3d -o option redirects into a run output folder (stored without a trailing slash)
[[nodiscard]] inline auto get_runoutputfolder_filepath(const std::string_view filename) -> std::string {
  return globals::runoutputfolder.empty() ? std::string(filename)
                                          : std::format("{}/{}", globals::runoutputfolder, filename);
}

// exactly match the generated per-rank output filenames: output_<rank>-<thread>.txt and the
// estimators/nlte/radfield/macroatom _<rank>.out files, possibly with a compression extension added by
// the post-processing scripts (e.g. exspec-after.sh runs zstd)
[[nodiscard]] inline auto is_rank_outfile_name(std::string_view filename) -> bool {
  const auto alldigits = [](const std::string_view str) {
    return !str.empty() && std::ranges::all_of(str, [](const char c) { return c >= '0' && c <= '9'; });
  };

  for (const std::string_view compressext : {".zst", ".gz", ".xz"}) {
    if (filename.ends_with(compressext)) {
      filename.remove_suffix(compressext.size());
      break;
    }
  }

  if (constexpr std::string_view logprefix = "output_"; filename.starts_with(logprefix) && filename.ends_with(".txt")) {
    const auto rank_thread = filename.substr(logprefix.size(), filename.size() - logprefix.size() - 4);
    const auto dashpos = rank_thread.find('-');
    return dashpos != std::string_view::npos && alldigits(rank_thread.substr(0, dashpos)) &&
           alldigits(rank_thread.substr(dashpos + 1));
  }

  if (filename.ends_with(".out")) {
    for (const std::string_view prefix : {"estimators_", "nlte_", "radfield_", "macroatom_"}) {
      if (filename.starts_with(prefix)) {
        return alldigits(filename.substr(prefix.size(), filename.size() - prefix.size() - 4));
      }
    }
  }

  return false;
}

[[nodiscard]] inline auto fopen_required(const std::string& filename, std::span<const char> mode) -> FILE* {
  if (mode[0] == 'r') {
    // search data folders in order to find file to read
    for (const auto& datadir : datafolders) {
      const auto datafolderfilename = std::format("{}{}", datadir, filename);
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

  fatal_crash("Could not open file '{}' for mode '{}'.", filename, mode.data());
}

[[nodiscard]] inline auto fopen_required_uniqueptr(const std::string& filename, std::span<const char> mode) {
  return std::unique_ptr<FILE, int (*)(FILE*)>(fopen_required(filename, mode),
                                               [](FILE* fp) -> int { return std::fclose(fp); });
}

[[nodiscard]] inline auto fstream_required(const std::string_view filename, std::ios::openmode mode) -> std::fstream {
  if (filename.empty()) {
    fatal_crash("Cannot open file with empty filename.");
  }

  if ((mode & std::ios::in) != 0U) {
    // search data folders in order to find file to read
    for (const auto& datadir : datafolders) {
      const auto datafolderfilename = std::format("{}{}", datadir, filename);
      auto file = std::fstream(datafolderfilename, mode);
      if (file.is_open()) {
        return file;
      }
    }
  } else {
    // don't prepend data folders when writing
    auto file = std::fstream(std::string(filename), mode);
    if (file.is_open()) {
      return file;
    }
  }

  fatal_crash("Could not open file '{}'", filename);
}

// open a per-rank output file such as estimators_0000.out for writing
[[nodiscard]] inline auto open_rank_outfile(const std::string_view basename) -> std::fstream {
  return fstream_required(get_runoutputfolder_filepath(std::format("{}_{:04d}.out", basename, globals::my_rank)),
                          std::ios::out | std::ios::trunc);
}

// padded to a full cache line in CPU multithreaded modes so that adjacent mutexes in an array don't false share
class ALIGNAS_AVOID_FALSE_SHARING PaddedMutex {
 private:
  int lock_{0};

  friend class ScopedMutex;

 public:
  constexpr PaddedMutex() = default;
  constexpr explicit PaddedMutex(const int lock) : lock_(lock) {}

  constexpr auto operator=(const int lock) noexcept -> PaddedMutex& {
    lock_ = lock;
    return *this;
  }
};

static_assert(!align_to_avoid_false_sharing || alignof(PaddedMutex) >= DESTRUCTIVE_INTERFERENCE_SIZE);
static_assert(!align_to_avoid_false_sharing || alignof(PaddedMutex) % DESTRUCTIVE_INTERFERENCE_SIZE == 0);
static_assert(!align_to_avoid_false_sharing || sizeof(PaddedMutex) >= DESTRUCTIVE_INTERFERENCE_SIZE);
static_assert(!align_to_avoid_false_sharing || sizeof(PaddedMutex) % DESTRUCTIVE_INTERFERENCE_SIZE == 0);

class ScopedMutex {
 private:
  PaddedMutex* lock_;

  static void mutex_lock(PaddedMutex& lock) {
    while (std::atomic_ref<int>(lock.lock_).exchange(1, std::memory_order_acquire) == 1) {
      std::atomic_ref<int>(lock.lock_).wait(1, std::memory_order_relaxed);
      // blocks until lock != 1 (i.e., someone called unlock->notify)
    }
  }

  static void mutex_unlock(PaddedMutex& lock) {
    std::atomic_ref<int>(lock.lock_).store(0, std::memory_order_release);
    std::atomic_ref<int>(lock.lock_).notify_one();  // wake one sleeping thread
  }

 public:
  explicit ScopedMutex(PaddedMutex& lock) : lock_(&lock) { mutex_lock(*lock_); }
  ~ScopedMutex() { mutex_unlock(*lock_); }

  // disable copying and moving to avoid accidentally sharing locks between threads
  ScopedMutex(const ScopedMutex&) = delete;
  auto operator=(const ScopedMutex&) -> ScopedMutex& = delete;
  ScopedMutex(ScopedMutex&&) = delete;
  auto operator=(ScopedMutex&&) -> ScopedMutex& = delete;
};

#endif  // MPI_LOGGING_H
