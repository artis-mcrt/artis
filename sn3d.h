// Small shared helpers for the simulation executables: the already-running lock check and
// inline numeric utilities (vector sizing, tolerance comparisons, cumulative sampling, and bin-index lookups).

#ifndef SN3D_H
#define SN3D_H

#include <unistd.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <iostream>
#include <iterator>
#include <print>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "mpi_logging.h"
#include "outputfilestream.h"
#include "run_lock.h"

inline auto get_run_lock() -> RunLock& {
  static RunLock run_lock;
  return run_lock;
}

inline void acquire_run_lock() {
  if (globals::my_rank == 0) {
    std::array<char, MPI_MAX_PROCESSOR_NAME> hostname{};
    int hostname_length = 0;
    assert_always(MPI_Get_processor_name(hostname.data(), &hostname_length) == MPI_SUCCESS);
    const auto contents = std::format("{}\n{}\n{}\n", getpid(), std::filesystem::current_path().generic_string(),
                                      std::string_view{hostname.data(), static_cast<size_t>(hostname_length)});
    const auto error = get_run_lock().acquire("artis.pid", contents);
    if (error == std::errc::file_exists) {
      fatal_crash(
          "artis.pid exists. Confirm that sn3d and exspec stopped in this folder on every host. "
          "Then remove the stale file before you start a new run.");
    }
    if (error) {
      fatal_crash("Could not create artis.pid: {}", error.message());
    }
  }
  MPI_Barrier_allranks();
}

inline void release_run_lock() {
  // Every rank must finish its output before another job can use the model folder.
  MPI_Barrier_allranks();
  if (globals::my_rank == 0) {
    if (const auto error = get_run_lock().release()) {
      fatal_crash("Could not release artis.pid: {}", error.message());
    }
  }
}

template <typename T>
constexpr void reserve_resize(std::vector<T>& vec, const size_t size) {
  // just resizing can (only with libstdc++?) allocate a larger capacity than needed
  vec.reserve(size);
  vec.resize(size);
}

template <double fractional_accuracy>
inline auto ftol(const double a, const double b) -> bool {
  return std::abs(a - b) <= (fractional_accuracy * std::min(std::abs(a), std::abs(b)));
}

// Return the first index whose array value is strictly greater than target.
// Using upper_bound consistently prevents zero-weight entries from being selected when
// target is exactly equal to a repeated cumulative value.
template <typename Range, typename Value>
[[nodiscard]] constexpr auto index_upperbound(const Range& cumulative_values, const Value& target) -> ptrdiff_t {
  return std::ranges::upper_bound(cumulative_values, target) - std::begin(cumulative_values);
}

template <typename Range, typename Value>
[[nodiscard]] constexpr auto int_index_upperbound(const Range& cumulative_values, const Value& target) -> int {
  return static_cast<int>(index_upperbound(cumulative_values, target));
}

// Return the first index whose array value is greater than or equal to target.
template <typename Range, typename Value>
[[nodiscard]] constexpr auto index_lowerbound(const Range& cumulative_values, const Value& target) -> ptrdiff_t {
  return std::ranges::lower_bound(cumulative_values, target) - std::begin(cumulative_values);
}

template <typename Range, typename Value>
[[nodiscard]] constexpr auto int_index_lowerbound(const Range& cumulative_values, const Value& target) -> int {
  return static_cast<int>(index_lowerbound(cumulative_values, target));
}

static_assert(index_upperbound(std::array{1., 2., 2., 3.}, 2.) == 3);  // skips over repeated (zero-weight) entries
static_assert(index_upperbound(std::array{1., 2., 2., 3.}, 0.5) == 0);
static_assert(index_upperbound(std::array{1., 2., 2., 3.}, 3.) == 4);  // target >= all values gives the past-end index
static_assert(index_lowerbound(std::array{1., 2., 2., 3.}, 2.) == 1);
static_assert(index_lowerbound(std::array{1., 2., 2., 3.}, 4.) == 4);

// Signed bin index of value on a grid spaced uniformly by binwidth, where bin 0 starts at minvalue.
// Values below minvalue give negative indices and values past the last bin give indices >= nbins, so
// the caller must range-check the result. Flooring matters here: a plain int cast would truncate toward
// zero and misplace values up to one bin width below minvalue into bin 0.
[[nodiscard]] constexpr auto get_linearbinindex(const double value, const double minvalue, const double binwidth)
    -> ptrdiff_t {
  const double fracindex = (value - minvalue) / binwidth;
  // floor implemented with constexpr-legal operations (std::floor is not constant-evaluable on all
  // of the supported standard libraries yet)
  const auto truncated = static_cast<ptrdiff_t>(fracindex);
  return (fracindex < static_cast<double>(truncated)) ? truncated - 1 : truncated;
}

static_assert(get_linearbinindex(1.5, 1., 1.) == 0);
static_assert(get_linearbinindex(3., 1., 1.) == 2);
static_assert(get_linearbinindex(1., 1., 1.) == 0);  // bins are left-closed
static_assert(get_linearbinindex(0.5, 1., 1.) == -1);  // below the grid must not land in bin 0
static_assert(get_linearbinindex(-5., 1., 2.) == -3);

// Bin index of value on a grid of nbins bins spaced uniformly in the log of the value, where bin 0
// starts at minvalue and each bin spans dlog in log space. The caller must ensure that value is within
// the grid range (minvalue < value < maximum); the clamp only guards against floating-point rounding
// placing an in-range value just outside [0, nbins-1].
[[nodiscard]] inline auto get_logbinindex(const double value, const double minvalue, const double dlog,
                                          const ptrdiff_t nbins) -> ptrdiff_t {
  return std::clamp(static_cast<ptrdiff_t>(std::floor((std::log(value) - std::log(minvalue)) / dlog)), 0Z, nbins - 1);
}

// Edge value of a log-uniformly spaced grid whose bin 0 starts at minvalue with logarithmic spacing
// dlog, i.e. bin i spans [get_loggrid_edge(minvalue, dlog, i), get_loggrid_edge(minvalue, dlog, i + 1)).
// Inverse of get_logbinindex().
[[nodiscard]] inline auto get_loggrid_edge(const double minvalue, const double dlog, const double index) -> double {
  return std::exp(std::log(minvalue) + (index * dlog));
}

#endif  // SN3D_H
