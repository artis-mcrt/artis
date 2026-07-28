// Small shared helpers for the simulation executables: the already-running lock check and
// inline numeric utilities (exact vector resize, fractional-tolerance comparison, bin-index lookups).

#ifndef SN3D_H
#define SN3D_H

#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <print>
#include <sstream>
#include <string>
#include <vector>

#include "mpi_logging.h"

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
      pidfile.close();
      if (is_pid_running(artispid_in) && std::filesystem::current_path().generic_string() == line) {
        std::println(stderr,
                     "\n[error] artis or exspec is already running in this folder with existing pid {}. Refusing to "
                     "start. (delete artis.pid if you are sure this is incorrect)",
                     artispid_in);
        std::abort();
      }
    }

    auto pidfile = std::fstream("artis.pid", std::ofstream::out | std::ofstream::trunc);
    std::println(pidfile, "{}", artispid);
    std::println(pidfile, "{}", std::filesystem::current_path().generic_string());
  }

  // make sure rank 0 checked for a pid file before we proceed
  MPI_Barrier_allranks();
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

#endif  // SN3D_H
