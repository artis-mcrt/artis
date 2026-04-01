#ifndef SN3D_H
#define SN3D_H

#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

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
  MPI_Barrier_allranks();
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
