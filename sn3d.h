#ifndef SN3D_H
#define SN3D_H

#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <cassert>
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
                     "\nERROR: artis or exspec is already running in this folder with existing pid {}. Refusing to "
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
constexpr void resize_exactly(std::vector<T>& vec, const size_t size) {
  // just resizing can (only with libstdc++?) allocate a larger capacity than needed
  vec.reserve(size);
  vec.resize(size);
}

template <double fractional_accuracy>
inline auto ftol(const double a, const double b) -> bool {
  return std::abs(a - b) <= (fractional_accuracy * std::min(std::abs(a), std::abs(b)));
}

#endif  // SN3D_H
