// Implementation of the log-file output routines declared in mpi_logging.h.

#include "mpi_logging.h"

#include <chrono>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <print>
#include <string_view>

#include "constants.h"
#include "outputfilestream.h"

namespace {

// GCC counts the definition of a variable with a constructor as a use, so the threadprivate pragma must come
// between this declaration and the definition below.
extern OutputFileStream output_file;
bool outputstartofline = true;

#ifdef _OPENMP
#ifndef GPU_ON
#pragma omp threadprivate(output_file, outputstartofline)
#endif
#endif

OutputFileStream output_file;

// Prepend an ISO-8601 timestamp when starting a new output line.
void print_line_start() noexcept {
  if (outputstartofline) {
    std::print(output_file, "{:%FT%TZ} ", std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now()));
  }
}
}  // anonymous namespace

// the logs stay plain, so tail -f and the job scripts read them. exspec-after.sh compresses them later.
void set_log_file(const std::string_view filename) noexcept { output_file = open_uncompressed_output_file(filename); }

void log_write(const std::string_view message, const bool add_newline) noexcept {
  print_line_start();
  if (add_newline) {
    std::println(output_file, "{}", message);
    outputstartofline = true;
  } else {
    std::print(output_file, "{}", message);
    if (!message.empty()) {
      outputstartofline = (message.back() == '\n');
    }
  }
  output_file.flush();
}

DEVICE_FUNC void report_assert_failure(const char* file, const int line, const char* expr, const char* func) noexcept {
  MY_IF_DEVICE(
      printf("\n[rank %d] %s:%d: failed assertion `%s` in function %s\n", globals::my_rank, file, line, expr, func););
  MY_IF_HOST(if (output_file) {
    std::println(output_file, "\n[rank {}] {}:{}: failed assertion `{}` in function {}", globals::my_rank, file, line,
                 expr, func);
    output_file.flush();
  } std::println(std::cerr, "\n[rank {}] {}:{}: failed assertion `{}` in function {}", globals::my_rank, file, line,
                 expr, func););
}

void report_fatal_error_and_abort(const char* file, const int line, const char* func, const char* message) noexcept {
  if (output_file) {
    std::println(output_file, "\n[rank {}] [error] {}:{} in {}: {}", globals::my_rank, file, line, func, message);
    output_file.flush();
  }
  std::println(std::cerr, "\n[rank {}] [error] {}:{} in {}: {}", globals::my_rank, file, line, func, message);
  std::abort();
}
