#include "mpi_logging.h"

#include <array>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <print>
#include <string_view>

namespace {

thread_local std::array<char, 1024> outputlinebuf{};
thread_local bool outputstartofline = true;

#ifdef _OPENMP
#ifndef GPU_ON
#pragma omp threadprivate(outputlinebuf, outputstartofline)
#endif
#endif

// Prepend an ISO-8601 timestamp when starting a new output line.
void print_line_start() noexcept {
  if (outputstartofline) {
    std::print(output_file, "{:%FT%TZ} ", std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now()));
  }
}
}  // anonymous namespace

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

#ifndef __NVCOMPILER_CUDA_ARCH__
void printout(const char* format, ...) noexcept {
  print_line_start();
  va_list args{};
  va_start(args, format);
  vsnprintf(outputlinebuf.data(), outputlinebuf.size(), format, args);
  va_end(args);

  const auto linebuflen = strlen(outputlinebuf.data());
  outputstartofline = (linebuflen == 0 || (outputlinebuf[linebuflen - 1] == '\n'));
  std::print(output_file, "{}", outputlinebuf.data());
  output_file.flush();
}

void report_assert_failure(const char* file, const int line, const char* expr, const char* func) noexcept {
  if (output_file) {
    std::println(output_file, "\n[rank {}] {}:{}: failed assertion `{}` in function {}", globals::my_rank, file, line,
                 expr, func);
    output_file.flush();
  }
  std::println(std::cerr, "\n[rank {}] {}:{}: failed assertion `{}` in function {}", globals::my_rank, file, line, expr,
               func);
}
#endif
