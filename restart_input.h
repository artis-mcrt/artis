#ifndef RESTART_INPUT_H
#define RESTART_INPUT_H

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <istream>
#include <map>
#include <sstream>
#include <string>
#include <string_view>

#include "mpi_logging.h"

namespace restart_input {
inline bool collect_files = false;
inline std::map<std::string, std::uint64_t> file_hashes;
inline std::uint64_t input_hash = 0;

// FNV-1a detects a change of the input bytes without a platform-specific library.
constexpr auto hash_text(const std::string_view text, std::uint64_t hash = 14695981039346656037ULL) -> std::uint64_t {
  for (const char value : text) {
    hash = (hash ^ static_cast<unsigned char>(value)) * 1099511628211ULL;
  }
  return hash;
}

inline void begin() {
  file_hashes.clear();
  collect_files = true;
}

inline void record_file(const std::string_view name, const std::filesystem::path& path) {
  if (!collect_files || globals::my_rank != 0 || name == "input.txt" || name == "input-newrun.txt" ||
      file_hashes.contains(std::string{name})) {
    return;
  }
  std::ifstream file(path, std::ios::binary);
  std::array<char, 65536> buffer{};
  auto hash = hash_text("");
  while (file.read(buffer.data(), buffer.size()) || file.gcount() > 0) {
    hash = hash_text(std::string_view{buffer.data(), static_cast<size_t>(file.gcount())}, hash);
  }
  if (!file.eof() || file.bad()) {
    fatal_crash("Could not read '{}' for the restart input check", path.string());
  }
  file_hashes.emplace(name, hash);
}

// Ignore comments, white space, and the rows that control the restart and the requested timestep range.
inline auto hash_parameters(std::istream& file, const std::array<int, 3>& mutable_rows) -> std::uint64_t {
  auto hash = hash_text("");
  int row = 0;
  for (std::string line; std::getline(file, line);) {
    line = line.substr(0, line.find('#'));
    std::istringstream values(line);
    std::string parameter_values;
    for (std::string value; values >> value;) {
      parameter_values.append(value).append(" ");
    }
    if (parameter_values.empty()) {
      continue;
    }
    if (row != mutable_rows[0] && row != mutable_rows[1] && row != mutable_rows[2]) {
      hash = hash_text(parameter_values, hash);
      hash = hash_text("\n", hash);
    }
    row++;
  }
  assert_always(file.eof() && !file.bad());
  return hash;
}

inline void finish() {
  collect_files = false;
  input_hash = hash_text("");
  if (globals::my_rank == 0) {
    for (const auto& [name, file_hash] : file_hashes) {
      input_hash = hash_text(name, input_hash);
      input_hash = hash_text("\n", input_hash);
      input_hash = hash_text(std::to_string(file_hash), input_hash);
      input_hash = hash_text("\n", input_hash);
    }
  }
  MPI_Bcast_safe(input_hash, 0, MPI_COMM_WORLD);
}

inline void check(const std::uint64_t saved_input_hash, const std::string_view filename) {
  if (saved_input_hash != input_hash) {
    fatal_crash("{} belongs to different input files. Restore the input files of the saved run.", filename);
  }
}
}  // namespace restart_input

#endif  // RESTART_INPUT_H
