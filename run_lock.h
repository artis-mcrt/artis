// An exclusive file lock for a simulation or spectrum calculation in a model folder.

#ifndef RUN_LOCK_H
#define RUN_LOCK_H

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <cstddef>
#include <filesystem>
#include <span>
#include <string_view>
#include <system_error>

class RunLock {
  int file_descriptor = -1;
  std::filesystem::path lock_path;

 public:
  RunLock() = default;
  RunLock(const RunLock&) = delete;
  auto operator=(const RunLock&) -> RunLock& = delete;
  RunLock(RunLock&&) = delete;
  auto operator=(RunLock&&) -> RunLock& = delete;

  ~RunLock() {
    // Keep the file after an error. A normal exit calls release() after all ranks finish their output.
    if (file_descriptor >= 0) {
      close(file_descriptor);
    }
  }

  [[nodiscard]] auto acquire(const std::filesystem::path& path, const std::string_view contents) -> std::error_code {
    if (file_descriptor >= 0) {
      return std::make_error_code(std::errc::operation_in_progress);
    }
    lock_path = std::filesystem::absolute(path);
    // The existence check and file creation form one operation, also for other hosts on the same filesystem.
    file_descriptor = open(lock_path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (file_descriptor < 0) {
      return {errno, std::generic_category()};
    }

    auto unwritten = std::span(contents);
    while (!unwritten.empty()) {
      const auto count = write(file_descriptor, unwritten.data(), unwritten.size());
      if (count < 0) {
        if (errno == EINTR) {
          continue;
        }
        return {errno, std::generic_category()};
      }
      if (count == 0) {
        return std::make_error_code(std::errc::io_error);
      }
      unwritten = unwritten.subspan(static_cast<size_t>(count));
    }
    return {};
  }

  [[nodiscard]] auto release() noexcept -> std::error_code {
    if (file_descriptor < 0) {
      return {};
    }
    struct stat owned_file{};
    struct stat named_file{};
    if (fstat(file_descriptor, &owned_file) != 0) {
      return {errno, std::generic_category()};
    }
    if (lstat(lock_path.c_str(), &named_file) != 0) {
      return {errno, std::generic_category()};
    }
    // Keep the descriptor open until removal. This prevents reuse of the inode while we own the lock.
    if (owned_file.st_dev != named_file.st_dev || owned_file.st_ino != named_file.st_ino) {
      return std::make_error_code(std::errc::state_not_recoverable);
    }
    if (unlink(lock_path.c_str()) != 0) {
      return {errno, std::generic_category()};
    }
    const int released_descriptor = file_descriptor;
    file_descriptor = -1;
    if (close(released_descriptor) != 0) {
      return {errno, std::generic_category()};
    }
    return {};
  }
};

// The lock of artis.pid that this process owns. The fatal error path of mpi_logging.cc also releases it.
inline RunLock run_lock;

#endif  // RUN_LOCK_H
