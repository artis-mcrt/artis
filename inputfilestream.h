// Input file streams that read a plain text file or, when libzstd is linked, a zstd compressed
// file with the same name and the extension .zst, e.g. model.txt.zst

#ifndef INPUTFILESTREAM_H
#define INPUTFILESTREAM_H

#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <istream>
#include <memory>
#include <streambuf>
#include <string_view>
#include <utility>

#ifdef USE_ZSTD
#include <iterator>
#include <string>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <zstd.h>
#pragma clang unsafe_buffer_usage end
#endif

#include "constants.h"
#include "mpi_logging.h"

#ifdef USE_ZSTD
// A stream buffer that decompresses a zstd file while it reads. A seek to an earlier position
// starts the decompression again from the start of the file. The readers seek back only over the
// header lines of a file, so the repeated work is small.
class ZstdInputBuffer final : public std::streambuf {
 public:
  explicit ZstdInputBuffer(std::string filename_in)
      : filename(std::move(filename_in)),
        compressedfile(filename, std::ios::in | std::ios::binary),
        dctx(ZSTD_createDCtx()),
        compressedbuf(ZSTD_DStreamInSize()),
        decompressedbuf(ZSTD_DStreamOutSize()) {
    assert_always(dctx != nullptr);
  }

  ZstdInputBuffer(const ZstdInputBuffer&) = delete;
  auto operator=(const ZstdInputBuffer&) -> ZstdInputBuffer& = delete;
  ZstdInputBuffer(ZstdInputBuffer&&) = delete;
  auto operator=(ZstdInputBuffer&&) -> ZstdInputBuffer& = delete;

  ~ZstdInputBuffer() override { ZSTD_freeDCtx(dctx); }

  [[nodiscard]] auto is_open() const -> bool { return compressedfile.is_open(); }

 protected:
  auto underflow() -> int_type override {
    if (gptr() < egptr()) {
      return traits_type::to_int_type(*gptr());
    }

    decompressedpos_bufstart += std::distance(eback(), egptr());

    ZSTD_outBuffer output{.dst = decompressedbuf.data(), .size = decompressedbuf.size(), .pos = 0};
    while (output.pos == 0) {
      if (compressedinput.pos == compressedinput.size) {
        compressedfile.read(compressedbuf.data(), static_cast<std::streamsize>(compressedbuf.size()));
        const auto bytesread = compressedfile.gcount();
        if (bytesread <= 0) {
          if (frame_incomplete) {
            fatal_crash("{} ends inside a zstd frame. The file is truncated.", filename);
          }
          setg(nullptr, nullptr, nullptr);
          return traits_type::eof();
        }
        compressedinput = {.src = compressedbuf.data(), .size = static_cast<size_t>(bytesread), .pos = 0};
      }

      const size_t decompress_result = ZSTD_decompressStream(dctx, &output, &compressedinput);
      if (ZSTD_isError(decompress_result) != 0U) {
        fatal_crash("{}: zstd cannot decompress the file: {}", filename, ZSTD_getErrorName(decompress_result));
      }
      // a result of zero means that a frame ended. The next call starts the next frame, if there is one.
      frame_incomplete = (decompress_result != 0);
    }

    setg(decompressedbuf.data(), decompressedbuf.data(),
         std::next(decompressedbuf.data(), static_cast<std::ptrdiff_t>(output.pos)));
    return traits_type::to_int_type(*gptr());
  }

  auto seekoff(const off_type off, const std::ios_base::seekdir dir, const std::ios_base::openmode which)
      -> pos_type override {
    if (dir == std::ios_base::cur) {
      return seekpos(decompressedpos_bufstart + std::distance(eback(), gptr()) + off, which);
    }
    if (dir == std::ios_base::beg) {
      return seekpos(off, which);
    }
    // the decompressed size is unknown, so a seek from the end is not possible
    return {static_cast<off_type>(-1)};
  }

  auto seekpos(const pos_type pos, const std::ios_base::openmode which) -> pos_type override {
    const off_type target = pos;
    if ((which & std::ios_base::in) == 0 || target < 0) {
      return {static_cast<off_type>(-1)};
    }

    if (target < decompressedpos_bufstart) {
      // start again from the start of the file
      ZSTD_DCtx_reset(dctx, ZSTD_reset_session_only);
      compressedfile.clear();
      compressedfile.seekg(0);
      compressedinput = {};
      decompressedpos_bufstart = 0;
      frame_incomplete = false;
      setg(nullptr, nullptr, nullptr);
    }

    while (true) {
      const auto buflen = std::distance(eback(), egptr());
      if (target <= decompressedpos_bufstart + buflen) {
        setg(eback(), std::next(eback(), static_cast<std::ptrdiff_t>(target - decompressedpos_bufstart)), egptr());
        return {target};
      }
      // use up the buffer, so that underflow() decompresses the next part
      gbump(static_cast<int>(buflen - std::distance(eback(), gptr())));
      if (traits_type::eq_int_type(underflow(), traits_type::eof())) {
        return {static_cast<off_type>(-1)};
      }
    }
  }

 private:
  std::string filename;
  std::ifstream compressedfile;
  ZSTD_DCtx* dctx;
  std::vector<char> compressedbuf;
  std::vector<char> decompressedbuf;
  ZSTD_inBuffer compressedinput{};
  std::streamoff decompressedpos_bufstart = 0;  // decompressed byte offset of eback()
  bool frame_incomplete = false;
};
#endif

// An input stream that owns its buffer: a std::filebuf for a plain file or a ZstdInputBuffer for
// a compressed file. The readers use it like any std::istream.
class InputFileStream : public std::istream {
 public:
  explicit InputFileStream(std::unique_ptr<std::streambuf> buffer_in)
      : std::istream(buffer_in.get()), buffer(std::move(buffer_in)) {}

  InputFileStream(const InputFileStream&) = delete;
  auto operator=(const InputFileStream&) -> InputFileStream& = delete;

  // std::istream::swap exchanges the stream state but not the buffer pointer
  InputFileStream(InputFileStream&& other) noexcept : std::istream(nullptr), buffer(std::move(other.buffer)) {
    std::istream::swap(other);
    set_rdbuf(buffer.get());
    other.set_rdbuf(nullptr);
  }

  auto operator=(InputFileStream&& other) noexcept -> InputFileStream& {
    std::istream::swap(other);
    buffer = std::move(other.buffer);
    set_rdbuf(buffer.get());
    other.set_rdbuf(nullptr);
    return *this;
  }

  ~InputFileStream() override = default;

 private:
  std::unique_ptr<std::streambuf> buffer;
};

// Open an input file for reading. The search covers the folders of datafolders in order. In each
// folder, the plain name comes before the compressed name, e.g. model.txt before model.txt.zst.
[[nodiscard]] inline auto istream_required(const std::string_view filename) -> InputFileStream {
  if (filename.empty()) {
    fatal_crash("Cannot open file with empty filename.");
  }

  for (const auto& datadir : datafolders) {
    const auto datafolderfilename = std::format("{}{}", datadir, filename);
    auto filebuf = std::make_unique<std::filebuf>();
    if (filebuf->open(datafolderfilename, std::ios::in) != nullptr) {
      return InputFileStream(std::move(filebuf));
    }

    const auto zstfilename = std::format("{}.zst", datafolderfilename);
#ifdef USE_ZSTD
    auto zstdbuf = std::make_unique<ZstdInputBuffer>(zstfilename);
    if (zstdbuf->is_open()) {
      return InputFileStream(std::move(zstdbuf));
    }
#else
    if (std::filesystem::exists(zstfilename)) {
      fatal_crash(
          "Found '{}', but this build has no libzstd. Build with libzstd (see the Makefile) or decompress the file.",
          zstfilename);
    }
#endif
  }

  fatal_crash("Could not open file '{}'", filename);
}

// True if the run folder holds the input file, in plain or in compressed form
[[nodiscard]] inline auto inputfile_exists(const std::string_view filename) -> bool {
  return std::filesystem::exists(std::filesystem::path(filename)) ||
         std::filesystem::exists(std::filesystem::path(std::format("{}.zst", filename)));
}

#endif  // INPUTFILESTREAM_H
