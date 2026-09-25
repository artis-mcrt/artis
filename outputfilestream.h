// Output file streams that write a plain text file or, with COMPRESS_OUTPUT_FILES, a zstd compressed
// file with the extension .zst, e.g. estimators_0000.out.zst

#ifndef OUTPUTFILESTREAM_H
#define OUTPUTFILESTREAM_H

#include <format>
#include <fstream>
#include <ios>
#include <memory>
#include <ostream>
#include <streambuf>
#include <string>
#include <string_view>
#include <utility>

#ifdef USE_ZSTD
#include <cstddef>
#include <iterator>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <zstd.h>
#pragma clang unsafe_buffer_usage end
#endif

#include "artisoptions.h"
#include "mpi_logging.h"

#ifndef USE_ZSTD
static_assert(!COMPRESS_OUTPUT_FILES, "COMPRESS_OUTPUT_FILES needs a build with libzstd. See the Makefile.");
#endif

// The level of a file that the program writes at once and then closes, e.g. a packet file. Level 13 is
// the level of exspec-after.sh, and one open stream at this level needs about 50 MB of memory.
constexpr int ZSTD_LEVEL_FILE_WRITTEN_ONCE = 13;

// The level of a file that stays open over the timesteps, e.g. an estimator file or a log. Each rank
// keeps several of them open, and one open stream at this level needs about 4 MB of memory.
constexpr int ZSTD_LEVEL_FILE_OPEN_DURING_RUN = 3;

#ifdef USE_ZSTD
// A stream buffer that compresses with zstd while it writes. A flush of the stream ends a zstd block, so
// the file on disk is readable up to the last flush, also when the program stops without an end of the
// frame. zstd -d then reports a premature end, but gives the content.
class ZstdOutputBuffer final : public std::streambuf {
 public:
  ZstdOutputBuffer(const std::string& filename, const int compression_level)
      : compressedfile(filename, std::ios::out | std::ios::trunc | std::ios::binary),
        cctx(ZSTD_createCCtx()),
        inbuf(ZSTD_CStreamInSize()),
        outbuf(ZSTD_CStreamOutSize()) {
    assert_always(cctx != nullptr);
    assert_always(ZSTD_isError(ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compression_level)) == 0U);
    setp(inbuf.data(), std::next(inbuf.data(), static_cast<std::ptrdiff_t>(inbuf.size())));
  }

  ZstdOutputBuffer(const ZstdOutputBuffer&) = delete;
  auto operator=(const ZstdOutputBuffer&) -> ZstdOutputBuffer& = delete;
  ZstdOutputBuffer(ZstdOutputBuffer&&) = delete;
  auto operator=(ZstdOutputBuffer&&) -> ZstdOutputBuffer& = delete;

  ~ZstdOutputBuffer() override {
    close();
    ZSTD_freeCCtx(cctx);
  }

  [[nodiscard]] auto is_open() const -> bool { return compressedfile.is_open(); }

  // end the zstd frame and close the file. False means that a write failed, e.g. on a full disk.
  auto close() -> bool {
    if (!compressedfile.is_open()) {
      return false;
    }
    const bool write_ok = compress_pending(ZSTD_e_end);
    compressedfile.close();
    return write_ok && !compressedfile.fail();
  }

 protected:
  auto overflow(const int_type ch) -> int_type override {
    if (!compress_pending(ZSTD_e_continue)) {
      return traits_type::eof();
    }
    if (!traits_type::eq_int_type(ch, traits_type::eof())) {
      *pptr() = traits_type::to_char_type(ch);
      pbump(1);
    }
    return traits_type::not_eof(ch);
  }

  auto sync() -> int override {
    if (!compress_pending(ZSTD_e_flush)) {
      return -1;
    }
    compressedfile.flush();
    return compressedfile.fail() ? -1 : 0;
  }

 private:
  // compress the content of the put area and write it to the file
  auto compress_pending(const ZSTD_EndDirective mode) -> bool {
    ZSTD_inBuffer input{.src = pbase(), .size = static_cast<size_t>(std::distance(pbase(), pptr())), .pos = 0};
    bool finished = false;
    while (!finished) {
      ZSTD_outBuffer output{.dst = outbuf.data(), .size = outbuf.size(), .pos = 0};
      const size_t remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
      if (ZSTD_isError(remaining) != 0U) {
        fatal_crash("zstd cannot compress the output: {}", ZSTD_getErrorName(remaining));
      }
      compressedfile.write(outbuf.data(), static_cast<std::streamsize>(output.pos));
      finished = (mode == ZSTD_e_continue) ? (input.pos == input.size) : (remaining == 0);
    }
    setp(inbuf.data(), std::next(inbuf.data(), static_cast<std::ptrdiff_t>(inbuf.size())));
    return !compressedfile.fail();
  }

  std::ofstream compressedfile;
  ZSTD_CCtx* cctx;
  std::vector<char> inbuf;
  std::vector<char> outbuf;
};
#endif

// An output stream that owns its buffer: a std::filebuf for a plain file or a ZstdOutputBuffer for a
// compressed file. The writers use it like a std::ofstream.
class OutputFileStream : public std::ostream {
 public:
  OutputFileStream() : std::ostream(nullptr) {}

  explicit OutputFileStream(std::unique_ptr<std::filebuf> filebuf_in)
      : std::ostream(filebuf_in.get()), filebuf(std::move(filebuf_in)) {}

#ifdef USE_ZSTD
  explicit OutputFileStream(std::unique_ptr<ZstdOutputBuffer> zstdbuf_in)
      : std::ostream(zstdbuf_in.get()), zstdbuf(std::move(zstdbuf_in)) {}
#endif

  OutputFileStream(const OutputFileStream&) = delete;
  auto operator=(const OutputFileStream&) -> OutputFileStream& = delete;

  // std::ostream::swap exchanges the stream state but not the buffer pointer
  OutputFileStream(OutputFileStream&& other) noexcept : std::ostream(nullptr) { take_buffers(other); }

  auto operator=(OutputFileStream&& other) noexcept -> OutputFileStream& {
    take_buffers(other);
    return *this;
  }

  ~OutputFileStream() override = default;

  [[nodiscard]] auto is_open() const -> bool {
#ifdef USE_ZSTD
    if (zstdbuf != nullptr) {
      return zstdbuf->is_open();
    }
#endif
    return filebuf != nullptr && filebuf->is_open();
  }

  // like std::ofstream::close(), a failed write or close sets the fail state
  void close() {
#ifdef USE_ZSTD
    if (zstdbuf != nullptr) {
      if (!zstdbuf->close()) {
        setstate(std::ios::failbit);
      }
      return;
    }
#endif
    if (filebuf == nullptr || filebuf->close() == nullptr) {
      setstate(std::ios::failbit);
    }
  }

 private:
  void take_buffers(OutputFileStream& other) {
    std::ostream::swap(other);
    filebuf = std::move(other.filebuf);
#ifdef USE_ZSTD
    zstdbuf = std::move(other.zstdbuf);
    if (zstdbuf != nullptr) {
      set_rdbuf(zstdbuf.get());
    } else {
      set_rdbuf(filebuf.get());
    }
#else
    set_rdbuf(filebuf.get());
#endif
    other.set_rdbuf(nullptr);
  }

  std::unique_ptr<std::filebuf> filebuf;
#ifdef USE_ZSTD
  std::unique_ptr<ZstdOutputBuffer> zstdbuf;
#endif
};

// the path of an output file: the given name, or the name with .zst when COMPRESS_OUTPUT_FILES is set
[[nodiscard]] inline auto output_filepath(const std::string_view filename) -> std::string {
  return COMPRESS_OUTPUT_FILES ? std::format("{}.zst", filename) : std::string(filename);
}

// open an output file that COMPRESS_OUTPUT_FILES does not apply to, e.g. input.txt or a restart file
[[nodiscard]] inline auto open_uncompressed_output_file(const std::string_view filename) -> OutputFileStream {
  if (filename.empty()) {
    fatal_crash("Cannot open file with empty filename.");
  }

  auto filebuf = std::make_unique<std::filebuf>();
  if (filebuf->open(std::string(filename), std::ios::out | std::ios::trunc) == nullptr) {
    fatal_crash("Could not open file '{}' for writing", filename);
  }
  return OutputFileStream(std::move(filebuf));
}

// open an output file for writing. With COMPRESS_OUTPUT_FILES, the file gets the extension .zst and the
// zstd compression at the given level.
[[nodiscard]] inline auto open_output_file(const std::string_view filename,
                                           [[maybe_unused]] const int compression_level = ZSTD_LEVEL_FILE_WRITTEN_ONCE)
    -> OutputFileStream {
#ifdef USE_ZSTD
  if constexpr (COMPRESS_OUTPUT_FILES) {
    const auto zstfilename = output_filepath(filename);
    auto zstdbuf = std::make_unique<ZstdOutputBuffer>(zstfilename, compression_level);
    if (!zstdbuf->is_open()) {
      fatal_crash("Could not open file '{}' for writing", zstfilename);
    }
    return OutputFileStream(std::move(zstdbuf));
  }
#endif
  return open_uncompressed_output_file(filename);
}

// open a per-rank output file such as estimators_0000.out in the job folder. The file stays open over
// the timesteps.
[[nodiscard]] inline auto open_rank_outfile(const std::string_view basename) -> OutputFileStream {
  return open_output_file(get_jobfolder_filepath(std::format("{}_{:04d}.out", basename, globals::my_rank)),
                          ZSTD_LEVEL_FILE_OPEN_DURING_RUN);
}

#endif  // OUTPUTFILESTREAM_H
