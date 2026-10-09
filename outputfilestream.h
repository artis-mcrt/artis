// Output file streams. A build with libzstd writes each output file zstd compressed with the extension
// .zst, e.g. estimators_0000.out.zst. A build without libzstd writes plain text files.

#ifndef OUTPUTFILESTREAM_H
#define OUTPUTFILESTREAM_H

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <memory>
#include <ostream>
#include <span>
#include <sstream>
#include <streambuf>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#ifdef USE_ZSTD
#include <iterator>

#pragma clang unsafe_buffer_usage begin
#include <zstd.h>
#pragma clang unsafe_buffer_usage end
#endif

#include "artisoptions.h"
#include "mpi_logging.h"

// The zstd level of the output files. Level 9 gives files that are 2 to 3 percent larger than level 13,
// at 6 times the speed. One open stream at this level needs about 16 MB of memory, and a rank keeps up
// to four streams open over the timesteps: the estimator, nlte, radfield, and macroatom files.
constexpr int ZSTD_LEVEL_DEFAULT = 9;

#ifdef USE_ZSTD
// A stream buffer that compresses with zstd while it writes. A flush of the stream ends a zstd frame,
// so every reader gets the content up to the last flush, also when the program stops without a close.
// Each frame carries a checksum of its content. With worker threads, zstd compresses the blocks of a
// large file in parallel. A libzstd without thread support ignores the request.
class ZstdOutputBuffer final : public std::streambuf {
 public:
  ZstdOutputBuffer(const std::string& filename, const int compression_level, const int worker_threads)
      : compressedfile(filename, std::ios::out | std::ios::trunc | std::ios::binary),
        cctx(ZSTD_createCCtx()),
        inbuf(ZSTD_CStreamInSize()),
        outbuf(ZSTD_CStreamOutSize()) {
    assert_always(cctx != nullptr);
    assert_always(ZSTD_isError(ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compression_level)) == 0U);
    assert_always(ZSTD_isError(ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1)) == 0U);
    if (worker_threads > 0) {
      ZSTD_CCtx_setParameter(cctx, ZSTD_c_nbWorkers, worker_threads);
    }
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
    // a file with no content gets one empty frame, so that every reader accepts it
    const bool write_ok = end_frame(!wrote_frame);
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
    if (!end_frame(false)) {
      return -1;
    }
    compressedfile.flush();
    return compressedfile.fail() ? -1 : 0;
  }

 private:
  // end the current frame, if the stream got content since the last frame end or if the caller asks
  // for an empty frame
  auto end_frame(const bool also_when_empty) -> bool {
    if (!also_when_empty && !frame_open && pptr() == pbase()) {
      return true;
    }
    return compress_pending(ZSTD_e_end);
  }

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
    frame_open = (mode != ZSTD_e_end);
    wrote_frame = wrote_frame || (mode == ZSTD_e_end);
    setp(inbuf.data(), std::next(inbuf.data(), static_cast<std::ptrdiff_t>(inbuf.size())));
    return !compressedfile.fail();
  }

  std::ofstream compressedfile;
  ZSTD_CCtx* cctx;
  std::vector<char> inbuf;
  std::vector<char> outbuf;
  bool frame_open = false;
  bool wrote_frame = false;
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
  OutputFileStream(OutputFileStream&&) = delete;

  // std::ostream::swap exchanges the stream state but not the buffer pointer. The moved-from stream
  // has no buffer, so it gets the bad state.
  auto operator=(OutputFileStream&& other) noexcept -> OutputFileStream& {
    if (this == &other) {
      return *this;
    }
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
    other.setstate(std::ios::badbit);
    return *this;
  }

  ~OutputFileStream() override = default;

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
  std::unique_ptr<std::filebuf> filebuf;
#ifdef USE_ZSTD
  std::unique_ptr<ZstdOutputBuffer> zstdbuf;
#endif
};

// the path of an output file: the name with .zst in a build with libzstd, else the given name
[[nodiscard]] inline auto output_filepath(const std::string_view filename) -> std::string {
#ifdef USE_ZSTD
  return std::format("{}.zst", filename);
#else
  return std::string(filename);
#endif
}

// Remove the file of the other form, plain or compressed. A reader opens the plain name first, so a
// stale file of an earlier build with the other form must not stay next to the new file.
inline void remove_other_output_form(const std::string_view filename) {
#ifdef USE_ZSTD
  const auto otherpath = std::filesystem::path(filename);
#else
  const auto otherpath = std::filesystem::path(std::format("{}.zst", filename));
#endif
  std::error_code ec;
  std::filesystem::remove(otherpath, ec);
}

// open an output file that stays plain in every build, e.g. input.txt, a log, or a restart file
[[nodiscard]] inline auto open_uncompressed_output_file(const std::string_view filename) -> OutputFileStream {
  if (filename.empty()) {
    fatal_crash("Cannot open file with empty filename.");
  }

  auto filebuf = std::make_unique<std::filebuf>();
  if (filebuf->open(std::string(filename), std::ios::out | std::ios::trunc) == nullptr) {
    fatal_crash("Could not open the output file '{}'", filename);
  }
  return OutputFileStream(std::move(filebuf));
}

// Open an output file. In a build with libzstd, the file gets the extension .zst and the zstd
// compression at the given level, with the given number of worker threads.
[[nodiscard]] inline auto open_output_file(const std::string_view filename,
                                           [[maybe_unused]] const int compression_level = ZSTD_LEVEL_DEFAULT,
                                           [[maybe_unused]] const int worker_threads = 0) -> OutputFileStream {
  if (filename.empty()) {
    fatal_crash("Cannot open file with empty filename.");
  }
  remove_other_output_form(filename);
#ifdef USE_ZSTD
  const auto zstfilename = output_filepath(filename);
  auto zstdbuf = std::make_unique<ZstdOutputBuffer>(zstfilename, compression_level, worker_threads);
  if (!zstdbuf->is_open()) {
    fatal_crash("Could not open the output file '{}'", zstfilename);
  }
  return OutputFileStream(std::move(zstdbuf));
#else
  return open_uncompressed_output_file(filename);
#endif
}

// open a per-rank output file such as estimators_0000.out in the job folder. The file stays open over
// the timesteps.
[[nodiscard]] inline auto open_rank_outfile(const std::string_view basename) -> OutputFileStream {
  return open_output_file(get_jobfolder_filepath(std::format("{}_{:04d}.out", basename, globals::my_rank)));
}

#ifdef USE_ZSTD
// Compress the text into one zstd frame with a checksum. zstd decompresses a sequence of frames into the
// sequence of their texts, so a file can hold the frames of many writers one after the other.
[[nodiscard]] inline auto compress_to_zstd_frame(const std::string_view text, const int compression_level)
    -> std::string {
  const std::unique_ptr<ZSTD_CCtx, decltype(&ZSTD_freeCCtx)> cctx(ZSTD_createCCtx(), &ZSTD_freeCCtx);
  assert_always(cctx != nullptr);
  assert_always(ZSTD_isError(ZSTD_CCtx_setParameter(cctx.get(), ZSTD_c_compressionLevel, compression_level)) == 0U);
  assert_always(ZSTD_isError(ZSTD_CCtx_setParameter(cctx.get(), ZSTD_c_checksumFlag, 1)) == 0U);
  std::string frame(ZSTD_compressBound(text.size()), '\0');
  const size_t framesize = ZSTD_compress2(cctx.get(), frame.data(), frame.size(), text.data(), text.size());
  if (ZSTD_isError(framesize) != 0U) {
    fatal_crash("zstd cannot compress the output: {}", ZSTD_getErrorName(framesize));
  }
  frame.resize(framesize);
  // the bound of the compressed size can be much larger than the frame, which the caller keeps until its transfer ends
  frame.shrink_to_fit();
  return frame;
}
#endif

// An output file of the job folder for the text of the ranks, e.g. the estimators. Without the option
// WRITE_COMBINED_ALLRANK_OUT_FILES, each rank writes its own file <basename>_<rank>.out. With the option, each rank
// writes into its own text buffer. At the end of each timestep, rank 0 then writes the texts of all ranks into one file
// <basename>_allranks.out, in the order of the ranks. In a build with libzstd, each rank compresses its own text into
// one zstd frame before the send, so all ranks share the compression. Rank 0 keeps the text of only one other rank at a
// time.
class JobFolderOutputFile {
 public:
  // Every rank must call this function. rank_has_text tells if this rank writes text into the file. The header starts
  // the file of each rank, or the file of all ranks.
  void open(const std::string_view basename, const bool rank_has_text, const std::string_view header) {
    is_open = true;
    if constexpr (!WRITE_COMBINED_ALLRANK_OUT_FILES) {
      if (rank_has_text) {
        assert_always(rankfile.rdbuf() == nullptr);
        rankfile = open_rank_outfile(basename);
        rankfile << header;
      }
      return;
    }

    if (globals::my_rank != 0) {
      return;
    }
    const auto filename = get_jobfolder_filepath(std::format("{}_allranks.out", basename));
    remove_other_output_form(filename);
    const auto filepath = output_filepath(filename);
    allranksfile.open(filepath, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!allranksfile.is_open()) {
      fatal_crash("Could not open the output file '{}'", filepath);
    }
#ifdef USE_ZSTD
    // The file starts with the frame of the header. Without a header, this is an empty frame, as
    // ZstdOutputBuffer::close() writes to a file with no content. A job that stops before its first write thus leaves a
    // file that every zstd reader accepts.
    const auto headerbytes = compress_to_zstd_frame(header, ZSTD_LEVEL_DEFAULT);
#else
    const auto headerbytes = header;
#endif
    allranksfile.write(headerbytes.data(), static_cast<std::streamsize>(headerbytes.size()));
    allranksfile.flush();
    if (allranksfile.fail()) {
      fatal_crash("Could not write to the output file '{}'", filepath);
    }
  }

  [[nodiscard]] auto stream() -> std::ostream& {
    if constexpr (WRITE_COMBINED_ALLRANK_OUT_FILES) {
      return ranktext;
    }
    return rankfile;
  }

  // Put the text of the timestep into the file. Every rank must call this function, because with the option
  // WRITE_COMBINED_ALLRANK_OUT_FILES, each rank sends its text to rank 0.
  void end_timestep() {
    if constexpr (!WRITE_COMBINED_ALLRANK_OUT_FILES) {
      // one flush per timestep, so each compressed file ends a frame that a reader can decode
      rankfile.flush();
      return;
    }
    // open() runs on every rank, so all ranks have the same value
    if (is_open) {
      write_all_ranks();
    }
  }

 private:
  // Write the text of all ranks to the file and clear the buffer of each rank. Every rank must call this function.
  void write_all_ranks() {
    // std::print does not throw, e.g. when the buffer cannot grow. It sets the bad state, and the text is then not
    // complete.
    if (ranktext.fail()) {
      fatal_crash("A write to the text buffer of rank {} failed", globals::my_rank);
    }
#ifdef USE_ZSTD
    auto rankbytes =
        ranktext.view().empty() ? std::string{} : compress_to_zstd_frame(ranktext.view(), ZSTD_LEVEL_DEFAULT);
#else
    auto rankbytes = std::string(ranktext.view());
#endif
    ranktext.str({});
    ranktext.clear();

    const auto rankbytecount = static_cast<std::int64_t>(rankbytes.size());
    std::vector<std::int64_t> bytecount_of_rank(globals::my_rank == 0 ? globals::nprocs : 0);
    // with libstdc++, data() gives long* and not std::int64_t*, and the clang-tidy check of the MPI types then fails
    std::int64_t* const bytecount_of_rank_data = bytecount_of_rank.data();
    assert_always(MPI_Gather(&rankbytecount, 1, MPI_INT64_T, bytecount_of_rank_data, 1, MPI_INT64_T, 0,
                             MPI_COMM_WORLD) == MPI_SUCCESS);

    if (globals::my_rank != 0) {
      transfer_bytes(std::span<char>{rankbytes}, 0, false);
      return;
    }

    allranksfile.write(rankbytes.data(), static_cast<std::streamsize>(rankbytes.size()));
    std::string receivedbytes;
    for (int rank = 1; rank < globals::nprocs; rank++) {
      receivedbytes.resize(static_cast<size_t>(bytecount_of_rank[static_cast<size_t>(rank)]));
      transfer_bytes(std::span<char>{receivedbytes}, rank, true);
      allranksfile.write(receivedbytes.data(), static_cast<std::streamsize>(receivedbytes.size()));
    }
    allranksfile.flush();
    if (allranksfile.fail()) {
      fatal_crash("Could not write to the output file of all ranks");
    }
  }

  // Send the bytes to the other rank or receive them from it. The MPI count is a 32-bit int, so the function sends a
  // large text in chunks. MPI keeps the order of the messages between two ranks.
  static void transfer_bytes(const std::span<char> bytes, const int otherrank, const bool receive) {
    const auto nchunks = get_chunk_count(std::ssize(bytes), MPI_COUNT_MAX);
    for (auto chunk = 0Z; chunk < nchunks; chunk++) {
      const auto [chunkstart, chunksize] = get_range_chunk(std::ssize(bytes), nchunks, chunk);
      const auto chunkbytes = bytes.subspan(static_cast<size_t>(chunkstart), static_cast<size_t>(chunksize));
      const auto int_chunksize = static_cast<int>(chunkbytes.size());
      if (receive) {
        assert_always(MPI_Recv(chunkbytes.data(), int_chunksize, MPI_CHAR, otherrank, 0, MPI_COMM_WORLD,
                               MPI_STATUS_IGNORE) == MPI_SUCCESS);
      } else {
        assert_always(MPI_Send(chunkbytes.data(), int_chunksize, MPI_CHAR, otherrank, 0, MPI_COMM_WORLD) ==
                      MPI_SUCCESS);
      }
    }
  }

  bool is_open{false};
  // the file of this rank, without WRITE_COMBINED_ALLRANK_OUT_FILES
  OutputFileStream rankfile;
  // the text buffer of this rank and the file of all ranks on rank 0, with WRITE_COMBINED_ALLRANK_OUT_FILES
  std::ostringstream ranktext;
  std::ofstream allranksfile;
};

#endif  // OUTPUTFILESTREAM_H
