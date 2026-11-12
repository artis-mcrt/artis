#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///
"""Check the files in the spectrum completion record before a job skips exspec."""

import argparse
import itertools
import shutil
import subprocess  # ruff: ignore[suspicious-subprocess-import]  # zstd supplies the file compression.
import sys
import typing as t
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

RECORD_PATH = Path(".artis-spectra-complete")


class RecordedFile(t.NamedTuple):
    """The path, size, and modification time of one completed file."""

    path: Path
    size: int
    timestamp: int


def read_record() -> list[RecordedFile]:
    """Read a complete record or stop with an error."""
    lines = RECORD_PATH.read_text(encoding="utf-8").splitlines()
    if len(lines) < 4 or lines[0] != "ARTIS spectra 1" or lines[-1] != "complete":
        msg = "The spectrum completion record is incomplete"
        raise ValueError(msg)
    files: list[RecordedFile] = []
    for line in lines[1:-1]:
        size, timestamp, filename = line.split(maxsplit=2)
        files.append(RecordedFile(Path(filename), int(size), int(timestamp)))
    return files


def file_matches_record(file: RecordedFile) -> bool:
    """Check the recorded size and modification time of one file."""
    stat = file.path.stat()
    return file.path.is_file() and stat.st_size == file.size and stat.st_mtime_ns == file.timestamp


def files_match_record(files: list[RecordedFile]) -> bool:
    """Require every file and the complete set of packet files to match the record."""
    if not all(file_matches_record(file) for file in files):
        return False
    recorded_paths = {file.path for file in files}
    packet_paths = {path for path in Path("packets").glob("packets00_*.out*") if path.is_file()}
    return (
        bool({Path("input.txt"), Path("input.txt.zst")} & recorded_paths)
        and bool(packet_paths)
        and packet_paths == {path for path in recorded_paths if path.parent == Path("packets")}
    )


def compress_recorded_file(file: RecordedFile, zstd: str) -> RecordedFile:
    """Keep a plain source until zstd succeeds and its original metadata still matches."""
    if file.size <= 200 * 1024 or file.path.suffix not in {".txt", ".out"}:
        return file
    if not file_matches_record(file):
        msg = f"{file.path} changed before compression"
        raise ValueError(msg)
    subprocess.run(  # ruff: ignore[subprocess-without-shell-equals-true]  # The file path is one argument.
        [zstd, "-T1", "-13", "-v", "-f", "--", str(file.path)], check=True
    )
    if not file_matches_record(file):
        msg = f"{file.path} changed during compression"
        raise ValueError(msg)
    compressed_path = file.path.with_name(file.path.name + ".zst")
    stat = compressed_path.stat()
    file.path.unlink()
    return RecordedFile(compressed_path, stat.st_size, stat.st_mtime_ns)


def compress_files(files: list[RecordedFile], workers: int) -> None:
    """Update a valid completion record only after its expected file conversions succeed."""
    zstd = shutil.which("zstd")
    if zstd is None:
        msg = "Could not find zstd"
        raise FileNotFoundError(msg)
    RECORD_PATH.unlink()
    with ThreadPoolExecutor(max_workers=workers) as executor:
        compressed_files = list(executor.map(compress_recorded_file, files, itertools.repeat(zstd)))
    if not files_match_record(compressed_files):
        msg = "A recorded file changed during compression"
        raise ValueError(msg)
    lines = ["ARTIS spectra 1"]
    lines.extend(f"{file.size} {file.timestamp} {file.path}" for file in compressed_files)
    lines.append("complete")
    partial_path = RECORD_PATH.with_name(RECORD_PATH.name + ".tmp")
    partial_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    partial_path.replace(RECORD_PATH)


def check_or_compress(workers: int | None) -> bool:
    """Check the record, then apply only the requested file conversions."""
    # Older runs have no record. The normal post-process commands can still compress their files.
    if workers is not None and not RECORD_PATH.exists():
        return True
    files = read_record()
    if not files_match_record(files):
        return False
    if workers is not None:
        compress_files(files, workers)
    return True


def main() -> None:
    """Return failure for a missing, invalid, or obsolete completion record."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compress", type=int, metavar="WORKERS", help="Compress the recorded plain files")
    args = parser.parse_args()
    try:
        complete = check_or_compress(args.compress)
    except (OSError, ValueError, subprocess.CalledProcessError):
        complete = False
    sys.exit(0 if complete else 1)


if __name__ == "__main__":
    main()
