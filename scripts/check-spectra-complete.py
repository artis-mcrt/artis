#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///
"""Check the files in the spectrum completion record before a job skips exspec."""

import sys
from pathlib import Path


def spectra_complete() -> bool:
    """Require every recorded file to keep its size and modification time."""
    lines = Path(".artis-spectra-complete").read_text(encoding="utf-8").splitlines()
    if len(lines) < 4 or lines[0] != "ARTIS spectra 1" or lines[-1] != "complete":
        return False
    recorded_paths: set[Path] = set()
    for line in lines[1:-1]:
        size, timestamp, filename = line.split(maxsplit=2)
        path = Path(filename)
        stat = path.stat()
        if not path.is_file() or stat.st_size != int(size) or stat.st_mtime_ns != int(timestamp):
            return False
        recorded_paths.add(path)
    packet_paths = {path for path in Path("packets").glob("packets00_*.out*") if path.is_file()}
    return (
        Path("input.txt") in recorded_paths
        and bool(packet_paths)
        and packet_paths == {path for path in recorded_paths if path.parent == Path("packets")}
    )


def main() -> None:
    """Return failure for a missing, invalid, or obsolete completion record."""
    try:
        complete = spectra_complete()
    except (OSError, ValueError):
        complete = False
    sys.exit(0 if complete else 1)


if __name__ == "__main__":
    main()
