#!/usr/bin/env python3
# /// script
# requires-python = ">=3.14"
# ///
"""Combine the estimator files of the ranks in each job folder into one file, estimators_allranks.out.zst.

sn3d writes the same text with WRITE_ESTIMATORS_ALLRANKS_FILE. The file holds the timesteps in their order, and in
each timestep it holds the text of the ranks in their order. Each text of one rank and one timestep is one zstd
frame, as sn3d writes it. The script keeps the files of the ranks.

Do not combine the folder of a job that still runs. The combined file then holds only the timesteps up to that time.

Run the script in the run folder, e.g. "uv run artis/scripts/combine_estimator_files.py". It then combines the files
of each job_from_ts* folder. The arguments can also name the folders. The module compression.zstd needs Python 3.14
or a later version, and uv gets that version from the metadata above.
"""

import argparse
import gzip
import lzma
import re
import tempfile
import typing as t
from collections.abc import Iterator
from pathlib import Path

from compression import zstd

# the zstd level of sn3d (ZSTD_LEVEL_DEFAULT in outputfilestream.h)
ZSTD_LEVEL = 9
ZSTD_OPTIONS = {zstd.CompressionParameter.compression_level: ZSTD_LEVEL, zstd.CompressionParameter.checksum_flag: 1}

ALLRANKS_FILENAME = "estimators_allranks.out"

# the order of the extensions that the readers use, e.g. find_estimator_file() of artistools
RANKFILE_EXTENSIONS = ("", ".zst", ".gz", ".xz")


def get_rank_files(folder: Path) -> dict[int, Path]:
    """Return the estimator file of each rank in the folder. A rank without cells has no file."""
    ranks = sorted(
        {
            int(match.group(1))
            for path in folder.glob("estimators_*.out*")
            if (match := re.fullmatch(r"estimators_(\d+)\.out(?:\.zst|\.gz|\.xz)?", path.name))
        }
    )
    rankfiles: dict[int, Path] = {}
    for rank in ranks:
        for extension in RANKFILE_EXTENSIONS:
            path = folder / f"estimators_{rank:04d}.out{extension}"
            if path.is_file():
                rankfiles[rank] = path
                break
    return rankfiles


def open_text(path: Path) -> t.TextIO:
    """Open a text file that is plain or compressed."""
    if path.suffix == ".zst":
        return zstd.open(path, "rt", encoding="utf-8")
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    if path.suffix == ".xz":
        return lzma.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def get_timestep_texts(path: Path) -> Iterator[tuple[int, str]]:
    """Give the text of each timestep of a rank file, in the order of the file.

    The line "timestep <n> modelgridindex <mgi> ..." starts the text of each cell. The text of a timestep continues
    to the first line of the next timestep, so it holds the empty line after each cell.
    """
    timestep: int | None = None
    lines: list[str] = []
    with open_text(path) as rankfile:
        for line in rankfile:
            if line.startswith("timestep "):
                linetimestep = int(line.split()[1])
                if linetimestep != timestep:
                    if timestep is not None:
                        yield timestep, "".join(lines)
                    timestep = linetimestep
                    lines = []
            elif timestep is None and line.strip():
                msg = f"{path}: the first line with content is not a timestep line: {line!r}"
                raise ValueError(msg)
            lines.append(line)
    if timestep is not None:
        yield timestep, "".join(lines)


def combine_folder(folder: Path) -> None:
    """Write the combined estimator file of one folder."""
    outpath = folder / f"{ALLRANKS_FILENAME}.zst"
    if outpath.exists() or (folder / ALLRANKS_FILENAME).exists():
        print(f"{folder}: {ALLRANKS_FILENAME} exists already. The script does not change it.")
        return

    rankfiles = get_rank_files(folder)
    if not rankfiles:
        print(f"{folder}: the folder has no estimator files of ranks.")
        return

    print(f"{folder}: the script combines the estimator files of {len(rankfiles)} ranks...")

    # Each timestep gets a temporary file, and the script adds the frame of each rank in the order of the ranks.
    # The script then joins the temporary files in the order of the timesteps. The memory thus holds only
    # the text of one rank and one timestep.
    with tempfile.TemporaryDirectory(dir=folder, prefix=".combine_estimators_") as tmpdirname:
        tmpdir = Path(tmpdirname)
        timestep_files: dict[int, t.BinaryIO] = {}
        try:
            for rank, rankfile in rankfiles.items():
                timesteps_of_rank: set[int] = set()
                for timestep, text in get_timestep_texts(rankfile):
                    if timestep in timesteps_of_rank:
                        msg = f"{rankfile}: timestep {timestep} occurs in two separate parts of the file"
                        raise ValueError(msg)
                    timesteps_of_rank.add(timestep)
                    if timestep not in timestep_files:
                        timestep_files[timestep] = (tmpdir / f"timestep_{timestep:05d}.zst").open("wb")
                    timestep_files[timestep].write(zstd.compress(text.encode("utf-8"), options=ZSTD_OPTIONS))
                print(f"  rank {rank}: {len(timesteps_of_rank)} timesteps")
        finally:
            for timestep_file in timestep_files.values():
                timestep_file.close()

        partialpath = tmpdir / f"{ALLRANKS_FILENAME}.zst"
        with partialpath.open("wb") as outfile:
            for timestep in sorted(timestep_files):
                with (tmpdir / f"timestep_{timestep:05d}.zst").open("rb") as timestep_file:
                    while chunk := timestep_file.read(1 << 24):
                        outfile.write(chunk)
        # the rename puts a complete file at the name in one step
        partialpath.replace(outpath)

    print(f"{folder}: wrote {outpath.name} with {len(timestep_files)} timesteps.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "folders", nargs="*", type=Path, help="The folders with the estimator files. The default is job_from_ts*."
    )
    args = parser.parse_args()

    folders: list[Path] = args.folders or sorted(folder for folder in Path().glob("job_from_ts*") if folder.is_dir())
    if not folders:
        print("There is no job_from_ts* folder here. Give the folders as arguments.")
    for folder in folders:
        combine_folder(folder)


if __name__ == "__main__":
    main()
