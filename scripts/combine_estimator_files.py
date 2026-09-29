#!/usr/bin/env python3
# /// script
# requires-python = ">=3.14"
# ///
"""Combine the estimator files of the ranks in each job folder into one file, estimators_allranks.out.zst.

sn3d writes the same text with WRITE_ESTIMATORS_ALLRANKS_FILE. The file holds the timesteps in their order, and in
each timestep it holds the text of the ranks in their order. Each text of one rank and one timestep is one zstd
frame, and the file starts with one empty frame, as sn3d writes it. The script keeps the files of the ranks.

Do not combine the folder of a job that still runs. The combined file then holds only the timesteps up to that time.
When a file of a rank is newer than the combined file, the script combines the files again.

The script combines a folder only when it holds the file of each rank with cells. modelgridrankassignments.out of
the run folder gives these ranks. Without that file, the ranks of the files must start at 0 and have no gap. The files
of all the ranks must also hold the same timesteps.

Run the script in the run folder, e.g. "uv run artis/scripts/combine_estimator_files.py". It then combines the files
of each job_from_ts* folder. The arguments can also name the folders. The module compression.zstd needs Python 3.14
or a later version, and uv gets that version from the metadata above.
"""

import argparse
import gzip
import lzma
import re
import shutil
import sys
import tempfile
import typing as t
from collections.abc import Iterator
from pathlib import Path

try:
    from compression import zstd
except ImportError:
    sys.exit(
        "combine_estimator_files.py needs Python 3.14 or a later version for the module compression.zstd."
        ' Run it with uv, e.g. "uv run artis/scripts/combine_estimator_files.py".'
    )

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


def get_ranks_with_cells(runfolder: Path) -> set[int] | None:
    """Return the ranks that write an estimator file, or None when the run folder has no modelgridrankassignments.out.

    sn3d writes the line "rank nstart ndo ndo_nonempty" for each rank. A rank with ndo > 0 writes the estimators.
    """
    for extension in RANKFILE_EXTENSIONS:
        path = runfolder / f"modelgridrankassignments.out{extension}"
        if path.is_file():
            with open_text(path) as assignmentfile:
                return {
                    int(columns[0])
                    for columns in (line.split() for line in assignmentfile if not line.startswith("#"))
                    if len(columns) >= 3 and int(columns[2]) > 0
                }
    return None


def get_timestep_texts(path: Path) -> Iterator[tuple[int, str]]:
    """Give the text of each timestep of a rank file, in the order of the file.

    The line "timestep <n> modelgridindex <mgi> ..." starts the text of each cell. The text of a timestep continues
    to the first line of the next timestep, so it holds the empty line after each cell.
    """
    timestep: int | None = None
    lines: list[str] = []
    try:
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
    except (EOFError, zstd.ZstdError, lzma.LZMAError, gzip.BadGzipFile) as err:
        # e.g. a job that stopped during a write leaves a file that ends inside a compressed frame
        msg = f"{path}: the compressed file is not complete: {err}"
        raise ValueError(msg) from err
    if timestep is not None:
        yield timestep, "".join(lines)


def combine_folder(folder: Path) -> None:
    """Write the combined estimator file of one folder."""
    rankfiles = get_rank_files(folder)
    outpath = folder / f"{ALLRANKS_FILENAME}.zst"
    # the readers take the plain file before the .zst file, thus a plain file is the combined file that they read
    plainpath = folder / ALLRANKS_FILENAME
    existingpath = plainpath if plainpath.exists() else outpath if outpath.exists() else None
    if existingpath is not None:
        if all(rankfile.stat().st_mtime < existingpath.stat().st_mtime for rankfile in rankfiles.values()):
            print(f"{folder}: {existingpath.name} exists already. The script does not change it.")
            return
        # e.g. the script ran while the job still ran, and the job then wrote more timesteps
        print(f"{folder}: a file of a rank is newer than {existingpath.name}. The script combines the files again.")

    if not rankfiles:
        print(f"{folder}: the folder has no estimator files of ranks.")
        return

    # artistools uses a newer combined file instead of the files of the ranks.
    # The combined file must thus hold the cells of each rank.
    ranks_with_cells = get_ranks_with_cells(folder.resolve().parent)
    if ranks_with_cells is None:
        ranks_with_cells = set(range(max(rankfiles) + 1))
    if missing_ranks := sorted(ranks_with_cells - rankfiles.keys()):
        msg = f"The folder has no estimator file of the ranks {missing_ranks}."
        raise ValueError(msg)

    print(f"{folder}: the script combines the estimator files of {len(rankfiles)} ranks...")

    # Each timestep gets a temporary file, and the script adds the frame of each rank in the order of the ranks.
    # The script then joins the temporary files in the order of the timesteps. The memory thus holds only
    # the text of one rank and one timestep. Each write opens and closes its file, so few files are open at a time.
    with tempfile.TemporaryDirectory(dir=folder, prefix=".combine_estimators_") as tmpdirname:
        tmpdir = Path(tmpdirname)
        timesteps: set[int] = set()
        timesteps_of_ranks: dict[int, set[int]] = {}
        for rank, rankfile in rankfiles.items():
            timesteps_of_rank: set[int] = set()
            for timestep, text in get_timestep_texts(rankfile):
                if timestep in timesteps_of_rank:
                    msg = f"{rankfile}: timestep {timestep} occurs in two separate parts of the file"
                    raise ValueError(msg)
                timesteps_of_rank.add(timestep)
                with (tmpdir / f"timestep_{timestep:05d}.zst").open("ab") as timestep_file:
                    timestep_file.write(zstd.compress(text.encode("utf-8"), options=ZSTD_OPTIONS))
            timesteps |= timesteps_of_rank
            timesteps_of_ranks[rank] = timesteps_of_rank
            print(f"  rank {rank}: {len(timesteps_of_rank)} timesteps")

        # sn3d writes each timestep for all ranks together. A job that stopped between the writes of two ranks
        # leaves files with different timesteps, and the combined file would then hold a timestep with some ranks only
        if incomplete_timesteps := sorted(
            timestep
            for timestep in timesteps
            if any(timestep not in rank_timesteps for rank_timesteps in timesteps_of_ranks.values())
        ):
            ranks_without = sorted(
                rank
                for rank, rank_timesteps in timesteps_of_ranks.items()
                if any(ts not in rank_timesteps for ts in incomplete_timesteps)
            )
            msg = f"The files of the ranks {ranks_without} do not hold the timesteps {incomplete_timesteps}."
            raise ValueError(msg)

        partialpath = tmpdir / f"{ALLRANKS_FILENAME}.zst"
        with partialpath.open("wb") as outfile:
            # the empty frame makes a valid zstd file also for a folder with no text
            outfile.write(zstd.compress(b"", options=ZSTD_OPTIONS))
            for timestep in sorted(timesteps):
                with (tmpdir / f"timestep_{timestep:05d}.zst").open("rb") as timestep_file:
                    shutil.copyfileobj(timestep_file, outfile, 1 << 24)
        # the rename puts a complete file at the name in one step
        partialpath.replace(outpath)

    # a stale plain file would hide the new .zst file from the readers
    if existingpath == plainpath:
        plainpath.unlink()
        print(f"{folder}: removed the stale {plainpath.name}.")

    print(f"{folder}: wrote {outpath.name} with {len(timesteps)} timesteps.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "folders", nargs="*", type=Path, help="The folders with the estimator files. The default is job_from_ts*."
    )
    args = parser.parse_args()

    folders: list[Path] = args.folders or sorted(folder for folder in Path().glob("job_from_ts*") if folder.is_dir())
    if not folders:
        print("There is no job_from_ts* folder here. Give the folders as arguments.")

    # an error in one folder stops only that folder, and the exit status then shows the error
    failedfolders: list[Path] = []
    for folder in folders:
        try:
            combine_folder(folder)
        except (OSError, ValueError) as err:
            print(f"{folder}: the script did not combine the folder. {err}", file=sys.stderr)
            failedfolders.append(folder)
    if failedfolders:
        sys.exit(f"The script did not combine these folders: {', '.join(map(str, failedfolders))}")


if __name__ == "__main__":
    main()
