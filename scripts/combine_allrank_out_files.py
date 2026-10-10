#!/usr/bin/env python3
# /// script
# requires-python = ">=3.14"
# ///
"""Combine the output files of the ranks in each job folder into one file of all ranks for each type of file.

The types are estimators, nlte, radfield, and macroatom. For example, the script combines nlte_0000.out,
nlte_0001.out, ... into nlte_allranks.out.zst. sn3d writes the same files with WRITE_COMBINED_ALLRANK_OUT_FILES.
A file of all ranks holds the header line of its type, then the timesteps in their order, and in each timestep the
text of the ranks in their order. The estimator files have no header line. Each text of one rank and one timestep
is one zstd frame, and the header line is the first frame, as sn3d writes it. Without a header, the first frame is
empty.

As zstd does with its source files, the script keeps the files of the ranks by default. With --rm, it removes the
files of the ranks of one type after a successful combination, i.e. when the combined file of that type is complete
and current. A type with an error keeps its files of the ranks.

Do not combine the folder of a job that still runs. The combined file then holds only the timesteps up to that time.
The combined file gets the time of the newest file of a rank that the script read. A later write to a file of a rank
thus makes that file newer, and the readers and this script then take the files of the ranks again.

The script combines the files of one type only when they hold the cells in the order of sn3d. sn3d gives each rank
a contiguous range of the cells, in the order of the ranks. The estimator files hold each model cell once in each
timestep, so the cells of a timestep must be 0, 1, 2, ... in the order of the ranks. modelgridrankassignments.out of
the model folder gives the number of model cells. Without that file, the script cannot find a missing estimator
file of the last rank. The estimator files of all the ranks must also hold the same timesteps. The nlte and
radfield files hold only some cells, so their cells must only increase over the ranks. The lines of the macroatom
files have the order of the events, so the script does not check their cells.

Run the script in the model folder, e.g. "uv run artis/scripts/combine_allrank_out_files.py". It then combines the
files of each job_from_ts* folder. The arguments can also name the folders. The module compression.zstd needs Python
3.14 or a later version, and uv gets that version from the metadata above.
"""

import argparse
import gzip
import lzma
import os
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
        "combine_allrank_out_files.py needs Python 3.14 or a later version for the module compression.zstd."
        ' Run it with uv, e.g. "uv run artis/scripts/combine_allrank_out_files.py".'
    )

# the zstd level of sn3d (ZSTD_LEVEL_DEFAULT in outputfilestream.h)
ZSTD_LEVEL = 9
ZSTD_OPTIONS = {zstd.CompressionParameter.compression_level: ZSTD_LEVEL, zstd.CompressionParameter.checksum_flag: 1}

# the types of the output files of the ranks (JOBFOLDER_OUTFILE_BASENAMES in mpi_logging.h)
BASENAMES = ("estimators", "nlte", "radfield", "macroatom")

# the errors of a compressed file that ends inside a frame, e.g. because a job stopped during a write
COMPRESSION_ERRORS = (EOFError, zstd.ZstdError, lzma.LZMAError, gzip.BadGzipFile)

# the order of the extensions that the readers use, e.g. find_estimator_file() of artistools
RANKFILE_EXTENSIONS = ("", ".zst", ".gz", ".xz")


def get_rank_files(folder: Path, basename: str) -> dict[int, Path]:
    """Return the file of each rank in the folder for one type. A rank without text has no file."""
    ranks = sorted(
        {
            int(match.group(1))
            for path in folder.glob(f"{basename}_*.out*")
            if (match := re.fullmatch(rf"{basename}_(\d+)\.out(?:\.zst|\.gz|\.xz)?", path.name))
        }
    )
    rankfiles: dict[int, Path] = {}
    for rank in ranks:
        for extension in RANKFILE_EXTENSIONS:
            path = folder / f"{basename}_{rank:04d}.out{extension}"
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


def get_npts_model(modelpath: Path) -> int | None:
    """Return the number of model cells, or None when the model folder has no modelgridrankassignments.out.

    sn3d writes this file again at the start of each job, with the line "rank nstart ndo ndo_nonempty" for each rank.
    A job with a different number of ranks gives other ranks, but the sum of ndo is always the number of model cells.
    """
    for extension in RANKFILE_EXTENSIONS:
        path = modelpath / f"modelgridrankassignments.out{extension}"
        if path.is_file():
            try:
                with open_text(path) as assignmentfile:
                    return sum(
                        int(columns[2])
                        for columns in (line.split() for line in assignmentfile if not line.startswith("#"))
                        if len(columns) >= 3
                    )
            except COMPRESSION_ERRORS as err:
                msg = f"{path}: the compressed file is not complete: {err}"
                raise ValueError(msg) from err
    return None


def get_header(path: Path, basename: str) -> str:
    """Return the header line of a file of a rank, with its newline. The estimator files have no header line."""
    if basename == "estimators":
        return ""
    try:
        with open_text(path) as rankfile:
            header = rankfile.readline()
    except COMPRESSION_ERRORS as err:
        msg = f"{path}: the compressed file is not complete: {err}"
        raise ValueError(msg) from err
    if not header.startswith("timestep modelgridindex ") or not header.endswith("\n"):
        msg = f"{path}: the first line is not the header line: {header!r}"
        raise ValueError(msg)
    return header


def get_timestep_and_cell(line: str, path: Path, basename: str) -> tuple[int, int]:
    """Return the timestep and the model cell of a line.

    An estimator line that starts a cell is "timestep <n> modelgridindex <mgi> ...". Each line of the other types
    starts with "<n> <mgi>".
    """
    columns = line.split()
    try:
        if basename == "estimators":
            return int(columns[1]), int(columns[3])
        return int(columns[0]), int(columns[1])
    except (IndexError, ValueError) as err:
        # e.g. a job that stopped during a write leaves a cut line at the end of the file
        msg = f"{path}: the line is not complete: {line!r}"
        raise ValueError(msg) from err


def get_line_timestep(line: str, path: Path, basename: str, ncolumns: int) -> int | None:
    """Return the timestep of a line, or None for an estimator line that continues a cell."""
    if basename == "estimators":
        return get_timestep_and_cell(line, path, basename)[0] if line.startswith("timestep ") else None
    # each line has the columns of the header line, thus a line with other columns is a cut line
    if not line.endswith("\n") or len(line.split()) != ncolumns:
        msg = f"{path}: the line does not have the {ncolumns} columns of the header line: {line!r}"
        raise ValueError(msg)
    return get_timestep_and_cell(line, path, basename)[0]


def get_timestep_texts(path: Path, basename: str, header: str) -> Iterator[tuple[int, str]]:
    """Give the text of each timestep of a file of a rank, in the order of the file.

    The text of a timestep continues to the first line of the next timestep. The text of an estimator timestep thus
    holds the empty line after each cell.
    """
    ncolumns = len(header.split())
    timestep: int | None = None
    lines: list[str] = []
    try:
        with open_text(path) as rankfile:
            if header:
                rankfile.readline()
            for line in rankfile:
                linetimestep = get_line_timestep(line, path, basename, ncolumns)
                if linetimestep is not None and linetimestep != timestep:
                    if timestep is not None:
                        yield timestep, "".join(lines)
                    timestep = linetimestep
                    lines = []
                elif timestep is None and line.strip():
                    msg = f"{path}: the first line with content is not a timestep line: {line!r}"
                    raise ValueError(msg)
                lines.append(line)
    except COMPRESSION_ERRORS as err:
        # e.g. a job that stopped during a write leaves a file that ends inside a compressed frame
        msg = f"{path}: the compressed file is not complete: {err}"
        raise ValueError(msg) from err
    if timestep is not None:
        yield timestep, "".join(lines)


def check_cell_order(
    text: str, timestep: int, path: Path, basename: str, last_cell_of_timestep: dict[int, int]
) -> None:
    """Check that the text of one rank and one timestep continues the cells of the lower ranks in the order of sn3d.

    last_cell_of_timestep holds the last cell of each timestep in the files of the lower ranks.
    """
    if basename == "macroatom":
        return
    last_cell_of_lower_ranks = last_cell_of_timestep.get(timestep, -1)
    previous_cell = last_cell_of_lower_ranks
    for line in text.splitlines():
        if basename == "estimators" and not line.startswith("timestep "):
            continue
        _, cell = get_timestep_and_cell(line, path, basename)
        if basename == "estimators":
            # each cell once, in the order of the cells
            cell_in_order = cell == previous_cell + 1
        else:
            # the lines of a cell together, in the order of the cells
            cell_in_order = cell > last_cell_of_lower_ranks and cell >= previous_cell
        if not cell_in_order:
            msg = (
                f"{path}: timestep {timestep} holds cell {cell} after cell {previous_cell}, which is not the order of"
                " sn3d. A file of a rank is missing or stale."
            )
            raise ValueError(msg)
        previous_cell = cell
    last_cell_of_timestep[timestep] = previous_cell


def combine_files(folder: Path, basename: str) -> list[Path]:
    """Write the combined file of one type in one folder, and return the files of the ranks that it holds."""
    rankfiles = get_rank_files(folder, basename)
    allranks_filename = f"{basename}_allranks.out"
    outpath = folder / f"{allranks_filename}.zst"
    # the readers take the plain file before the .zst file, thus a plain file is the combined file that they read
    plainpath = folder / allranks_filename
    existingpath = plainpath if plainpath.exists() else outpath if outpath.exists() else None
    if existingpath is not None:
        # the combined file has the time of the newest file of a rank that it holds, thus a file of a rank with the same
        # time holds no later text. artistools applies the same rule to the estimators
        if all(rankfile.stat().st_mtime_ns <= existingpath.stat().st_mtime_ns for rankfile in rankfiles.values()):
            print(f"{folder}: {existingpath.name} exists already. The script does not change it.")
            return list(rankfiles.values())
        # e.g. the script ran while the job still ran, and the job then wrote more timesteps
        print(f"{folder}: a file of a rank is newer than {existingpath.name}. The script combines the files again.")

    if not rankfiles:
        return []

    headers = {get_header(rankfile, basename) for rankfile in rankfiles.values()}
    if len(headers) != 1:
        msg = f"The {basename} files of the ranks have different header lines: {sorted(headers)}"
        raise ValueError(msg)
    header = headers.pop()

    # artistools uses a newer combined estimator file instead of the files of the ranks. That file must thus hold
    # each model cell once in each timestep, and no cell of a stale file of a rank that the job did not have
    npts_model = get_npts_model(folder.resolve().parent) if basename == "estimators" else None

    # sn3d can add a timestep to the files of the ranks while the script runs, e.g. during the join below. The time of
    # the text before the read goes to the combined file, thus such a write makes a file of a rank newer than it
    text_mtime_ns = max(rankfile.stat().st_mtime_ns for rankfile in rankfiles.values())

    print(f"{folder}: the script combines the {basename} files of {len(rankfiles)} ranks...")

    # Each timestep gets a temporary file, and the script adds the frame of each rank in the order of the ranks.
    # The script then joins the temporary files in the order of the timesteps. The memory thus holds only
    # the text of one rank and one timestep. Each write opens and closes its file, so few files are open at a time.
    with tempfile.TemporaryDirectory(dir=folder, prefix=f".combine_{basename}_") as tmpdirname:
        tmpdir = Path(tmpdirname)
        timesteps: set[int] = set()
        timesteps_of_ranks: dict[int, set[int]] = {}
        last_cell_of_timestep: dict[int, int] = {}
        for rank, rankfile in rankfiles.items():
            timesteps_of_rank: set[int] = set()
            for timestep, text in get_timestep_texts(rankfile, basename, header):
                if timestep in timesteps_of_rank:
                    msg = f"{rankfile}: timestep {timestep} occurs in two separate parts of the file"
                    raise ValueError(msg)
                timesteps_of_rank.add(timestep)
                # sn3d writes an empty line after each cell, thus a text without it ends inside a cell
                if basename == "estimators" and not text.endswith("\n\n"):
                    msg = (
                        f"{rankfile}: timestep {timestep} does not end with the empty line after a cell."
                        " The file is not complete, e.g. because the job stopped during a write."
                    )
                    raise ValueError(msg)
                check_cell_order(text, timestep, rankfile, basename, last_cell_of_timestep)
                with (tmpdir / f"timestep_{timestep:05d}.zst").open("ab") as timestep_file:
                    timestep_file.write(zstd.compress(text.encode("utf-8"), options=ZSTD_OPTIONS))
            timesteps |= timesteps_of_rank
            timesteps_of_ranks[rank] = timesteps_of_rank
            print(f"  rank {rank}: {len(timesteps_of_rank)} timesteps")

        # sn3d writes the estimators of each timestep for all ranks together. A job that stopped between the writes of
        # two ranks leaves files with different timesteps, and the combined file would then hold a timestep with some
        # ranks only. A rank of the other types writes no text in a timestep without a record for its cells.
        if basename == "estimators" and (
            incomplete_timesteps := sorted(
                timestep
                for timestep in timesteps
                if any(timestep not in rank_timesteps for rank_timesteps in timesteps_of_ranks.values())
            )
        ):
            ranks_without = sorted(
                rank
                for rank, rank_timesteps in timesteps_of_ranks.items()
                if any(ts not in rank_timesteps for ts in incomplete_timesteps)
            )
            msg = f"The files of the ranks {ranks_without} do not hold the timesteps {incomplete_timesteps}."
            raise ValueError(msg)

        if npts_model is not None and (
            short_timesteps := sorted(
                timestep for timestep, last_cell in last_cell_of_timestep.items() if last_cell + 1 != npts_model
            )
        ):
            msg = (
                f"The timesteps {short_timesteps} do not hold exactly the {npts_model} model cells, e.g. because"
                " the file of the last rank is missing."
            )
            raise ValueError(msg)

        partialpath = tmpdir / f"{allranks_filename}.zst"
        with partialpath.open("wb") as outfile:
            # the frame of the header, or an empty frame, makes a valid zstd file also for a folder with no text
            outfile.write(zstd.compress(header.encode("utf-8"), options=ZSTD_OPTIONS))
            for timestep in sorted(timesteps):
                with (tmpdir / f"timestep_{timestep:05d}.zst").open("rb") as timestep_file:
                    shutil.copyfileobj(timestep_file, outfile, 1 << 24)
        os.utime(partialpath, ns=(partialpath.stat().st_atime_ns, text_mtime_ns))
        # the rename puts a complete file at the name in one step
        partialpath.replace(outpath)

    # a stale plain file would hide the new .zst file from the readers
    if existingpath == plainpath:
        plainpath.unlink()
        print(f"{folder}: removed the stale {plainpath.name}.")

    print(f"{folder}: wrote {outpath.name} with {len(timesteps)} timesteps.")
    return list(rankfiles.values())


def remove_rank_files(folder: Path, basename: str, rankfiles: list[Path]) -> None:
    """Remove the files of the ranks of one type when the combined file holds all their text."""
    allranks_filename = f"{basename}_allranks.out"
    combinedpath = next(
        path for path in (folder / allranks_filename, folder / f"{allranks_filename}.zst") if path.exists()
    )
    # sn3d can add to a file of a rank during the combination. That file is then newer than the combined file, and it
    # holds text that the combined file lacks
    combined_mtime_ns = combinedpath.stat().st_mtime_ns
    if newerfiles := [rankfile.name for rankfile in rankfiles if rankfile.stat().st_mtime_ns > combined_mtime_ns]:
        msg = (
            f"{', '.join(newerfiles)} changed after the combination, thus the script keeps the files of the ranks."
            " Run the script again after the job ends."
        )
        raise ValueError(msg)
    for rankfile in rankfiles:
        rankfile.unlink()
    print(f"{folder}: removed the {len(rankfiles)} {basename} files of the ranks.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "folders", nargs="*", type=Path, help="The folders with the files of the ranks. The default is job_from_ts*."
    )
    # the options work as the same options of zstd: the default keeps the source files
    rmgroup = parser.add_mutually_exclusive_group()
    rmgroup.add_argument(
        "--rm",
        action="store_true",
        help="Remove the files of the ranks of one type after a successful combination, i.e. when the combined file of"
        " that type is complete and current.",
    )
    rmgroup.add_argument("-k", "--keep", action="store_true", help="Keep the files of the ranks (the default).")
    args = parser.parse_args()

    folders: list[Path] = args.folders or sorted(folder for folder in Path().glob("job_from_ts*") if folder.is_dir())
    if not folders:
        sys.exit("There is no job_from_ts* folder here. Give the folders as arguments.")

    # an error in one type of one folder stops only that type, and the exit status then shows the error
    failures: list[str] = []
    for folder in folders:
        if not folder.is_dir():
            print(f"{folder}: the folder does not exist.", file=sys.stderr)
            failures.append(str(folder))
            continue
        if not any(get_rank_files(folder, basename) for basename in BASENAMES):
            print(f"{folder}: the folder has no output files of ranks.")
            continue
        for basename in BASENAMES:
            try:
                rankfiles = combine_files(folder, basename)
                if rankfiles and args.rm:
                    remove_rank_files(folder, basename, rankfiles)
                elif rankfiles:
                    print(
                        f"{folder}: kept the {len(rankfiles)} {basename} files of the ranks. Give --rm to remove them."
                    )
            except (OSError, ValueError) as err:
                print(f"{folder}: the script did not combine the {basename} files. {err}", file=sys.stderr)
                failures.append(f"{folder} ({basename})")
    if failures:
        sys.exit(f"The script did not combine these files: {', '.join(failures)}")


if __name__ == "__main__":
    main()
