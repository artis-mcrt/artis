#!/usr/bin/env python3
# /// script
# requires-python = ">=3.14"
# ///
"""Combine the output files of the ranks in each job folder into one combined file for each type.

The script combines the files of the ranks of these types, e.g. nlte_0000.out, nlte_0001.out, ... into
nlte_allranks.out.zst:
- estimators;
- nlte;
- radfield;
- macroatom.

sn3d writes the same combined files with WRITE_COMBINED_ALLRANK_OUT_FILES. A combined file holds the header line of its
type and then the timesteps in their order. Each timestep holds the text of the ranks in their order. The estimator
files have no header line. Each text of one rank and one timestep is one zstd frame. The header line is the first frame,
as sn3d writes it, and without a header the first frame is empty.

As zstd does with its source files, the script keeps the files of the ranks by default. With --rm, the script combines
the files again, and then removes the files of the ranks of each type with a successful combination. --rm needs
modelgridrankassignments.out, because only that file shows that the folder holds the file of each rank. A type with an
error keeps its files of the ranks.

Do not combine the folder of a job that still runs. The combined file then holds only the timesteps up to that time.
The combined file gets the time of the newest file of a rank that the script read. A later write to a file of a rank
thus makes that file newer. The readers and this script then take the files of the ranks again.

The script combines the files of one type only when the folder holds the file of each rank that writes one.
modelgridrankassignments.out of the model folder gives these ranks:
- each rank with cells writes an estimator file;
- each rank with a cell that holds matter writes an nlte file and a radfield file;
- each rank writes a macroatom file.
sn3d writes modelgridrankassignments.out again at the start of each job. Each job of a simulation has the same number
of ranks, so the file is correct for each job folder. Without that file, the script cannot find a missing file of a
rank. An older file without the column ndo_nonempty gives no check of the nlte and radfield files.

The script also checks the order of the cells. sn3d gives each rank a contiguous range of the cells, in the order of
the ranks. The estimator files hold each model cell once in each timestep, so the cells of a timestep must be 0, 1,
2, ... in the order of the ranks. The estimator files of the ranks must also hold the same timesteps.

The nlte and radfield files hold only some cells, and a rank can write its radfield cells in the order of its OpenMP
threads. Thus each cell of these files must only come after the cells of the lower ranks. The lines of the macroatom
files have the order of the events, so the script does not check their cells.

Run the script in the model folder, e.g. "uv run artis/scripts/combine_allrank_out_files.py". It then combines the
files of each job_from_ts* folder. The arguments can also name the folders. The module compression.zstd needs Python
3.14 or a later version, and uv gets that version from the metadata above.
"""

import argparse
import contextlib
import gzip
import lzma
import os
import re
import signal
import sys
import tempfile
import typing as t
import zlib
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

# the errors of a compressed file that is damaged or ends inside a frame, e.g. because a job stopped during a write
COMPRESSION_ERRORS = (EOFError, zlib.error, zstd.ZstdError, lzma.LZMAError, gzip.BadGzipFile)

# the order of the extensions that the readers use, e.g. find_estimator_file() of artistools
RANKFILE_EXTENSIONS = ("", ".zst", ".gz", ".xz")


def get_rank_files(folder: Path, basename: str) -> dict[int, Path]:
    """Return the file of each rank in the folder for one type. A rank without text has no file."""
    names = {path.name for path in folder.glob(f"{basename}_*.out*")}
    ranks = sorted(
        {
            int(match.group(1))
            for name in names
            if (match := re.fullmatch(rf"{basename}_(\d+)\.out(?:\.zst|\.gz|\.xz)?", name))
        }
    )
    rankfiles: dict[int, Path] = {}
    for rank in ranks:
        for extension in RANKFILE_EXTENSIONS:
            if (name := f"{basename}_{rank:04d}.out{extension}") in names:
                rankfiles[rank] = folder / name
                break
    return rankfiles


@contextlib.contextmanager
def open_text(path: Path) -> Iterator[t.TextIO]:
    """Open a text file that is plain or compressed. A compressed file that is damaged or cut gives a ValueError."""
    opener = {".zst": zstd.open, ".gz": gzip.open, ".xz": lzma.open}.get(path.suffix, open)
    with opener(path, "rt", encoding="utf-8") as textfile:
        try:
            yield textfile
        except COMPRESSION_ERRORS as err:
            # e.g. a job that stopped during a write leaves a file that ends inside a compressed frame
            msg = f"{path}: the compressed file is damaged or not complete: {err}"
            raise ValueError(msg) from err


def get_rank_assignments(modelpath: Path) -> list[tuple[int, ...]] | None:
    """Return the numbers "rank nstart ndo ndo_nonempty" of each rank, or None without modelgridrankassignments.out.

    sn3d writes this file again at the start of each job. An older file has no column ndo_nonempty.
    """
    for extension in RANKFILE_EXTENSIONS:
        path = modelpath / f"modelgridrankassignments.out{extension}"
        if path.is_file():
            with open_text(path) as assignments_file:
                return [
                    tuple(int(value) for value in line.split())
                    for line in assignments_file
                    if line.strip() and not line.startswith("#")
                ]
    return None


def get_ranks_that_write(assignments: list[tuple[int, ...]], basename: str) -> set[int] | None:
    """Return the ranks that write a file of the type, or None when the assignments do not tell it.

    sn3d opens the estimator file of a rank with cells, the nlte and radfield files of a rank with a cell that holds
    matter, and the macroatom file of each rank.
    """
    if basename == "macroatom":
        return {assignment[0] for assignment in assignments}
    if basename == "estimators":
        return {assignment[0] for assignment in assignments if assignment[2] > 0}
    if all(len(assignment) >= 4 for assignment in assignments):
        return {assignment[0] for assignment in assignments if assignment[3] > 0}
    return None


def read_header(textfile: t.TextIO, path: Path, basename: str) -> str:
    """Read the header line of a file of a rank, with its newline. The estimator files have no header line."""
    if basename == "estimators":
        return ""
    header = textfile.readline()
    if not header.startswith("timestep modelgridindex ") or not header.endswith("\n"):
        msg = f"{path}: the first line is not the header line: {header!r}"
        raise ValueError(msg)
    return header


def get_timestep_and_cell(line: str, path: Path, basename: str, ncolumns: int) -> tuple[int, int] | None:
    """Return the timestep and the model cell of a line, or None for an estimator line that continues a cell.

    An estimator line that starts a cell is "timestep <n> modelgridindex <mgi> ...". Each line of the other types
    starts with "<n> <mgi>" and has the columns of the header line.
    """
    if basename == "estimators":
        if not line.startswith("timestep "):
            return None
        # the values after the words "timestep" and "modelgridindex"
        values = line.split(maxsplit=4)[1:4:2]
    else:
        values = line.split()
        # a cut line, e.g. of a job that stopped during a write, has fewer columns than the header line
        if len(values) != ncolumns or not line.endswith("\n"):
            msg = f"{path}: the line does not have the {ncolumns} columns of the header line: {line!r}"
            raise ValueError(msg)
    try:
        return int(values[0]), int(values[1])
    except (IndexError, ValueError) as err:
        msg = f"{path}: the line is not complete: {line!r}"
        raise ValueError(msg) from err


def get_timestep_texts(
    textfile: t.TextIO, path: Path, basename: str, ncolumns: int
) -> Iterator[tuple[int, str, list[int]]]:
    """Give the text of each timestep of a file of a rank, in the order of the file, and the cells of the text.

    The text of a timestep continues to the first line of the next timestep. The text of an estimator timestep thus
    holds the empty line after each cell. The cells hold the cell of each line that names one.
    """
    timestep: int | None = None
    lines: list[str] = []
    cells: list[int] = []
    for line in textfile:
        timestep_and_cell = get_timestep_and_cell(line, path, basename, ncolumns)
        if timestep_and_cell is not None:
            line_timestep, cell = timestep_and_cell
            if line_timestep != timestep:
                if timestep is not None:
                    yield timestep, "".join(lines), cells
                timestep = line_timestep
                lines = []
                cells = []
            cells.append(cell)
        elif timestep is None and line.strip():
            msg = f"{path}: the first line with content is not a timestep line: {line!r}"
            raise ValueError(msg)
        lines.append(line)
    if timestep is not None:
        yield timestep, "".join(lines), cells


def check_cell_order(
    cells: list[int], timestep: int, path: Path, basename: str, highest_cell_of_timestep: dict[int, int]
) -> None:
    """Check that the cells of one rank and one timestep come after the cells of the lower ranks, as sn3d writes them.

    highest_cell_of_timestep holds the highest cell of each timestep in the files of the lower ranks.
    """
    if basename == "macroatom":
        return
    highest_cell_of_lower_ranks = highest_cell_of_timestep.get(timestep, -1)
    previous_cell = highest_cell_of_lower_ranks
    for cell in cells:
        # each estimator cell once, in the order of the cells. A rank can write its nlte and radfield cells in any order
        cell_in_order = cell == previous_cell + 1 if basename == "estimators" else cell > highest_cell_of_lower_ranks
        if not cell_in_order:
            msg = (
                f"{path}: timestep {timestep} holds cell {cell} after cell {previous_cell}, and the lower ranks end at"
                f" cell {highest_cell_of_lower_ranks}. This is not the order of sn3d. A file of a rank is missing or"
                " stale."
            )
            raise ValueError(msg)
        previous_cell = cell
    highest_cell_of_timestep[timestep] = max(cells)


def combine_files(
    folder: Path,
    basename: str,
    rankfiles: dict[int, Path],
    assignments: list[tuple[int, ...]] | None,
    before_removal: bool,
) -> Path:
    """Write the combined file of one type in one folder, and return its path.

    The caller removes the files of the ranks after a call with before_removal. The function then combines the files
    also when the combined file is current, and it needs the rank assignments to find a missing file of a rank.
    """
    ranks_that_write = get_ranks_that_write(assignments, basename) if assignments is not None else None
    if before_removal and ranks_that_write is None:
        msg = (
            "--rm needs modelgridrankassignments.out of the model folder, because that file shows that the folder holds"
            " the file of each rank. Run the script without --rm to combine the files."
        )
        raise ValueError(msg)

    allranks_filename = f"{basename}_allranks.out"
    outpath = folder / f"{allranks_filename}.zst"
    # the readers take the plain file before the .zst file, thus a plain file is the combined file that they read
    plainpath = folder / allranks_filename
    existingpath = plainpath if plainpath.exists() else outpath if outpath.exists() else None
    if existingpath is not None and not before_removal:
        existing_mtime_ns = existingpath.stat().st_mtime_ns
        # the combined file has the time of the newest file of a rank that it holds, thus a file of a rank with the same
        # time holds no later text. artistools applies the same rule to the estimators
        if all(rankfile.stat().st_mtime_ns <= existing_mtime_ns for rankfile in rankfiles.values()):
            print(f"{folder}: {existingpath.name} exists already. The script does not change it.")
            return existingpath
        # e.g. the script ran while the job still ran, and the job then wrote more timesteps
        print(f"{folder}: a file of a rank is newer than {existingpath.name}. The script combines the files again.")

    # artistools uses a newer combined file instead of the files of the ranks. That file must thus hold the text of each
    # rank that writes one, and the estimator file must hold each model cell once in each timestep
    if ranks_that_write is not None and (missing_ranks := sorted(ranks_that_write - rankfiles.keys())):
        msg = (
            f"The {basename} files of the ranks {missing_ranks} are missing. modelgridrankassignments.out of the"
            " model folder shows that these ranks write one."
        )
        raise ValueError(msg)
    npts_model = (
        sum(assignment[2] for assignment in assignments)
        if assignments is not None and basename == "estimators"
        else None
    )

    # sn3d can add a timestep to the files of the ranks while the script runs, e.g. during the join below. The time of
    # the text before the read goes to the combined file, thus such a write makes a file of a rank newer than it
    text_mtime_ns = max(rankfile.stat().st_mtime_ns for rankfile in rankfiles.values())

    print(f"{folder}: the script combines the {basename} files of {len(rankfiles)} ranks...")

    # The script writes the frame of each rank and timestep into one temporary file, in the order of the ranks. It then
    # copies the frames of each timestep in the order of the timesteps. The memory thus holds only the text of one rank
    # and one timestep.
    with (
        tempfile.TemporaryDirectory(dir=folder, prefix=f".combine_{basename}_") as tmpdirname,
        (Path(tmpdirname) / "frames.zst").open("w+b") as frames_file,
    ):
        # the position and the size of the frame of each rank in each timestep, in the order of the ranks
        frames_of_timestep: dict[int, list[tuple[int, int]]] = {}
        timesteps_of_ranks: dict[int, set[int]] = {}
        highest_cell_of_timestep: dict[int, int] = {}
        headers: set[str] = set()
        for rank, rankfile in rankfiles.items():
            timesteps_of_rank: set[int] = set()
            with open_text(rankfile) as textfile:
                headers.add(read_header(textfile, rankfile, basename))
                if len(headers) != 1:
                    msg = f"The {basename} files of the ranks have different header lines: {sorted(headers)}"
                    raise ValueError(msg)
                ncolumns = len(next(iter(headers)).split())
                for timestep, text, cells in get_timestep_texts(textfile, rankfile, basename, ncolumns):
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
                    check_cell_order(cells, timestep, rankfile, basename, highest_cell_of_timestep)
                    frame = zstd.compress(text.encode("utf-8"), options=ZSTD_OPTIONS)
                    frames_of_timestep.setdefault(timestep, []).append((frames_file.tell(), len(frame)))
                    frames_file.write(frame)
            timesteps_of_ranks[rank] = timesteps_of_rank
            print(f"  rank {rank}: {len(timesteps_of_rank)} timesteps")

        # sn3d writes the estimators of each timestep for all ranks together. A job that stopped between the writes of
        # two ranks leaves files with different timesteps, and the combined file would then hold a timestep with some
        # ranks only. A rank of the other types writes no text in a timestep without a record for its cells.
        if basename == "estimators" and (
            incomplete_timesteps := sorted(
                timestep
                for timestep in frames_of_timestep
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
                timestep
                for timestep, highest_cell in highest_cell_of_timestep.items()
                if highest_cell + 1 != npts_model
            )
        ):
            msg = (
                f"The timesteps {short_timesteps} do not hold exactly the {npts_model} model cells, e.g. because"
                " the file of the last rank is missing."
            )
            raise ValueError(msg)

        partialpath = Path(tmpdirname) / f"{allranks_filename}.zst"
        with partialpath.open("wb") as outfile:
            # the frame of the header, or an empty frame, makes a valid zstd file also for a folder with no text
            outfile.write(zstd.compress(headers.pop().encode("utf-8"), options=ZSTD_OPTIONS))
            for timestep in sorted(frames_of_timestep):
                for position, size in frames_of_timestep[timestep]:
                    frames_file.seek(position)
                    outfile.write(frames_file.read(size))
        os.utime(partialpath, ns=(partialpath.stat().st_atime_ns, text_mtime_ns))
        # the rename puts a complete file at the name in one step
        partialpath.replace(outpath)

    # a stale plain file would hide the new .zst file from the readers
    if existingpath == plainpath:
        plainpath.unlink()
        print(f"{folder}: removed the stale {plainpath.name}.")

    print(f"{folder}: wrote {outpath.name} with {len(frames_of_timestep)} timesteps.")
    return outpath


def remove_rank_files(basename: str, rankfiles: dict[int, Path], combined_path: Path) -> None:
    """Remove the files of the ranks of one type when the combined file holds all their text."""
    # sn3d can add to a file of a rank during the combination. That file is then newer than the combined file, and it
    # holds text that the combined file lacks
    combined_mtime_ns = combined_path.stat().st_mtime_ns
    if newer_files := [path.name for path in rankfiles.values() if path.stat().st_mtime_ns > combined_mtime_ns]:
        msg = (
            f"{', '.join(newer_files)} changed after the combination, thus the script keeps the files of the ranks."
            " Run the script again after the job ends."
        )
        raise ValueError(msg)
    for path in rankfiles.values():
        path.unlink()
    print(f"{combined_path.parent}: removed the {len(rankfiles)} {basename} files of the ranks.")


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
        help="Combine the files again, and then remove the files of the ranks of each type with a successful"
        " combination. This needs modelgridrankassignments.out of the model folder.",
    )
    rmgroup.add_argument("-k", "--keep", action="store_true", help="Keep the files of the ranks (the default).")
    args = parser.parse_args()

    folders: list[Path] = args.folders or sorted(folder for folder in Path().glob("job_from_ts*") if folder.is_dir())
    if not folders:
        sys.exit("There is no job_from_ts* folder here. Give the folders as arguments.")

    # sys.exit() ends the with statements, so a stop by a signal, e.g. at the time limit of Slurm, removes the
    # temporary folder
    signal.signal(signal.SIGTERM, lambda signum, _frame: sys.exit(128 + signum))

    # an error in one type of one folder stops only that type, and the exit status then shows the error
    failures: list[str] = []
    for folder in folders:
        if not folder.is_dir():
            print(f"{folder}: the folder does not exist.", file=sys.stderr)
            failures.append(str(folder))
            continue
        rankfiles_of_type = {basename: get_rank_files(folder, basename) for basename in BASENAMES}
        if not any(rankfiles_of_type.values()):
            print(f"{folder}: the folder has no output files of ranks.")
            continue
        try:
            assignments = get_rank_assignments(folder.resolve().parent)
        except (OSError, ValueError) as err:
            print(f"{folder}: the script did not combine the files. {err}", file=sys.stderr)
            failures.append(str(folder))
            continue
        if assignments is None:
            print(
                f"{folder}: the model folder has no modelgridrankassignments.out, thus the script cannot find a missing"
                " file of a rank."
            )
        for basename, rankfiles in rankfiles_of_type.items():
            if not rankfiles:
                continue
            try:
                combined_path = combine_files(folder, basename, rankfiles, assignments, before_removal=args.rm)
                if args.rm:
                    remove_rank_files(basename, rankfiles, combined_path)
                else:
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
