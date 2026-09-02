#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import artistools as at
from artistools.packets.packets import type_ids as packet_type_ids

TYPE_ESCAPE = str(packet_type_ids["TYPE_ESCAPE"])
TYPE_RPKT = str(packet_type_ids["TYPE_RPKT"])

GLOB_PATTERN = "**/packets/packets00_*.out*"


def get_type_escapetype(line: str, type_id_index: int, escape_type_id_index: int) -> tuple[str, str]:
    # the columns after escape_type_id stay in one string. A packet line has more than 30 columns,
    # and this function reads two of them.
    linesplit = line.split(maxsplit=escape_type_id_index + 1)
    return linesplit[type_id_index], linesplit[escape_type_id_index]


def is_escaped_rpkt(line: str, type_id_index: int, escape_type_id_index: int) -> bool:
    """Test if the filter keeps this packet. This function is the only definition of the filter."""
    return get_type_escapetype(line, type_id_index, escape_type_id_index) == (TYPE_ESCAPE, TYPE_RPKT)


def format_size(size_bytes: float) -> str:
    """Give a size in bytes as a short string, e.g. '1.2 GB'."""
    for unit in ("B", "kB", "MB", "GB"):
        if abs(size_bytes) < 1000:
            return f"{size_bytes:.0f} {unit}" if unit == "B" else f"{size_bytes:.1f} {unit}"
        size_bytes /= 1000
    return f"{size_bytes:.1f} TB"


def get_column_indices(firstline: str) -> tuple[int, int]:
    """Give the column index of type_id and of escape_type_id. An old file has no header."""
    if firstline.startswith("#"):
        header_cols = firstline.split()
        return header_cols.index("type_id"), header_cols.index("escape_type_id")
    return 2, 15


def find_packets_files(paths: list[Path]) -> list[Path]:
    """Find the packets files under each of the paths, most recent last."""
    found = {
        file: None
        for path in paths
        for file in path.glob(GLOB_PATTERN)
        if "parquet" not in file.name and ".tmp" not in file.name
    }
    return sorted(found, key=lambda p: p.stat().st_mtime)


def main() -> None:
    assert sys.version_info >= (3, 14), "This script requires Python 3.14 or higher."
    import compression.zstd

    parser = argparse.ArgumentParser(
        description="Filter packets files to keep only escaped rpkts.",
        epilog="The script shows what the filter would remove. Add -f to do the filtering.",
    )
    parser.add_argument("paths", nargs="*", type=Path, default=[Path()], help="folders to search (default: .)")
    parser.add_argument("-f", "--filter", action="store_true", dest="dofilter", help="do the filtering")
    parser.add_argument("--rm", action="store_true", help="remove the original files instead of a backup")
    args = parser.parse_args()

    matching_files = find_packets_files(args.paths)
    if not matching_files:
        print(f"No matching files found for pattern {GLOB_PATTERN}")
        return

    files_to_filter = 0
    total_packets = 0
    total_removed = 0
    total_size_in = 0
    total_size_out = 0.0

    for filein in matching_files:
        print(f"\n{filein}")
        linesin = at.zopen(filein).readlines()
        type_id_index, escape_type_id_index = get_column_indices(linesin[0])

        # one pass gives the statistics and the lines of the output file, in the order of the input
        outlines: list[str] = []
        packets_in = 0
        kept_packets = 0
        chars_in = 0
        chars_kept = 0
        for line in linesin:
            if line.startswith("#"):
                outlines.append(line)
                continue
            packets_in += 1
            chars_in += len(line)
            if is_escaped_rpkt(line, type_id_index, escape_type_id_index):
                outlines.append(line)
                kept_packets += 1
                chars_kept += len(line)

        if kept_packets == packets_in:
            print("  contains only escaped rpkts, skipping...")
            continue

        size_in = filein.stat().st_size
        removed_packets = packets_in - kept_packets
        removed_fraction = 1 - chars_kept / chars_in
        files_to_filter += 1
        total_packets += packets_in
        total_removed += removed_packets
        total_size_in += size_in
        print(
            f"  the filter removes {removed_packets} of {packets_in} packets "
            f"({removed_packets / packets_in * 100:.1f}%). It removes {removed_fraction * 100:.1f}% of the "
            f"uncompressed text, or about {format_size(removed_fraction * size_in)} of {format_size(size_in)}."
        )
        if not args.dofilter:
            total_size_out += (1 - removed_fraction) * size_in
            continue

        fileout_rpkt = Path(
            *filein.parts[:-1],
            filein.parts[-1].removesuffix(".zst").removesuffix(".gz").removesuffix(".xz") + ".zst",
        )
        fileout_rpkt_temp = Path(
            *fileout_rpkt.parts[:-1],
            f"{fileout_rpkt.parts[-1]}.partial.tmp.nosync",
        )
        print("  filtering rpkts...", end="", flush=True)

        with compression.zstd.open(fileout_rpkt_temp, "wt", level=12) as foutrpkt:
            foutrpkt.writelines(outlines)
        size_out = fileout_rpkt_temp.stat().st_size
        total_size_out += size_out
        print(
            f" kept {kept_packets} of {packets_in} packets ({kept_packets / packets_in * 100:.1f}%)."
            f" The size changed from {format_size(size_in)} to {format_size(size_out)}"
            f" ({size_out / size_in * 100:.1f}%)."
        )
        if args.rm:
            filein.unlink()
            print(f"  removed {filein}")
        else:
            backupfolder = filein.parent / "packets_beforefilter"
            backupfolder.mkdir(exist_ok=True)
            filein.rename(backupfolder / filein.name)
            print(f"  backed up {filein} to {backupfolder / filein.name}")

        fileout_rpkt_temp.rename(fileout_rpkt)
        print(f"  wrote {fileout_rpkt}")

    print()
    if files_to_filter == 0:
        print(f"Checked {len(matching_files)} files. Every file contains only escaped rpkts.")
        return

    removed_percent = total_removed / total_packets * 100
    if args.dofilter:
        print(
            f"Filtered {files_to_filter} of {len(matching_files)} files. The filter removed "
            f"{total_removed} of {total_packets} packets ({removed_percent:.1f}%). "
            f"The size changed from {format_size(total_size_in)} to {format_size(total_size_out)}."
        )
    else:
        print(
            f"{files_to_filter} of {len(matching_files)} files need the filter. It removes "
            f"{total_removed} of {total_packets} packets ({removed_percent:.1f}%). "
            f"The size changes from {format_size(total_size_in)} to about {format_size(total_size_out)}."
        )
        print("Use -f to do the filtering. Add --rm to remove the original files.")


if __name__ == "__main__":
    main()
