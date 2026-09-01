#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import artistools as at
from artistools.packets.packets import type_ids as packet_type_ids


def get_type_escapetype(line: str, type_id_index: int, escape_type_id_index: int) -> tuple[str, str]:
    linesplit = line.split()
    type_id = linesplit[type_id_index]
    escape_type_id = linesplit[escape_type_id_index]
    return type_id, escape_type_id


def main() -> None:
    assert sys.version_info >= (3, 14), "This script requires Python 3.14 or higher."
    import compression.zstd

    parser = argparse.ArgumentParser(description="Filter packets files to keep only escaped rpkts")
    parser.add_argument("--rm", action="store_true", help="Remove original files after processing")
    parser.add_argument("-f", action="store_true", help="Confirm performing the filtering")
    args = parser.parse_args()
    TYPE_ESCAPE = str(packet_type_ids["TYPE_ESCAPE"])
    TYPE_RPKT = str(packet_type_ids["TYPE_RPKT"])
    glob_pattern = "**/packets/packets00_*.out*"
    matching_files = sorted(Path().glob(glob_pattern), key=lambda p: p.stat().st_mtime)
    if not matching_files:
        print(f"No matching files found for pattern {glob_pattern}")
        return

    for filein in matching_files:
        if "parquet" in filein.name or ".tmp" in filein.name:
            continue
        print(f"\n{filein}")
        linesin = at.zopen(filein).readlines()
        if linesin[0].startswith("#"):
            header_cols = linesin[0].split()
            type_id_index = header_cols.index("type_id")
            escape_type_id_index = header_cols.index("escape_type_id")
        else:
            type_id_index, escape_type_id_index = 2, 15

        packets_in = 0
        kept_packets = 0
        bytes_in = 0
        bytes_kept = 0
        for line in linesin:
            if line.startswith("#"):
                continue
            packets_in += 1
            bytes_in += len(line)
            if get_type_escapetype(line, type_id_index, escape_type_id_index) == (TYPE_ESCAPE, TYPE_RPKT):
                kept_packets += 1
                bytes_kept += len(line)

        if kept_packets == packets_in:
            print("  contains only escaped rpkts, skipping...")
            continue

        removed_packets = packets_in - kept_packets
        data_reduction = (1 - bytes_kept / bytes_in) * 100
        print(
            f"  the filter removes {removed_packets} of {packets_in} packets "
            f"({removed_packets / packets_in * 100:.2f}%). Approximate data reduction: {data_reduction:.2f}%."
        )
        if not args.f:
            print("  (not filtering. Use -f to confirm filtering and --rm to remove original files)")
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
            for line in linesin:
                type_id, escape_type_id = get_type_escapetype(line, type_id_index, escape_type_id_index)
                if line.startswith("#"):
                    assert type_id == "type_id"
                    assert escape_type_id == "escape_type_id"

                    foutrpkt.write(line)
                    continue

                if (type_id, escape_type_id) == (TYPE_ESCAPE, TYPE_RPKT):
                    foutrpkt.write(line)
        size_factor = fileout_rpkt_temp.stat().st_size / filein.stat().st_size
        print(
            f" kept {kept_packets} of {packets_in} packets ({kept_packets / packets_in * 100:.2f}%)"
            f" new/old size {size_factor * 100:.2f}%"
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


if __name__ == "__main__":
    main()
