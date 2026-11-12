#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///
"""Merge the files of the direction bins into light_curve_res.out, spec_res.out, and specpol_res.out."""

import re
import shutil
import sys
import tempfile
from pathlib import Path


def resfilepath(basefolder: Path | str, fileprefix: str, abin: int) -> Path:
    """Return the path of the file of one direction bin. The .zst file is the fallback of the plain file."""
    plainpath = Path(basefolder, f"{fileprefix}_res_{abin:02d}.out")
    zstpath = plainpath.with_name(plainpath.name + ".zst")
    return zstpath if not plainpath.is_file() and zstpath.is_file() else plainpath


def get_mabins() -> int:
    """Give the number of direction bins."""
    # one output file per direction bin (MABINS = NPHIBINS * NCOSTHETABINS), with the expected count
    # parsed from the spectrum_lightcurve.h of the source folder that this script belongs to
    spectrum_lightcurve_h = (Path(__file__).resolve().parent.parent / "spectrum_lightcurve.h").read_text(
        encoding="utf-8"
    )
    bincounts = {
        name: int(value)
        for name, value in re.findall(r"constexpr int (NPHIBINS|NCOSTHETABINS) = (\d+);", spectrum_lightcurve_h)
    }
    return bincounts["NPHIBINS"] * bincounts["NCOSTHETABINS"]


def merge_files(fileprefix: str, input_files: list[Path]) -> None:
    """Replace a merged file only after all input files copy successfully."""
    if not all(input_file.is_file() and input_file.stat().st_size > 0 for input_file in input_files):
        msg = f"Some {fileprefix} direction bins are missing or empty in {input_files[0].parent}"
        raise ValueError(msg)

    compressed = input_files[0].suffix == ".zst"
    if any((input_file.suffix == ".zst") != compressed for input_file in input_files):
        msg = f"The {fileprefix} direction bins in {input_files[0].parent} mix plain and .zst files"
        raise ValueError(msg)

    outfile = Path(f"{fileprefix}_res.out.zst" if compressed else f"{fileprefix}_res.out")
    otherform = Path(f"{fileprefix}_res.out" if compressed else f"{fileprefix}_res.out.zst")
    print(f"Merging {fileprefix} direction bins into {outfile}")
    # A sequence of zstd frames forms a valid zstd file. Copy the bytes unchanged.
    # The temporary folder shares the destination filesystem, so replace() is atomic.
    with tempfile.TemporaryDirectory(dir=Path(), prefix=f".merge_{fileprefix}_") as tempfolder:
        partialpath = Path(tempfolder) / outfile.name
        with partialpath.open("wb") as fout:
            for input_file in input_files:
                with input_file.open("rb") as fin:
                    shutil.copyfileobj(fin, fout)
        partialpath.replace(outfile)
    otherform.unlink(missing_ok=True)


def main() -> None:
    """Merge the direction bins and report a failure for each incomplete file set."""
    mabins = get_mabins()
    failures: list[str] = []
    for basefolder in [Path(), Path("speclc_angle_res")]:
        if not basefolder.is_dir():
            continue
        for fileprefix in ["light_curve", "spec", "specpol"]:
            input_files = [resfilepath(basefolder, fileprefix, abin) for abin in range(mabins)]
            # An absent output type is optional. Any present bin requires the complete set, including bin zero.
            if not any(input_file.is_file() for input_file in input_files):
                continue
            try:
                merge_files(fileprefix, input_files)
            except (OSError, ValueError) as err:
                print(str(err), file=sys.stderr)
                failures.append(f"{basefolder}/{fileprefix}")
    if failures:
        sys.exit(f"Could not merge these direction bins: {', '.join(failures)}")


if __name__ == "__main__":
    main()
