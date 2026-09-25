#!/usr/bin/env python3

import re
import subprocess
from pathlib import Path


def resfilepath(basefolder: Path | str, fileprefix: str, abin: int) -> Path:
    """Return the path of the file of one direction bin. The .zst file is the fallback of the plain file."""
    plainpath = Path(basefolder, f"{fileprefix}_res_{abin:02d}.out")
    zstpath = plainpath.with_name(plainpath.name + ".zst")
    return zstpath if not plainpath.is_file() and zstpath.is_file() else plainpath


def read_text(filepath: Path) -> str:
    """Return the content of a plain or a zstd compressed file."""
    if filepath.suffix == ".zst":
        return subprocess.run(["zstd", "-dc", str(filepath)], check=True, capture_output=True, text=True).stdout
    return filepath.read_text(encoding="utf-8")


def get_mabins() -> int:
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


def main() -> None:
    mabins = get_mabins()
    for basefolder in [Path(), Path("speclc_angle_res")]:
        if not basefolder.is_dir():
            continue
        for fileprefix in ["light_curve", "spec", "specpol"]:
            if not resfilepath(basefolder, fileprefix, 0).is_file():
                continue

            input_files = [resfilepath(basefolder, fileprefix, abin) for abin in range(mabins)]
            if all(input_file.is_file() and input_file.stat().st_size > 0 for input_file in input_files):
                # the merged file gets the form of the input files: plain, or zstd compressed
                compressed = input_files[0].suffix == ".zst"
                Path(f"{fileprefix}_res.out").unlink(missing_ok=True)
                Path(f"{fileprefix}_res.out.zst").unlink(missing_ok=True)
                outfile = Path(f"{fileprefix}_res.out.zst" if compressed else f"{fileprefix}_res.out")
                print(f"Merging {fileprefix}_res_??.out into {outfile}")
                merged = "".join(read_text(input_file) for input_file in input_files)
                if compressed:
                    subprocess.run(["zstd", "-q", "-13", "-o", str(outfile)], check=True, input=merged, text=True)
                else:
                    outfile.write_text(merged, encoding="utf-8")
            else:
                print(f"Some {fileprefix}_res_??.out files are missing or empty in {basefolder}")


if __name__ == "__main__":
    main()
