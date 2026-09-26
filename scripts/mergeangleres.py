#!/usr/bin/env python3

import re
from pathlib import Path


def resfilepath(basefolder: Path | str, fileprefix: str, abin: int) -> Path:
    """Return the path of the file of one direction bin. The .zst file is the fallback of the plain file."""
    plainpath = Path(basefolder, f"{fileprefix}_res_{abin:02d}.out")
    zstpath = plainpath.with_name(plainpath.name + ".zst")
    return zstpath if not plainpath.is_file() and zstpath.is_file() else plainpath


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
                # The merged file gets the form of the input files: plain, or zstd compressed. A join of
                # zstd frames is a valid zstd file, so the bytes of the files go together unchanged.
                compressed = input_files[0].suffix == ".zst"
                if any((input_file.suffix == ".zst") != compressed for input_file in input_files):
                    print(f"The {fileprefix}_res_??.out files in {basefolder} mix plain and .zst files. No merge.")
                    continue
                Path(f"{fileprefix}_res.out").unlink(missing_ok=True)
                Path(f"{fileprefix}_res.out.zst").unlink(missing_ok=True)
                outfile = Path(f"{fileprefix}_res.out.zst" if compressed else f"{fileprefix}_res.out")
                print(f"Merging {fileprefix}_res_??.out into {outfile}")
                with outfile.open("wb") as fout:
                    for input_file in input_files:
                        fout.write(input_file.read_bytes())
            else:
                print(f"Some {fileprefix}_res_??.out files are missing or empty in {basefolder}")


if __name__ == "__main__":
    main()
