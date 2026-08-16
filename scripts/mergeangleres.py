#!/usr/bin/env python3

import re
from pathlib import Path


def resfilepath(basefolder: Path | str, fileprefix: str, abin: int) -> Path:
    return Path(basefolder, f"{fileprefix}_res_{abin:02d}.out")


def get_mabins() -> int:
    # one output file per direction bin (MABINS = NPHIBINS * NCOSTHETABINS), with the expected count
    # parsed from the exspec.h of the source folder that this script belongs to
    exspec_h = (Path(__file__).resolve().parent.parent / "exspec.h").read_text(encoding="utf-8")
    bincounts = {name: int(value) for name, value in re.findall(r"constexpr int (NPHIBINS|NCOSTHETABINS) = (\d+);", exspec_h)}
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
                outfile = Path(f"{fileprefix}_res.out")
                outfile.unlink(missing_ok=True)
                Path(f"{fileprefix}_res.out.zst").unlink(missing_ok=True)
                print(f"Merging {fileprefix}_res_??.out into {outfile}")
                with outfile.open("wt", encoding="utf-8") as fout:
                    for input_file in input_files:
                        fout.writelines(input_file.open("rt", encoding="utf-8").readlines())
            else:
                print(f"Some {fileprefix}_res_??.out files are missing or empty in {basefolder}")


if __name__ == "__main__":
    main()
