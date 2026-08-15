#!/usr/bin/env python3

from pathlib import Path


def resfilepath(basefolder: Path | str, fileprefix: str, abin: int) -> Path:
    return Path(basefolder, f"{fileprefix}_res_{abin:02d}.out")


def main() -> None:
    for basefolder in [Path(), Path("speclc_angle_res")]:
        if not basefolder.is_dir():
            continue
        for fileprefix in ["light_curve", "spec", "specpol"]:
            if not resfilepath(basefolder, fileprefix, 0).is_file():
                continue

            # one file per direction bin (MABINS in exspec.h). Detect the bin count from the files on
            # disk rather than hardcoding it, so a changed bin count cannot silently break the merge
            files_present = sorted(basefolder.glob(f"{fileprefix}_res_*.out"))
            input_files = [resfilepath(basefolder, fileprefix, abin) for abin in range(len(files_present))]
            if files_present == input_files and all(input_file.stat().st_size > 0 for input_file in input_files):
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
