#!/usr/bin/env python3

from pathlib import Path


def resfilepath(fileprefix: str, abin: int) -> Path:
    return Path(f"{fileprefix}_res_{abin:02d}.out")


def main() -> None:
    for fileprefix in ["light_curve", "spec", "specpol"]:
        if not resfilepath(fileprefix, 0).is_file():
            continue

        outfile = Path(f"{fileprefix}_res.out")
        outfile.unlink(missing_ok=True)

        Path(f"{fileprefix}_res.out.xz").unlink(missing_ok=True)

        print(f"Merging {fileprefix}_res_??.out into {outfile}")
        with outfile.open("wt", encoding="utf-8") as fout:
            for abin in range(100):
                fout.writelines(resfilepath(fileprefix, abin).open("rt", encoding="utf-8").readlines())


if __name__ == "__main__":
    main()
