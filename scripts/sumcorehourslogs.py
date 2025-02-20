#!/usr/bin/env python3

from pathlib import Path


def main() -> None:
    sn3dlogfiles = sorted(Path().glob("**/output_0-0.txt"))

    last_line = ""
    total_core_hours = 0.0
    for logfile in sn3dlogfiles:
        prev_last_line = last_line
        with logfile.open("rt", encoding="utf-8") as flog:
            last_line = flog.readlines()[-1].strip()
        print(f"{str(logfile)}: ", end="")
        if prev_last_line != last_line:
            if "CPU hours" in last_line:
                job_core_hours = float(last_line.split("CPU hours")[0].split()[-1])
                total_core_hours += job_core_hours
                print(f"{job_core_hours:8.1f} core-h")
            elif "core hours" in last_line:
                job_core_hours = float(last_line.split("core hours")[0].split()[-1])
                total_core_hours += job_core_hours
                print(f"{job_core_hours:8.1f} core-h")
            else:
                print("  WARNING: sn3d didn't finish cleanly. Manually check log to get CPU time consumed.")
            print(f"  {last_line}")
        else:
            print("\n  ignored duplicate log")
        print()

    print(f"CPU time: {total_core_hours:12.3f}  core-h")
    print(f"CPU time: {total_core_hours / 1000:12.3f}  k core-h")


if __name__ == "__main__":
    main()
