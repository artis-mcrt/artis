#!/usr/bin/env python3

import argparse
import json
import subprocess
from datetime import datetime
from pathlib import Path

try:
    import zstandard
except ModuleNotFoundError:
    zstandard = None


def get_job_procs_threads(loglines: list[str]) -> tuple[int, int] | None:
    """Get the process count and the thread count from the lines of an sn3d log."""
    nprocs: int | None = None
    nthreads = 1
    for line in loglines:
        if "with nprocs " in line:
            nprocs = int(line.split("with nprocs ")[1].split()[0].rstrip(")"))
        elif "OpenMP parallelisation is active" in line:
            nthreads = int(line.split("(max ")[1].split(")")[0])
        elif "stdpar" in line and " CPU threads" in line:
            nthreads = int(line.split(" CPU threads")[0].split()[-1])
    return None if nprocs is None else (nprocs, nthreads)


def get_log_elapsed_hours(loglines: list[str]) -> float | None:
    """Get the elapsed hours between the first and the last timestamp of an sn3d log."""
    try:
        time_start = datetime.strptime(loglines[0].split(maxsplit=1)[0], "%Y-%m-%dT%H:%M:%SZ")
        time_end = datetime.strptime(loglines[-1].split(maxsplit=1)[0], "%Y-%m-%dT%H:%M:%SZ")
    except ValueError:
        return None
    return (time_end - time_start).total_seconds() / 3600.0


def get_timestep_range(last_line: str) -> tuple[int, int] | None:
    """Get the first and the last timestep of the job from the final log line."""
    if "pktprop ts " not in last_line:
        return None
    # example: "... (job: pktprop ts 0 to ts 9 grid-preprop, 45.909 wallclock hours ..."
    tokens = last_line.split("pktprop ts ")[1].split()
    return (int(tokens[0]), int(tokens[3]))


def get_log_timestep_range(loglines: list[str]) -> tuple[int, int] | None:
    """Get the first and the last propagated timestep from the lines of an sn3d log."""
    ts_first: int | None = None
    ts_last: int | None = None
    for line in loglines:
        # example: "2026-08-20T21:31:02Z timestep 21: start update_packets"
        if ts_first is None and ": start update_packets" in line and " timestep " in line:
            ts_first = int(line.split(" timestep ")[1].split(":")[0])
        # example: "2026-08-22T03:26:29Z timestep 63: time after update grid on all processes (...)"
        if ": time after update grid" in line and " timestep " in line:
            ts_last = int(line.split(" timestep ")[1].split(":")[0])
    if ts_first is None or ts_last is None:
        return None
    return (ts_first, ts_last)


def read_loglines(logfile: Path) -> list[str]:
    """Read the lines of a log file. A file with the suffix .zst is decompressed first."""
    if logfile.suffix == ".zst":
        if zstandard is None:
            proc = subprocess.run(["zstd", "-dc", str(logfile)], check=True, capture_output=True, text=True)
            return proc.stdout.splitlines()
        with zstandard.open(logfile, "rt", encoding="utf-8") as flog:
            return flog.readlines()
    with logfile.open("rt", encoding="utf-8") as flog:
        return flog.readlines()


def main() -> None:
    parser = argparse.ArgumentParser(description="Sum the core hours of the sn3d logs in the given run folders.")
    parser.add_argument(
        "runfolders", nargs="*", type=Path, default=[Path()], help="run folders to scan (default: the current folder)"
    )
    parser.add_argument("--json", action="store_true", help="write the results as JSON instead of a table")
    args = parser.parse_args()

    sn3dlogfiles = sorted(
        logfile
        for runfolder in args.runfolders
        for pattern in ("**/output_0-0.txt", "**/output_0-0.txt.zst")
        for logfile in runfolder.glob(pattern)
    )
    col1width = max((len(str(logfile)) for logfile in sn3dlogfiles), default=0) + 1
    verbose = not args.json

    jobrows: list[dict[str, str | float | int | bool | None]] = []
    last_line = ""
    total_core_hours = 0.0
    ntasks: int | None = None
    timestep_ranges: list[tuple[int, int]] = []
    for logfile in sn3dlogfiles:
        prev_last_line = last_line
        loglines = read_loglines(logfile)
        last_line = loglines[-1].strip()
        jobrow: dict[str, str | float | int | bool | None] = {
            "logfile": str(logfile),
            "core_hours": None,
            "estimate": False,
            "duplicate": False,
        }
        jobrows.append(jobrow)
        if verbose:
            print(f"{str(logfile) + ':':{col1width}s} ", end="")
        if prev_last_line == last_line:
            jobrow["duplicate"] = True
            if verbose:
                print("\n  ignored duplicate log")
                print()
            continue

        job_core_hours: float | None = None
        estimate = False
        if "CPU hours" in last_line:
            job_core_hours = float(last_line.split("CPU hours")[0].split()[-1])
        elif "core hours" in last_line:
            job_core_hours = float(last_line.split("core hours")[0].split()[-1])
            if " processes " in last_line and " threads " in last_line:
                job_ntasks = int(last_line.split(" processes")[0].split()[-1]) * int(
                    last_line.split(" threads")[0].split()[-1]
                )
                # make sure number of CPUs is the same for all jobs
                assert ntasks is None or ntasks == job_ntasks
                ntasks = job_ntasks
        log_ts_range: tuple[int, int] | None = None
        estimate_line: str | None = None
        if job_core_hours is None:
            procs_threads = get_job_procs_threads(loglines)
            elapsed_hours = get_log_elapsed_hours(loglines)
            log_ts_range = get_log_timestep_range(loglines)
            if procs_threads is not None and elapsed_hours is not None:
                nprocs, nthreads = procs_threads
                job_core_hours = elapsed_hours * nprocs * nthreads
                estimate = True
                if log_ts_range is not None:
                    estimate_line = (
                        f"WARNING: sn3d did not finish. Estimated: "
                        f"(job: pktprop ts {log_ts_range[0]} to ts {log_ts_range[1]} grid-preprop, "
                        f"{elapsed_hours:.3f} wallclock hours * {nprocs} processes * {nthreads} threads "
                        f"= {job_core_hours:.3f} core hours)"
                    )

        ts_range = get_timestep_range(last_line)
        if ts_range is None:
            ts_range = log_ts_range
        if ts_range is not None:
            timestep_ranges.append(ts_range)
            jobrow["timestep_first"], jobrow["timestep_last"] = ts_range
        if job_core_hours is not None:
            total_core_hours += job_core_hours
            jobrow["core_hours"] = job_core_hours
            jobrow["estimate"] = estimate

        if verbose:
            if job_core_hours is not None:
                print(f"{job_core_hours:8.1f} core-h")
                if estimate and estimate_line is None:
                    print("  WARNING: sn3d did not finish cleanly. The value is an estimate from the log timestamps.")
            else:
                print("  WARNING: sn3d didn't finish cleanly. Manually check log to get CPU time consumed.")
            print(f"  {loglines[0].strip()}")
            print(f"  {last_line}")
            if estimate_line is not None:
                print(f"  {estimate_line}")
            print()

    timestep_ranges.sort()
    wallclock_hours = total_core_hours / ntasks if ntasks is not None else None

    if args.json:
        print(
            json.dumps(
                {
                    "jobs": jobrows,
                    "tasks": ntasks,
                    "wallclock_hours": wallclock_hours,
                    "total_core_hours": total_core_hours,
                    "timestep_ranges": [list(ts_range) for ts_range in timestep_ranges],
                },
                indent=2,
            )
        )
        return

    if ntasks is not None and wallclock_hours is not None:
        print(f"{'Tasks:':15s} {ntasks:8d}")
        print(f"{'Wallclock time:':15s} {wallclock_hours:12.3f}  hours")

    print(f"{'CPU time:':15s} {total_core_hours:12.3f}  core-h")
    print(f"{'CPU time:':15s} {total_core_hours / 1000:12.3f}  k core-h")


if __name__ == "__main__":
    main()
