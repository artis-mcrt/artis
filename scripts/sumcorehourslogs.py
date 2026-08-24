#!/usr/bin/env python3

from datetime import datetime
from pathlib import Path


def get_job_tasks(loglines: list[str]) -> int | None:
    """Get the process count times the thread count from the lines of an sn3d log."""
    nprocs: int | None = None
    nthreads = 1
    for line in loglines:
        if "with nprocs " in line:
            nprocs = int(line.split("with nprocs ")[1].split()[0].rstrip(")"))
        elif "OpenMP parallelisation is active" in line:
            nthreads = int(line.split("(max ")[1].split(")")[0])
        elif "stdpar" in line and " CPU threads" in line:
            nthreads = int(line.split(" CPU threads")[0].split()[-1])
    return None if nprocs is None else nprocs * nthreads


def get_log_elapsed_hours(loglines: list[str]) -> float | None:
    """Get the elapsed hours between the first and the last timestamp of an sn3d log."""
    try:
        time_start = datetime.strptime(loglines[0].split(maxsplit=1)[0], "%Y-%m-%dT%H:%M:%SZ")
        time_end = datetime.strptime(loglines[-1].split(maxsplit=1)[0], "%Y-%m-%dT%H:%M:%SZ")
    except ValueError:
        return None
    return (time_end - time_start).total_seconds() / 3600.0


def get_slurm_elapsed_hours(logfile: Path) -> float | None:
    """Get the elapsed hours from the srun start to the slurm error, for the job of an sn3d log."""
    jobid = logfile.parent.name.removesuffix(".slurm")
    slurmoutfile = logfile.parent.parent / f"slurm-{jobid}.out"
    if not slurmoutfile.is_file():
        return None
    time_start = None
    time_error = None
    with slurmoutfile.open("rt", encoding="utf-8") as fslurmlog:
        for line in fslurmlog:
            if "before srun sn3d" in line:
                # example: "Sat Aug 15 21:28:07 CEST 2026: before srun sn3d. hours left: 47.97"
                # drop the timezone name, because the error time below is also a local time
                tokens = line.split(": before srun sn3d")[0].split()
                time_start = datetime.strptime(" ".join(tokens[:4] + tokens[5:]), "%a %b %d %H:%M:%S %Y")
            elif line.startswith("[") and "error:" in line:
                # example: "[2026-08-17T21:26:39.005] error: *** JOB 12401250 ... CANCELLED ..."
                time_error = datetime.fromisoformat(line[1:].split("]", maxsplit=1)[0])
    if time_start is None or time_error is None:
        return None
    return (time_error - time_start).total_seconds() / 3600.0


def main() -> None:
    sn3dlogfiles = sorted(Path().glob("**/output_0-0.txt"))
    col1width = max((len(str(logfile)) for logfile in sn3dlogfiles), default=0) + 1

    last_line = ""
    total_core_hours = 0.0
    ntasks: int | None = None
    for logfile in sn3dlogfiles:
        prev_last_line = last_line
        with logfile.open("rt", encoding="utf-8") as flog:
            loglines = flog.readlines()
        last_line = loglines[-1].strip()
        print(f"{str(logfile) + ':':{col1width}s} ", end="")
        if prev_last_line != last_line:
            if "CPU hours" in last_line:
                job_core_hours = float(last_line.split("CPU hours")[0].split()[-1])
                total_core_hours += job_core_hours
                print(f"{job_core_hours:8.1f} core-h")
            elif "core hours" in last_line:
                job_core_hours = float(last_line.split("core hours")[0].split()[-1])
                total_core_hours += job_core_hours
                print(f"{job_core_hours:8.1f} core-h")
                if " processes " in last_line and " threads " in last_line:
                    job_ntasks = int(last_line.split(" processes")[0].split()[-1]) * int(
                        last_line.split(" threads")[0].split()[-1]
                    )
                    # make sure number of CPUs is the same for all jobs
                    assert ntasks is None or ntasks == job_ntasks
                    ntasks = job_ntasks
            else:
                job_tasks = get_job_tasks(loglines)
                elapsed_hours = get_slurm_elapsed_hours(logfile)
                time_source = "the error time in the slurm log"
                if elapsed_hours is None:
                    elapsed_hours = get_log_elapsed_hours(loglines)
                    time_source = "the sn3d log timestamps"
                if job_tasks is not None and elapsed_hours is not None:
                    job_core_hours = elapsed_hours * job_tasks
                    total_core_hours += job_core_hours
                    print(f"{job_core_hours:8.1f} core-h")
                    print(f"  WARNING: sn3d did not finish cleanly. The value is an estimate from {time_source}.")
                else:
                    print("  WARNING: sn3d didn't finish cleanly. Manually check log to get CPU time consumed.")
            print(f"  {loglines[0].strip()}")
            print(f"  {last_line}")
        else:
            print("\n  ignored duplicate log")
        print()

    if ntasks is not None:
        print(f"{'Tasks:':15s} {ntasks:8d}")
        print(f"{'Wallclock time:':15s} {total_core_hours / ntasks:12.3f}  hours")

    print(f"{'CPU time:':15s} {total_core_hours:12.3f}  core-h")
    print(f"{'CPU time:':15s} {total_core_hours / 1000:12.3f}  k core-h")


if __name__ == "__main__":
    main()
