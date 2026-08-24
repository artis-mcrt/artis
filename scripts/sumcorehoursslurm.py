#!/usr/bin/env python3

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError


def main() -> None:
    parser = argparse.ArgumentParser(description="Sum the core hours of the slurm logs in the given run folders.")
    parser.add_argument(
        "runfolders",
        nargs="*",
        type=Path,
        default=[Path()],
        help="run folders to scan (default: the current folder)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="write the results as JSON instead of a table",
    )
    args = parser.parse_args()

    slurmoutfiles = sorted(
        slurmoutfile for runfolder in args.runfolders for slurmoutfile in runfolder.glob("slurm-*.out")
    )
    jobs: list[dict[str, Path | str | float | int | datetime]] = [
        {
            "slurmoutfile": slurmoutfile,
            "jobid": str(slurmoutfile.name).removeprefix("slurm-").removesuffix(".out"),
        }
        for slurmoutfile in slurmoutfiles
    ]

    for jobdict in jobs:
        outfilepath = jobdict["slurmoutfile"]
        assert isinstance(outfilepath, Path)
        with outfilepath.open("rt", encoding="utf-8") as fslurmlog:
            var_vals: dict[str, str] = {}
            for line in fslurmlog:
                if "before srun sn3d" in line or "before exspec" in line:
                    # example: "Sat Aug 15 21:28:07 CEST 2026: before srun sn3d. hours left: 47.97"
                    # drop the timezone name. The error time below is also a local time.
                    if "before srun sn3d" in line:
                        jobdict["sn3d_started"] = True
                        datepart = line.split(": before srun sn3d")[0]
                    else:
                        jobdict["exspec_started"] = True
                        datepart = line.split(": before exspec")[0]
                    tokens = datepart.split()
                    try:
                        jobdict["time_run_start"] = datetime.strptime(  # noqa: DTZ007
                            " ".join(tokens[:4] + tokens[5:]), "%a %b %d %H:%M:%S %Y"
                        )
                        # the abbreviation resolves an ambiguous time in the repeated hour
                        jobdict["tzabbrev"] = tokens[4]
                    except ValueError:
                        # a shell trace line or a different date format gives no start time
                        pass
                if "after srun sn3d" in line or "after exspec finished" in line:
                    jobdict["run_finished"] = True
                if "error:" in line:
                    # example: "[2026-08-17T21:26:39.005] error: *** JOB 12401250 ... CANCELLED AT 2026-08-17T21:26:39 ..."
                    # the bracketed timestamp is absent on some clusters, so also read the CANCELLED AT time
                    timestr = None
                    if line.startswith("["):
                        timestr = line[1:].split("]", maxsplit=1)[0]
                    elif " CANCELLED AT " in line:
                        timestr = line.split(" CANCELLED AT ")[1].split()[0]
                    if timestr is not None:
                        try:
                            jobdict["time_error"] = datetime.fromisoformat(timestr)
                        except ValueError:
                            # a bracketed line from another tool has no timestamp
                            pass
                if line.startswith(
                    (
                        "ntasks:",
                        "cpus-per-task:",
                        "nodes:",
                        "timezone:",
                        "wallclock hrs:",
                        "CPU core hrs:",
                    )
                ):
                    for var_val in line.split(" -> "):
                        var, _, val = var_val.strip().partition(": ")
                        if val:
                            var_vals[var] = val
            jobdict.update(var_vals)

    col1width = max((len(str(jobdict["slurmoutfile"])) for jobdict in jobs), default=6)

    # the shared core count of the sn3d jobs, as the fallback for a crashed job of an old job script
    ncores_seen: set[int] = set()
    nnodes_seen: set[int] = set()
    nnodes_all_known = True
    for jobdict in jobs:
        if not jobdict.get("sn3d_started", False):
            continue
        if "ntasks" in jobdict:
            ncores_seen.add(int(str(jobdict["ntasks"])) * int(str(jobdict.get("cpus-per-task", "1"))))
        if "nodes" in jobdict:
            nnodes_seen.add(int(str(jobdict["nodes"])))
        else:
            # an old log has no node count, so there is no shared node count
            nnodes_all_known = False
    # with a mix of core counts, the summary has no single task count and omits the wallclock time
    ncores = next(iter(ncores_seen)) if len(ncores_seen) == 1 else None
    nnodes = next(iter(nnodes_seen)) if nnodes_all_known and len(nnodes_seen) == 1 else None

    verbose = not args.json
    jobrows: list[dict[str, str | float | bool | None]] = []
    total_sn3d_core_hours = 0.0
    total_exspec_core_hours = 0.0
    total_wallclock_hours = 0.0
    wallclock_complete = True
    total_node_hours = 0.0
    node_hours_complete = True
    for jobdict in jobs:
        sn3d_started = bool(jobdict.get("sn3d_started", False))
        exspec_started = bool(jobdict.get("exspec_started", False))
        jobtype = "sn3d" if sn3d_started else ("exspec" if exspec_started else "unknown")
        if "ntasks" in jobdict:
            job_ncores = int(str(jobdict["ntasks"])) * int(str(jobdict.get("cpus-per-task", "1")))
        else:
            # a slurm log of an old job script logs ntasks only after a clean finish
            job_ncores = ncores
        run_finished = bool(jobdict.get("run_finished", False))
        time_run_start = jobdict.get("time_run_start")
        time_error = jobdict.get("time_error")
        if isinstance(time_run_start, datetime) and isinstance(time_error, datetime):
            # attach the logged timezone to a naive time.
            # an aware pair then subtracts correctly across a daylight-saving change.
            tzname = jobdict.get("timezone")
            if isinstance(tzname, str):
                try:
                    tzone = ZoneInfo(tzname)
                    if time_run_start.tzinfo is None:
                        time_run_start = time_run_start.replace(tzinfo=tzone)
                        # a start in the repeated hour of a daylight-saving change is ambiguous.
                        # the abbreviation from the date line selects the correct fold.
                        tzabbrev = jobdict.get("tzabbrev")
                        if isinstance(tzabbrev, str) and time_run_start.tzname() != tzabbrev:
                            time_fold1 = time_run_start.replace(fold=1)
                            if time_fold1.tzname() == tzabbrev:
                                time_run_start = time_fold1
                    if time_error.tzinfo is None:
                        time_error = time_error.replace(tzinfo=tzone)
                except (ZoneInfoNotFoundError, ValueError):
                    # an unknown timezone name keeps the times naive
                    pass
            if (time_run_start.tzinfo is None) != (time_error.tzinfo is None):
                # a naive time and an aware time have no valid difference
                time_error = None
            elif time_run_start.tzinfo is not None and time_error.tzinfo is not None:
                # convert to UTC. A subtraction of two times with one shared tzinfo ignores the offsets.
                time_run_start = time_run_start.astimezone(UTC)
                time_error = time_error.astimezone(UTC)

        job_core_hours: float | None = None
        estimate = False
        if "CPU core hrs" in jobdict:
            try:
                job_core_hours = float(str(jobdict["CPU core hrs"]))
            except ValueError:
                job_core_hours = None
        if (
            job_core_hours is None
            and (sn3d_started or exspec_started)
            and not run_finished
            and job_ncores is not None
            and isinstance(time_run_start, datetime)
            and isinstance(time_error, datetime)
            and time_error > time_run_start
        ):
            # slurm can add a grace period after the time limit, so the elapsed time has no upper clamp
            elapsed_hours = (time_error - time_run_start).total_seconds() / 3600.0
            job_core_hours = elapsed_hours * job_ncores
            estimate = True

        if job_core_hours is not None:
            if jobtype == "exspec":
                total_exspec_core_hours += job_core_hours
            else:
                total_sn3d_core_hours += job_core_hours
                # the per-job wallclock and node hours also work with a mix of allocations
                if job_ncores is not None:
                    job_wallclock_hours = job_core_hours / job_ncores
                    total_wallclock_hours += job_wallclock_hours
                    if "nodes" in jobdict:
                        total_node_hours += job_wallclock_hours * int(str(jobdict["nodes"]))
                    else:
                        node_hours_complete = False
                else:
                    wallclock_complete = False
                    node_hours_complete = False
        elif sn3d_started:
            # a started sn3d job without a time value makes the wallclock sum incomplete
            wallclock_complete = False
            node_hours_complete = False
        jobrows.append(
            {
                "slurmoutfile": str(jobdict["slurmoutfile"]),
                "jobtype": jobtype,
                "core_hours": job_core_hours,
                "estimate": estimate,
            }
        )

        if verbose:
            progname = "exspec" if jobtype == "exspec" else "sn3d"
            print(f"{jobdict['slurmoutfile']!s:{col1width}s}  ", end="")
            if job_core_hours is not None:
                note = "  (exspec)" if jobtype == "exspec" else ""
                if estimate:
                    note += f"  (WARNING: {progname} did not finish. Estimated from the error time.)"
                print(f"{job_core_hours:7.1f} core-h{note}")
            elif sn3d_started or exspec_started:
                print(f"{'?.?':>7s} core-h  (Unknown because {progname} started but did not finish. Check the log.)")
            else:
                print(f"{'?.?':>7s} core-h  (Unknown because sn3d did not start. exspec job?)")

    total_core_hours = total_sn3d_core_hours + total_exspec_core_hours
    wallclock_hours = total_wallclock_hours if wallclock_complete and total_sn3d_core_hours > 0.0 else None
    node_hours = total_node_hours if node_hours_complete and wallclock_hours is not None else None

    if args.json:
        print(
            json.dumps(
                {
                    "jobs": jobrows,
                    "tasks": ncores,
                    "nodes": nnodes,
                    "wallclock_hours": wallclock_hours,
                    "node_hours": node_hours,
                    "sn3d_core_hours": total_sn3d_core_hours,
                    "exspec_core_hours": total_exspec_core_hours,
                    "total_core_hours": total_core_hours,
                },
                indent=2,
            )
        )
        return

    print()
    print(f"{'Total:':{col1width}s}  {total_core_hours:7.1f} core-h")
    print()
    if len(ncores_seen) > 1:
        print("WARNING: the sn3d jobs use different core counts. The summary omits the task count.")
    if ncores is not None:
        print(f"{'Tasks:':15s} {ncores:8d}")
    if nnodes is not None:
        print(f"{'Nodes:':15s} {nnodes:8d}")
    if wallclock_hours is not None:
        print(f"{'Wallclock time:':15s} {wallclock_hours:12.3f}  hours")
        if node_hours is not None:
            print(f"{'Node time:':15s} {node_hours:12.3f}  node-h")

    if total_exspec_core_hours > 0.0:
        print(f"{'sn3d time:':15s} {total_sn3d_core_hours:12.3f}  core-h")
        print(f"{'exspec time:':15s} {total_exspec_core_hours:12.3f}  core-h")

    print(f"{'CPU time:':15s} {total_core_hours:12.3f}  core-h")
    print(f"{'CPU time:':15s} {total_core_hours / 1000:12.3f}  k core-h")


if __name__ == "__main__":
    main()
