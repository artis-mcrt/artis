#!/usr/bin/env python3

import argparse
import json
from datetime import datetime
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Sum the core hours of the slurm logs in the given run folders.")
    parser.add_argument(
        "runfolders", nargs="*", type=Path, default=[Path()], help="run folders to scan (default: the current folder)"
    )
    parser.add_argument("--json", action="store_true", help="write the results as JSON instead of a table")
    args = parser.parse_args()

    slurmoutfiles = sorted(
        slurmoutfile for runfolder in args.runfolders for slurmoutfile in runfolder.glob("slurm-*.out")
    )
    jobs: list[dict[str, Path | str | float | int | datetime]] = [
        {"slurmoutfile": slurmoutfile, "jobid": str(slurmoutfile.name).removeprefix("slurm-").removesuffix(".out")}
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
                    # drop the timezone name, because the error time below is also a local time
                    if "before srun sn3d" in line:
                        jobdict["sn3d_started"] = True
                        datepart = line.split(": before srun sn3d")[0]
                    else:
                        jobdict["exspec_started"] = True
                        datepart = line.split(": before exspec")[0]
                    tokens = datepart.split()
                    jobdict["time_run_start"] = datetime.strptime(
                        " ".join(tokens[:4] + tokens[5:]), "%a %b %d %H:%M:%S %Y"
                    )
                if "after srun sn3d" in line or "after exspec finished" in line:
                    jobdict["run_finished"] = True
                if line.startswith("[") and "error:" in line:
                    # example: "[2026-08-17T21:26:39.005] error: *** JOB 12401250 ... CANCELLED ..."
                    jobdict["time_error"] = datetime.fromisoformat(line[1:].split("]", maxsplit=1)[0])
                if line.startswith(("ntasks:", "cpus-per-task:", "nodes:", "wallclock hrs:", "CPU core hrs:")):
                    var_vals.update(var_val.strip().split(": ", maxsplit=1) for var_val in line.split(" -> "))
            jobdict.update(var_vals)

    col1width = max((len(str(jobdict["slurmoutfile"])) for jobdict in jobs), default=6)

    # the shared core count of the sn3d jobs, as the fallback for a crashed job of an old job script
    ncores: int | None = None
    nnodes: int | None = None
    for jobdict in jobs:
        if "ntasks" in jobdict and jobdict.get("sn3d_started", False):
            job_ncores = int(str(jobdict["ntasks"])) * int(str(jobdict.get("cpus-per-task", "1")))
            # make sure number of CPUs is the same for all sn3d jobs
            assert ncores is None or ncores == job_ncores
            ncores = job_ncores
            if "nodes" in jobdict:
                job_nnodes = int(str(jobdict["nodes"]))
                assert nnodes is None or nnodes == job_nnodes
                nnodes = job_nnodes

    verbose = not args.json
    jobrows: list[dict[str, str | float | bool | None]] = []
    total_sn3d_core_hours = 0.0
    total_exspec_core_hours = 0.0
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

        job_core_hours: float | None = None
        estimate = False
        if "CPU core hrs" in jobdict:
            job_core_hours = float(str(jobdict["CPU core hrs"]))
        elif (
            (sn3d_started or exspec_started)
            and not run_finished
            and job_ncores is not None
            and isinstance(time_run_start, datetime)
            and isinstance(time_error, datetime)
        ):
            job_core_hours = (time_error - time_run_start).total_seconds() / 3600.0 * job_ncores
            estimate = True

        if job_core_hours is not None:
            if jobtype == "exspec":
                total_exspec_core_hours += job_core_hours
            else:
                total_sn3d_core_hours += job_core_hours
        jobrows.append(
            {
                "slurmoutfile": str(jobdict["slurmoutfile"]),
                "jobtype": jobtype,
                "core_hours": job_core_hours,
                "estimate": estimate,
            }
        )

        if verbose:
            print(f"{str(jobdict['slurmoutfile']):{col1width}s}  ", end="")
            if job_core_hours is not None:
                note = "  (exspec)" if jobtype == "exspec" else ""
                if estimate:
                    progname = "exspec" if jobtype == "exspec" else "sn3d"
                    note += f"  (WARNING: {progname} did not finish. Estimated from the error time.)"
                print(f"{job_core_hours:7.1f} core-h{note}")
            elif sn3d_started or exspec_started:
                progname = "exspec" if jobtype == "exspec" else "sn3d"
                print(f"{'?.?':>7s} core-h  (Unknown because {progname} started but didn't finish. Check output log.)")
            else:
                print(f"{'?.?':>7s} core-h  (Unknown because sn3d didn't start. exspec job?)")

    total_core_hours = total_sn3d_core_hours + total_exspec_core_hours
    wallclock_hours = total_sn3d_core_hours / ncores if ncores is not None else None
    node_hours = wallclock_hours * nnodes if wallclock_hours is not None and nnodes is not None else None

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
    if ncores is not None and wallclock_hours is not None:
        print(f"{'Tasks:':15s} {ncores:8d}")
        if nnodes is not None:
            print(f"{'Nodes:':15s} {nnodes:8d}")
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
