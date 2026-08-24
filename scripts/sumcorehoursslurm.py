#!/usr/bin/env python3

from datetime import datetime
from pathlib import Path


def main() -> None:
    slurmoutfiles = sorted(Path().glob("slurm-*.out"))
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
                if "before srun sn3d" in line:
                    jobdict["sn3d_started"] = True
                    # example: "Sat Aug 15 21:28:07 CEST 2026: before srun sn3d. hours left: 47.97"
                    # drop the timezone name, because the error time below is also a local time
                    tokens = line.split(": before srun sn3d")[0].split()
                    jobdict["time_srun_start"] = datetime.strptime(
                        " ".join(tokens[:4] + tokens[5:]), "%a %b %d %H:%M:%S %Y"
                    )
                if "after srun sn3d" in line:
                    jobdict["sn3d_finished"] = True
                if line.startswith("[") and "error:" in line:
                    # example: "[2026-08-17T21:26:39.005] error: *** JOB 12401250 ... CANCELLED ..."
                    jobdict["time_error"] = datetime.fromisoformat(line[1:].split("]", maxsplit=1)[0])
                if line.startswith("ntasks:"):
                    var_vals = dict(
                        var_val.strip().split(": ", maxsplit=1) for var_val in line.split(" -> ", maxsplit=1)
                    )
            jobdict.update(var_vals)

    col1width = max((len(str(jobdict["slurmoutfile"])) for jobdict in jobs), default=6)

    ntasks: str | None = None
    for jobdict in jobs:
        if "ntasks" in jobdict:
            # make sure number of CPUs is the same for all jobs
            assert ntasks is None or ntasks == jobdict["ntasks"]
            assert isinstance(jobdict["ntasks"], str)
            ntasks = jobdict["ntasks"]

    total_core_hours = 0.0
    for jobdict in jobs:
        print(f"{str(jobdict['slurmoutfile']):{col1width}s} ", end="")
        if "CPU core hrs" in jobdict:
            assert isinstance(jobdict["CPU core hrs"], str)
            job_core_hours = float(jobdict["CPU core hrs"])
            total_core_hours += job_core_hours
            print(f"{job_core_hours:7.1f} core-h")
        else:
            sn3d_started = jobdict.get("sn3d_started", False)
            sn3d_finished = jobdict.get("sn3d_finished", False)
            time_srun_start = jobdict.get("time_srun_start")
            time_error = jobdict.get("time_error")
            if (
                sn3d_started
                and not sn3d_finished
                and ntasks is not None
                and isinstance(time_srun_start, datetime)
                and isinstance(time_error, datetime)
            ):
                job_core_hours = (time_error - time_srun_start).total_seconds() / 3600.0 * float(ntasks)
                total_core_hours += job_core_hours
                print(
                    f"{job_core_hours:7.1f} core-h  (WARNING: sn3d did not finish. Estimated from the error time.)"
                )
            elif sn3d_started and not sn3d_finished:
                print(f"{'?.?':>7s} core-h  (Unknown because sn3d started but didn't finish. Check output log.)")
            else:
                print(f"{'?.?':>7s} core-h  (Unknown because sn3d didn't start. exspec job?)")

    print()
    print(f"{'Total:':{col1width}s} {total_core_hours:7.1f} core-h")
    print()
    if ntasks is not None:
        print(f"{'Tasks:':15s} {int(ntasks):8d}")
        print(f"{'Wallclock time:':15s} {total_core_hours / float(ntasks):12.3f}  hours")

    print(f"{'CPU time:':15s} {total_core_hours:12.3f}  core-h")
    print(f"{'CPU time:':15s} {total_core_hours / 1000:12.3f}  k core-h")


if __name__ == "__main__":
    main()
