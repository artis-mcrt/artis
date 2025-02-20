#!/usr/bin/env python3

from pathlib import Path


def main() -> None:
    slurmoutfiles = sorted(Path().glob("slurm-*.out"))
    jobs: list[dict[str, Path | str | float | int]] = [
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
                if "after srun sn3d" in line:
                    jobdict["sn3d_finished"] = True
                if line.startswith("ntasks:"):
                    var_vals = dict(
                        var_val.strip().split(": ", maxsplit=1) for var_val in line.split(" -> ", maxsplit=1)
                    )
            jobdict.update(var_vals)

    ntasks: str | None = None
    total_core_hours = 0.0
    for jobdict in jobs:
        print(f"{str(jobdict['slurmoutfile']):20s} ", end="")
        if "ntasks" in jobdict:
            # make sure number of CPUs is the same for all jobs
            assert ntasks is None or ntasks == jobdict["ntasks"]
            assert isinstance(jobdict["ntasks"], str)
            ntasks = jobdict["ntasks"]
        if "CPU core hrs" in jobdict:
            assert isinstance(jobdict["CPU core hrs"], str)
            job_core_hours = float(jobdict["CPU core hrs"])
            total_core_hours += job_core_hours
            print(f"{job_core_hours:7.1f} core-h")
        else:
            sn3d_started = jobdict.get("sn3d_started", False)
            sn3d_finished = jobdict.get("sn3d_finished", False)
            if sn3d_started and not sn3d_finished:
                print(f"{'?.?':>7s} core-h  (Unknown because sn3d started but didn't finish. Check output log.)")
            else:
                print(f"{'?.?':>7s} core-h  (Unknown because sn3d didn't start. exspec job?)")

    print()
    print(f"{'Total:':20s} {total_core_hours:7.1f} core-h")
    print()
    if ntasks is not None:
        print(f"CPUs: {ntasks}")
        print(f"Wallclock time: {total_core_hours / float(ntasks):.1f} hours")

    print(f"CPU time: {total_core_hours:12.3f}  core-h")
    print(f"CPU time: {total_core_hours / 1000:12.3f}  k core-h")


if __name__ == "__main__":
    main()
