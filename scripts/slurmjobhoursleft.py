#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///
"""Print the hours that remain before the time limit of a Slurm job. Usage: slurmjobhoursleft.py JOBID."""

import datetime as dt
import subprocess
import sys


def get_hours_left(strendtime: str, now: dt.datetime) -> float:
    """Give the hours from now to the end time that squeue prints. now must be an aware local time."""
    # squeue prints a local time without a timezone. In the repeated hour of a daylight-saving change, such a time has
    # two meanings. The job still runs, so its end is the first meaning after now
    endtime_naive = dt.datetime.fromisoformat(strendtime)
    endtimes = sorted({endtime_naive.replace(fold=fold).astimezone() for fold in (0, 1)})
    endtime = next((endtime for endtime in endtimes if endtime > now), endtimes[-1])
    return (endtime - now).total_seconds() / 3600


def main() -> None:
    """Print the hours between now and the end time that squeue gives for the job."""
    assert len(sys.argv) == 2
    jobid = int(sys.argv[1])
    cmd = ["squeue", "-j", str(jobid), "--noheader", "--Format", "EndTime"]
    try:
        cmdendtime = subprocess.run(cmd, capture_output=True, check=True, text=True)
        print(get_hours_left(cmdendtime.stdout.strip(), dt.datetime.now().astimezone()))
    except (subprocess.CalledProcessError, ValueError) as err:
        # the job script then runs sn3d without a wall time limit
        sys.exit(f"slurmjobhoursleft.py: squeue gave no end time for job {jobid}: {err}")


if __name__ == "__main__":
    main()
