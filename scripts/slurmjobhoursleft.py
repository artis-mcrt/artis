#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///
"""Print the hours that remain before the time limit of a Slurm job. Usage: slurmjobhoursleft.py JOBID."""

import datetime as dt
import subprocess
import sys


def main() -> None:
    """Print the hours between now and the end time that squeue gives for the job."""
    assert len(sys.argv) == 2
    jobid = int(sys.argv[1])
    cmd = ["squeue", "-j", str(jobid), "--noheader", "--Format", "EndTime"]
    cmdendtime = subprocess.run(cmd, capture_output=True, check=True, text=True)
    strendtime = cmdendtime.stdout.strip()

    # squeue prints a local time without a timezone. astimezone() attaches the local timezone to both times
    endtime = dt.datetime.fromisoformat(strendtime).astimezone()
    total_sec = (endtime - dt.datetime.now().astimezone()).total_seconds()
    print(total_sec / 3600)


if __name__ == "__main__":
    main()
