#!/usr/bin/env python3

import datetime
import subprocess
import sys


def main() -> None:
    # usage:
    # python3 slurmjobhoursleft.py [JOBID]

    assert len(sys.argv) == 2
    jobid = int(sys.argv[1])
    cmd = f"squeue -j {jobid} --noheader --Format EndTime"
    cmdendtime = subprocess.run(cmd, capture_output=True, shell=True, check=True, text=True)
    strendtime = cmdendtime.stdout.strip()

    endtime = datetime.datetime.fromisoformat(strendtime)
    # squeue prints a local time without a timezone, so the current time must also be naive
    total_sec = (endtime - datetime.datetime.now()).total_seconds()  # noqa: DTZ005
    print(total_sec / 3600)


if __name__ == "__main__":
    main()
