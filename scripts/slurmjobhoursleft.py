#!/usr/bin/env python3

import datetime
import sys
import subprocess


def main() -> None:
    # usage:
    # python3 slurmjobhoursleft.py [JOBID]

    assert len(sys.argv) == 2
    jobid = int(sys.argv[1])
    cmd = f"squeue -j {jobid} --noheader --Format EndTime"
    cmdendtime = subprocess.run(cmd, capture_output=True, shell=True, check=True)
    strendtime = cmdendtime.stdout.decode().strip()

    endtime = datetime.datetime.fromisoformat(strendtime)
    total_sec = (endtime - datetime.datetime.now()).total_seconds()
    print(total_sec / 3600)


if __name__ == "__main__":
    main()
