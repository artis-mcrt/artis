#!/usr/bin/env bash

# the job scripts source this file before the run
# it logs the CPU allocation and it records the start time

echo "ntasks: ${SLURM_NTASKS:-1}"
echo "cpus-per-task: ${SLURM_CPUS_PER_TASK:-1}"
echo "nodes: ${SLURM_JOB_NUM_NODES:-1}"
# /etc/localtime is usually a symlink into the zoneinfo folder. On a host where it is
# a regular file, /etc/timezone or TZ can hold the zone name instead.
tzname=$(readlink /etc/localtime 2>/dev/null | sed 's|.*zoneinfo/||')
if [ -z "$tzname" ]; then
    tzname=$(cat /etc/timezone 2>/dev/null)
fi
echo "timezone: ${tzname:-$TZ}"
corehours_starttime=$(date +%s)
