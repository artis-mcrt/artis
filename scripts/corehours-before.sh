#!/usr/bin/env bash

# the job scripts source this file before the run
# it logs the CPU allocation and it records the start time

echo "ntasks: ${SLURM_NTASKS:-1}"
echo "cpus-per-task: ${SLURM_CPUS_PER_TASK:-1}"
echo "nodes: ${SLURM_JOB_NUM_NODES:-1}"
echo "timezone: $(readlink /etc/localtime 2>/dev/null | sed 's|.*zoneinfo/||')"
corehours_starttime=$(date +%s)
