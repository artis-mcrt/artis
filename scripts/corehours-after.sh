#!/usr/bin/env bash

# the job scripts source this file after the run
# it logs the wallclock hours and the core hours, which sumcorehoursslurm.py reads

hourselapsed=$(awk "BEGIN{print ($(date +%s) - ${corehours_starttime:-$(date +%s)}) / 3600}")
echo "wallclock hrs: $hourselapsed"
cpuhrs=$(awk "BEGIN{print ${SLURM_NTASKS:-1} * ${SLURM_CPUS_PER_TASK:-1} * $hourselapsed}")
echo "CPU core hrs: $cpuhrs"
