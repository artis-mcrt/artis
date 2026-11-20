#!/bin/bash -x

# A job script sources this file, and a sourced file sees the arguments of the job script.
if [ "${BASH_SOURCE:-$0}" = "$0" ] && [ $# -gt 0 ]; then
    echo "Usage: artis/scripts/exspec-zip-virgo-submit.sh with no argument, in the model folder. The script submits the exspec job on Virgo with sbatch." >&2
    exit 1
fi
sbatch -J "exspec_${PWD##*/}" --ntasks=1 --cpus-per-task 64 --mem-per-cpu=4000MB --partition=long --time=48:00:00 --constraint="[9654|9555]" --mail-type=ALL --mail-user="${EMAIL:-${USER}@gsi.de}" -- artis/scripts/exspec-zip-virgo-slurmjob.sh
