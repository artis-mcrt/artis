#!/bin/bash -x
# Submit the exspec job on Leonardo. Use this script in place of a direct sbatch.
# See artis-leonardo-submit.sh for the reason of the --mail-user argument.

# A job script sources this file, and a sourced file sees the arguments of the job script.
if [ "${BASH_SOURCE:-$0}" = "$0" ] && [ $# -gt 0 ]; then
    echo "Usage: artis/scripts/exspec-zip-leonardo-submit.sh with no argument, in the model folder. The script submits the exspec job on Leonardo with sbatch." >&2
    exit 1
fi

sbatch -J "exspec_${PWD##*/}" ${EMAIL:+--mail-user="$EMAIL"} -- artis/scripts/exspec-zip-leonardo-slurmjob.sh
