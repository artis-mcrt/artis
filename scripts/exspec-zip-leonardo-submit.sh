#!/bin/bash -x
# Submit the exspec job on Leonardo. Use this script in place of a direct sbatch.
# See artis-leonardo-submit.sh for the reason of the --mail-user argument.

sbatch -J "exspec_${PWD##*/}" ${EMAIL:+--mail-user="$EMAIL"} -- artis/scripts/exspec-zip-leonardo.sh
