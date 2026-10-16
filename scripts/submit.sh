#!/bin/bash
# Submit the ARTIS job script of the current cluster. Run this script in the model folder.
# The script finds the cluster from the fully qualified host name.

host_fqdn=$(hostname -f 2>/dev/null || hostname)

case "$host_fqdn" in
    *.cosma.*) cluster=cosma8 ;;
    *kelvin*) cluster=kelvin2 ;;
    *.leonardo.local) cluster=leonardo ;;
    ccsub* | lxbk* | *.gsi.de) cluster=virgo ;;
    *)
        echo "submit.sh: no cluster has the host name '$host_fqdn'. Submit the job script of your cluster with sbatch." >&2
        exit 1
        ;;
esac

if [ ! -d artis/scripts ]; then
    echo "submit.sh: the folder 'artis/scripts' is absent. Run this script in the model folder." >&2
    exit 1
fi

echo "submit.sh: host name '$host_fqdn' is on the cluster $cluster"

# Some clusters have a submit script, because sbatch must get the options on the command line.
if [ -f "artis/scripts/artis-$cluster-submit.sh" ]; then
    source "artis/scripts/artis-$cluster-submit.sh"
else
    sbatch -J "${PWD##*/}" ${EMAIL:+--mail-user="$EMAIL"} -- "artis/scripts/artis-$cluster.sh"
fi
