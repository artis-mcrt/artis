#!/bin/bash -x
export APPTAINER_SHARENS=true
export APPTAINER_CONFIGDIR=/tmp/$USER
sbatch -J $(basename $(exec pwd)) --ntasks=1920 --ntasks-per-node=192 --mem-per-cpu=4000MB --partition=long --time=24:00:00 --constraint=amd,epyc,9654 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh

# AMD EPYC 9654 nodes have 192 real cores per node and 4 GB/core. march=znver4
# other nodes have 128 real cores per node
