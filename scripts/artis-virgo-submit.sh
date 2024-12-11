#!/bin/bash -x
export APPTAINER_SHARENS=true
export APPTAINER_CONFIGDIR=/tmp/$USER
sbatch -J $(basename $(exec pwd)) --ntasks=1920 --ntasks-per-node=128 --mem-per-cpu=2000MB --partition=long --time=24:00:00 --constraint=amd,epyc,7713 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh

# AMD EPYC 9654 with 192 real cores per node and 4 GB/core. march=znver4
#sbatch -J $(basename $(exec pwd)) --ntasks=960 --ntasks-per-node=192 --mem-per-cpu=2000MB --partition=new --time=24:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh
