#!/bin/bash -x
sbatch -J $(basename $(pwd)) --ntasks=960 --ntasks-per-node=192 --mem-per-cpu=2000MB --exclusive --partition=long --time=48:00:00 --constraint=amd,epyc,9654 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh

# AMD EPYC 9654 (Zen 4) nodes have 192 real cores per node
# AMD EPYC 9555 (Zen 5) nodes have 128 real cores per node
