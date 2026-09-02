#!/bin/bash -x
sbatch -J $(basename $(pwd)) --ntasks=960 --ntasks-per-node=176 --mem-per-cpu=2000MB --partition=long --time=48:00:00 --constraint=amd,epyc,9654 --mail-type=ALL --mail-user=${EMAIL} --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh
# AMD EPYC 9654 (Zen 4) nodes have 192 real cores per node (use 192 and --exclusive for benchmarking)

#sbatch -J $(basename $(pwd)) --ntasks=1024 --ntasks-per-node=128 --mem-per-cpu=2000MB --exclusive --partition=long --time=48:00:00 --constraint=amd,epyc,9555 --mail-type=ALL --mail-user=${EMAIL} --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh
# AMD EPYC 9555 (Zen 5) nodes have 128 real cores per node

# you should add a line like this to your .bashrc to set the EMAIL env variable for notifications:
# export EMAIL=your_email_address
