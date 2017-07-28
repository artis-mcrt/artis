#!/bin/bash

#PBS -P fm5
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=1008GB
#PBS -l ncpus=1008
#PBS -l wd

# ncpus must be a factor of the cores per node (16)
# mem is total memory (all cores combined)
# on raijin, max 128GB/node = 8GB per core and MPI task (2% of nodes)
# max 64GB/node = 4GB per core (31% of nodes)
# max 32GB/node = 2GB per core (66% of nodes)
# normal queue walltime limits
# 48 hrs for 1-255 cores
# 24 hrs for 256-511 cores
# 10 hrs for 512-1023 cores
# 5 hours for 1024-56960 cores

ulimit -l 2097152

mpirun ./sn3d > out.txt

mkdir ${PBS_JOBID}
./movefiles.sh ${PBS_JOBID}

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    qsub artis-raijin.sh
fi
