#!/bin/bash

#PBS -P fm5
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=768GB
#PBS -l ncpus=496
#PBS -l wd
#PBS -m abe
#PBS -M luke.shingles@gmail.com

# ncpus must be a factor of the cores per node
# mem is total memory (all cores combined)

# on raijin normal queue Sandy Bridge:
#   16 cores per node
#     66% of nodes  32GB = 2GB per core
#     31% of nodes  64GB = 4GB per core
#      2% of nodes 128GB = 8GB per core

# raijin normalbw queue Broadwell:
#   28 cores per node
#     536 nodes 128GB = 4.57 GB per core
#     268 nodes 256GB = 9.14 GB per core
#      10 nodes   1TB = 36.57 GB per core
#
# raijin normal/normlbw walltime limits:
#    48 hours for 1-255 cores
#    24 hours for 256-511 cores
#    10 hours for 512-1023 cores
#     5 hours for >=1024 cores

module load gsl/2.5
module load intel-mpi/2018.3.222
module load intel-cc/2018.3.222

ulimit -l 2097152

mpirun ./sn3d > out.txt

mkdir ${PBS_JOBID}
./artis/scripts/movefiles.sh ${PBS_JOBID}

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    qsub $PBS_JOBNAME
fi
