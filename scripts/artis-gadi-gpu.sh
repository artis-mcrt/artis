#!/bin/bash

#PBS -P fm5
#PBS -q gpuvolta
#PBS -l walltime=24:00:00
#PBS -l mem=8GB
#PBS -q gpuvolta
#PBS -l ngpus = 1, minimum ngpus request is 1.
#PBS -l ncpus = 12, minimum ncpus request is 12, in the multiple of 12, and 12 x ngpus
#PBS -l wd
#PBS -m abe
#PBS -M luke.shingles@gmail.com

# gadi gpuvolta queue
# 2x 24 core (Intel Xeon Cascade Lake Platinum 8268, 2.9 GHz) in 160 compute nodes
# 4 x Nvidia Tesla Volta V100 Accelerator on each node
# 384 GBytes of RAM on CPU
# Charge rate: 3 SU per resource-hour
# For this queue, SUs are based on the higher of CPU request, or memory request divided by 8GB
#  2x24 (48) cores per node
#  192GB per node (4GB per core)

module load gsl
module load intel-compiler
module unload intel-mpi
module load openmpi
module load cuda

# mpirun ./sn3d -w 24 > out.txt
./sn3d -w 24 > out.txt

mkdir ${PBS_JOBID}
./artis/scripts/movefiles.sh ${PBS_JOBID}

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    qsub $PBS_JOBNAME
fi

if [ -f packets00_0000.out ]; then
    qsub ./artis/scripts/exspec-gzip-gadi_raijin.sh
fi