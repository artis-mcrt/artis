#!/bin/bash

#PBS -P fm5
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=1920GB
#PBS -l ncpus=960
#PBS -l wd
##PBS -m abe
##PBS -M luke.shingles@gmail.com

# ncpus must be a factor of the cores per node
# mem is total memory (all cores combined)

# gadi normal queue Cascade Lake
#  2x24 (48) cores per node
#  192GB per node (4GB per core)

# gadi normal queue
# Charge rate: 2 SU per resource-hour (walltime)
# "one resource" is the higher of CPU request, or memory request divided by 4GB

module load gsl
module load intel-compiler
module unload intel-mpi
module load openmpi

mpirun ./sn3d -w 24 > out.txt

mkdir ${PBS_JOBID}
./artis/scripts/movefiles.sh ${PBS_JOBID}

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    qsub $PBS_JOBNAME
fi

if [ -f packets00_0000.out ]; then
    qsub ./artis/scripts/exspec-gzip-gadi_raijin.sh
fi