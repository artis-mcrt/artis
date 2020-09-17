#!/bin/bash

#PBS -A dp033
#PBS -l walltime=24:00:00
#PBS -l nodes=27:ppn=36
##PBS -m abe
##PBS -M luke.shingles@gmail.com

# Make the Intel compiler and MPI libs available
module load gsl/intel/2.4
module load intel/compilers/18.0.3
module load intel/mpi/18.0.3

cd $PBS_O_WORKDIR

mpirun ./sn3d -w 24 > out.txt

mkdir ${PBS_JOBID}
./artis/scripts/movefiles.sh ${PBS_JOBID}

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    qsub $PBS_JOBNAME
fi

if [ -f packets00_0000.out ]; then
    qsub ./artis/scripts/exspec-gzip-dirac.sh
fi