#!/bin/bash

#PBS -A dp033
#PBS -l walltime=24:00:00
#PBS -l nodes=13:ppn=36

# Make the Intel compiler and MPI libs available
module load gsl/intel/2.4
module load intel/compilers/18.0.1
module load intel/mpi/18.0.1

# Run the program
cd $PBS_O_WORKDIR
mpirun ./sn3d > out.txt

mkdir ${PBS_JOBID}
./movefiles.sh ${PBS_JOBID}

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    qsub $PBS_JOBNAME
fi