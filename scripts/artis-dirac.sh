#!/bin/bash
#
#PBS -N artis
#PBS -j oe
#PBS -m bea
#PBS -M l.shingles@qub.ac.uk
#PBS -l walltime=24:00:00
#PBS -l nodes=1024:ppn=1
#PBS -A dp033

# Make the Intel compiler and MPI libs available
#module load intel/compilers/13.0.0
module load gsl
module load openmp
module load openmpi
#module load intel/impi/4.1.5

# Run the program
cd $PBS_O_WORKDIR
mpirun ./sn3d
