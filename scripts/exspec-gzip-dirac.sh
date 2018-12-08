#!/bin/bash

#PBS -A dp033
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1

module load gsl/intel/2.4
module load intel/compilers/18.0.3
module load intel/mpi/18.0.3

cd $PBS_O_WORKDIR

./artis/scripts/exspec-before.sh

./exspec

./artis/scripts/exspec-after.sh
