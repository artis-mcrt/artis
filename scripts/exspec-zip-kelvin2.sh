#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=k2-math-physics
#SBATCH --mail-type=ALL
##SBATCH --mail-user=fmcneill07@qub.ac.uk

module load libs/gcc/14.1.0
module load compilers/gcc/14.1.0
module load mpi/openmpi/5.0.3/gcc-14.1.0
module load apps/python3/3.12.4/gcc-14.1.0

module list

cd "$SLURM_SUBMIT_DIR" || exit

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

source ./artis/scripts/run-exspec-if-needed.sh

source ./artis/scripts/exspec-after.sh
