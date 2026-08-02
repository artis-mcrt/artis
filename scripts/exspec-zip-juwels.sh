#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=batch
##SBATCH --account=rtsn22
#SBATCH --account=knrt23
#SBATCH --mail-type=ALL
##SBATCH --mail-user=luke.shingles@gmail.com

module load Stages/2025
module load ParaStationMPI
# GSL is no longer used by artis, but is still loaded so that these scripts work with older versions too
module load GSL
module load zstd/.1.5.6
module list

cd $SLURM_SUBMIT_DIR

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

source ./artis/scripts/run-exspec-if-needed.sh

source ./artis/scripts/exspec-after.sh
