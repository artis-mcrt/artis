#!/bin/bash -x
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=mem192
#SBATCH --account=hmu14
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luke.shingles@gmail.com

module load Intel
module load ParaStationMPI
module load GSL

cd $SLURM_SUBMIT_DIR

./artis/scripts/exspec-before.sh

./exspec

./artis/scripts/exspec-after.sh
