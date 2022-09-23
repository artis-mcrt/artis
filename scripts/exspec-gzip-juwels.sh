#!/bin/bash -x
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --account=rtsn22
##SBATCH --mail-type=ALL
##SBATCH --mail-user=luke.shingles@gmail.com

module load Stages/2022 Intel/2021.4.0 ParaStationMPI/5.5.0-1
module load mpi-settings/plain
module load GSL

cd $SLURM_SUBMIT_DIR

if [ ! -f emission.out.xz ]; then
  ./artis/scripts/exspec-before.sh
  ./exspec
else
  echo 'Not running exspec because emission.out.xz was found'
fi

./artis/scripts/exspec-after.sh
