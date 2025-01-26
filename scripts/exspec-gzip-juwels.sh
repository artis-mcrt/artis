#!/bin/bash -x
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=batch
##SBATCH --account=rtsn22
#SBATCH --account=knrt23
#SBATCH --mail-type=ALL
##SBATCH --mail-user=luke.shingles@gmail.com

module load Stages/2025 ParaStationMPI GSL

cd $SLURM_SUBMIT_DIR

if [[ -f emission.out || -f emission.out.zst ]]; then
  ./artis/scripts/exspec-before.sh
  ./exspec
else
  echo 'Not running exspec because emission.out[.zst] was found'
fi

./artis/scripts/exspec-after.sh
