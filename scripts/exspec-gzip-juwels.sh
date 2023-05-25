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

module load Stages/2023 GCC ParaStationMPI GSL

export PATH="/p/software/juwels/stages/2022/software/zstd/1.5.0-GCCcore-11.2.0/bin/:$PATH"

cd $SLURM_SUBMIT_DIR

if [ ! -f emission.out.zst ]; then
  ./artis/scripts/exspec-before.sh
  ./exspec
else
  echo 'Not running exspec because emission.out.zst was found'
fi

./artis/scripts/exspec-after.sh
