#!/bin/bash -x
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --account=hmu14
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luke.shingles@gmail.com

module load Intel
module load ParaStationMPI
module load GSL

cd $SLURM_SUBMIT_DIR

if [ ! -f spec.out ]; then
  ./artis/scripts/exspec-before.sh
  ./exspec
else
  echo 'Not running exspec because spec.out was found'
fi

./artis/scripts/exspec-after.sh
