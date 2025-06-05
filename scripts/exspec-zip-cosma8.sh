#!/bin/bash -l

#SBATCH --ntasks 128
#SBATCH -J exspec-zip
#SBATCH -p cosma8
#SBATCH -A dp033
#SBATCH --exclusive
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=f.callan@qub.ac.uk
##SBATCH --mail-user=luke.shingles@gmail.com

module purge
module load cosma
module load gnu_comp
module load gsl
module load openmpi
module load python

module list

cd $SLURM_SUBMIT_DIR

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

if [[ -f emission.out || -f emission.out.zst ]]; then
  echo 'Not running exspec because emission.out[.zst] was found'
else
  ./artis/scripts/exspec-before.sh
  ./exspec
fi

./artis/scripts/exspec-after.sh
