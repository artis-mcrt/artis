#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --account=p200999
#SBATCH --mail-type=ALL
##SBATCH --mail-user=luke.shingles@gmail.com

module load env/staging/2024.1 OpenMPI/5.0.3-GCC-13.3.0 zstd GSL Python/3.12.3-GCCcore-13.3.0
module list

cd $SLURM_SUBMIT_DIR

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

if [[ -f emission.out || -f emission.out.zst ]] && [[ -f exspec.txt ]]; then
  echo 'Not running exspec because emission.out[.zst] and exspec.txt were found'
else
  ./artis/scripts/exspec-before.sh
  ./artis/exspec
fi

./artis/scripts/exspec-after.sh
