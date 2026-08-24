#!/bin/bash -l

#SBATCH --ntasks 32
#SBATCH --time=48:00:00
#SBATCH --partition=cosma8-serial
#SBATCH --account=dp033
##SBATCH --account=dp385
#SBATCH --exclusive
#SBATCH -o slurm-%J.out
#SBATCH -e slurm-%J.out
#SBATCH --mail-type=ALL
##SBATCH --mail-user=f.callan@qub.ac.uk
##SBATCH --mail-user=luke.shingles@gmail.com

module purge
module load cosma
module load gnu_comp
# GSL is no longer used by artis, but is still loaded so that these scripts work with older versions too
module load gsl
module load openmpi
module load python

module list

export PATH=/cosma/local/intel/oneAPI_2021.3.0/intelpython/python3.7/pkgs/zstd-1.4.5-h2daa505_0/bin:$PATH

cd $SLURM_SUBMIT_DIR

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

source ./artis/scripts/corehours-before.sh
echo "$(date): before exspec"

source ./artis/scripts/run-exspec-if-needed.sh

source ./artis/scripts/exspec-after.sh

echo "$(date): after exspec finished"
source ./artis/scripts/corehours-after.sh
