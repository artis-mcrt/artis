#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# the memory of a node (1030000 MB) divided by its 128 CPUs
#SBATCH --mem-per-cpu=8046M
#SBATCH --time=48:00:00
#SBATCH --partition=cosma8-serial
#SBATCH --account=dp033
##SBATCH --account=dp385
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

cd "${SLURM_SUBMIT_DIR:?}" || exit 1

# Use the exspec that the sn3d job script built from the same source as sn3d. A build here could read
# a newer checkout.

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

source ./artis/scripts/corehours-before.sh
echo "$(date): before exspec"

source ./artis/scripts/run-exspec-if-needed.sh

source ./artis/scripts/exspec-after.sh

echo "$(date): after exspec finished"
source ./artis/scripts/corehours-after.sh
