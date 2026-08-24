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

# GSL is no longer used by artis, but is still loaded so that these scripts work with older versions too
module load env/release/2025.1 gompi/2025a zstd GSL git Python
module list

cd $SLURM_SUBMIT_DIR

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

source ./artis/scripts/corehours-before.sh
echo "$(date): before exspec"

source ./artis/scripts/run-exspec-if-needed.sh

source ./artis/scripts/exspec-after.sh

echo "$(date): after exspec finished"
source ./artis/scripts/corehours-after.sh
