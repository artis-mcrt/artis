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

echo "ntasks: ${SLURM_NTASKS:-1}"
echo "cpus-per-task: ${SLURM_CPUS_PER_TASK:-1}"
echo "nodes: ${SLURM_JOB_NUM_NODES:-1}"
starttime=$(date +%s)
echo "$(date): before exspec"

source ./artis/scripts/run-exspec-if-needed.sh

source ./artis/scripts/exspec-after.sh

echo "$(date): after exspec finished"
hourselapsed=$(awk "BEGIN{print ($(date +%s) - $starttime) / 3600}")
echo "wallclock hrs: $hourselapsed"
cpuhrs=$(awk "BEGIN{print ${SLURM_NTASKS:-1} * ${SLURM_CPUS_PER_TASK:-1} * $hourselapsed}")
echo "CPU core hrs: $cpuhrs"
