#!/bin/bash
#SBATCH --ntasks=960
##SBATCH --ntasks=1920
#SBATCH --ntasks-per-node=48
#SBATCH --time=24:00:00
#SBATCH --partition=batch
##SBATCH --partition=mem192
##SBATCH --account=rtsn22
#SBATCH --account=knrt23
#SBATCH --mail-type=ALL
##SBATCH --mail-user=luke.shingles@gmail.com

# GSL is no longer used by artis, but is still loaded so that these scripts work with older versions too
module load Stages/2025 GCC ParaStationMPI GSL
module load UCX-settings/plain

module list

cd $SLURM_SUBMIT_DIR

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

# decompress any zipped input files
source ./artis/scripts/exspec-before.sh

hoursleft=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
source ./artis/scripts/corehours-before.sh
echo "$(date): before srun sn3d. hours left: $hoursleft"
time srun --hint=nomultithread -- ./artis/sn3d -w $hoursleft -o ${SLURM_JOB_ID}.slurm > out.txt
hoursleftafter=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): after srun sn3d finished. hours left: $hoursleftafter"
source ./artis/scripts/corehours-after.sh

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    sbatch --job-name="$SLURM_JOB_NAME" ./artis/scripts/artis-juwels.sh
else
    # post-processing can remove restart files, so only queue it when no continuation job was submitted
    if [ -f packets00_0000.out ]; then
        sbatch --job-name="exspec_$SLURM_JOB_NAME" ./artis/scripts/exspec-zip-juwels.sh
    fi
fi
