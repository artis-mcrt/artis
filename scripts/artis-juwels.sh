#!/bin/bash -x
#SBATCH --ntasks=480
#SBATCH --ntasks-per-node=48
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luke.shingles@gmail.com

module load Intel
module load Intel ParaStationMPI
module load GSL

cd $SLURM_SUBMIT_DIR

srun ./sn3d

mkdir slurm-${SLURM_JOBID}
./movefiles.sh slurm-${SLURM_JOBID}

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    sbatch $SLURM_JOB_NAME
fi

