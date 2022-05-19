#!/bin/bash -x
#SBATCH --ntasks=960
#SBATCH --ntasks-per-node=48
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --account=rtsn22
##SBATCH --mail-type=ALL
##SBATCH --mail-user=luke.shingles@gmail.com

module load GCC
module load OpenMPI
module load GSL

cd $SLURM_SUBMIT_DIR

srun ./sn3d -w 24 > out.txt

mkdir ${SLURM_JOBID}.slurm
./artis/scripts/movefiles.sh ${SLURM_JOBID}.slurm

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    sbatch ./artis/scripts/artis-juwels.sh
    # sbatch $SLURM_JOB_NAME
fi

if [ -f packets00_0000.out ]; then
    sbatch ./artis/scripts/exspec-gzip-juwels.sh
fi
