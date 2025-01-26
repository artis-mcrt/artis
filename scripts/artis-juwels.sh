#!/bin/bash -x
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

module load Stages/2025 GCC ParaStationMPI GSL
module load UCX-settings/plain

module list

cd $SLURM_SUBMIT_DIR

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

hoursleft=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): before srun sn3d. hours left: $hoursleft"
time srun --hint=nomultithread -- ./sn3d -w $hoursleft > out.txt
hoursleftafter=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): after srun sn3d finished. hours left: $hoursleftafter"
hourselapsed=$(python3 -c "print($hoursleft - $hoursleftafter)")
echo "hours of runtime: $hourselapsed"
cpuhrs=$(python3 -c "print($SLURM_NTASKS * $hourselapsed)")
echo "ntasks: $SLURM_NTASKS -> CPU core hrs: $cpuhrs"

mkdir ${SLURM_JOB_ID}.slurm
./artis/scripts/movefiles.sh ${SLURM_JOB_ID}.slurm

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    sbatch --job-name="$SLURM_JOB_NAME" ./artis/scripts/artis-juwels.sh
fi

if [ -f packets00_0000.out ]; then
    sbatch --job-name="exspec_$SLURM_JOB_NAME" ./artis/scripts/exspec-gzip-juwels.sh
fi
