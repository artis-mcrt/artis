#!/bin/bash -l

#SBATCH --ntasks 960
##SBATCH --ntasks 1024
##SBATCH --ntasks 1920
##SBATCH --ntasks 3072
#SBATCH -J artis
#SBATCH -p cosma8
#SBATCH -A dp033
#SBATCH --exclusive
#SBATCH -t 70:00:00
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

hoursleft=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): before srun sn3d. hours left: $hoursleft"
time mpirun -- ./sn3d -w $hoursleft > out.txt
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
    sbatch ./artis/scripts/artis-cosma8.sh
    # sbatch $SLURM_JOB_NAME
fi

if [ -f packets00_0000.out ]; then
    sbatch ./artis/scripts/exspec-zip-cosma8.sh
fi
