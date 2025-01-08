#!/bin/bash -x
#SBATCH --ntasks=1024
#SBATCH --ntasks-per-node=128
#SBATCH --time=24:00:00
#SBATCH --partition=k2-math-physics
#SBATCH --mail-type=ALL
##SBATCH --mail-user=fmcneill07@qub.ac.uk

module load libs/gcc/14.1.0
module load gsl/2.8/gcc-14.1.0
module load compilers/gcc/14.1.0
module load mpi/openmpi/5.0.3/gcc-14.1.0
module load apps/python3/3.12.4/gcc-14.1.0

module list

cd $SLURM_SUBMIT_DIR

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

hoursleft=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): before srun sn3d. hours left: $hoursleft"
# time srun --hint=nomultithread -- ./sn3d -w $hoursleft > out.txt
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
    sbatch ./artis/scripts/artis-kelvin2.sh
    # sbatch $SLURM_JOB_NAME
fi

if [ -f packets00_0000.out ]; then
    sbatch ./artis/scripts/exspec-gzip-kelvin2.sh
fi
