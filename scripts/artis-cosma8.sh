#!/bin/bash -l

#SBATCH --ntasks 960
##SBATCH --ntasks 1024
##SBATCH --ntasks 1920
##SBATCH --ntasks 3072
#SBATCH --time=70:00:00
#SBATCH --partition=cosma8-rome
##SBATCH --partition=cosma8-milan
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

cd artis
make sn3d
cd ..

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

# decompress any zipped input files
source ./artis/scripts/exspec-before.sh

hoursleft=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): before srun sn3d. hours left: $hoursleft"
time mpirun -- ./artis/sn3d -w $hoursleft -o ${SLURM_JOB_ID}.slurm > out.txt
hoursleftafter=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): after srun sn3d finished. hours left: $hoursleftafter"
hourselapsed=$(python3 -c "print($hoursleft - $hoursleftafter)")
echo "hours of runtime: $hourselapsed"
cpuhrs=$(python3 -c "print($SLURM_NTASKS * $hourselapsed)")
echo "ntasks: $SLURM_NTASKS -> CPU core hrs: $cpuhrs"

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    sbatch -J artis_$(basename $(pwd)) ./artis/scripts/artis-cosma8.sh
    # sbatch $SLURM_JOB_NAME
fi

if [ -f packets00_0000.out ]; then
    sbatch -J exspec_$(basename $(pwd)) ./artis/scripts/exspec-zip-cosma8.sh
fi
