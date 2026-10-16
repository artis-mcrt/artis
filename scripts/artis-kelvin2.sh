#!/bin/bash -x
#SBATCH --ntasks=1024
#SBATCH --ntasks-per-node=128
#SBATCH --time=24:00:00
#SBATCH --partition=k2-math-physics
#SBATCH --mail-type=ALL
##SBATCH --mail-user=fmcneill07@qub.ac.uk

module load libs/gcc/14.1.0
# GSL is no longer used by artis, but is still loaded so that these scripts work with older versions too
module load gsl/2.8/gcc-14.1.0
module load compilers/gcc/14.1.0
module load mpi/openmpi/5.0.3/gcc-14.1.0
module load apps/python3/3.12.4/gcc-14.1.0

module list

cd "${SLURM_SUBMIT_DIR:?}" || exit 1

export MAKEFLAGS="--check-symlink-times --jobs=$(nproc)"
cd artis
# Build exspec together with sn3d, so that exspec reads the packet files with the same source version.
# The exspec job scripts can run exspec without a build.
make sn3d exspec || exit 1
cd ..

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"


hoursleft=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
source ./artis/scripts/corehours-before.sh
echo "$(date): before srun sn3d. hours left: $hoursleft"
time mpirun -- ./artis/sn3d -w $hoursleft
mpirun_status=$?
hoursleftafter=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): after srun sn3d finished. hours left: $hoursleftafter"
source ./artis/scripts/corehours-after.sh

# sn3d gives 0 also when it writes RESTART_NEEDED, so a non-zero status is a crash.
if [ $mpirun_status -ne 0 ]; then
    echo "$(date): mpirun sn3d gave status $mpirun_status, so this job submits nothing"
    exit $mpirun_status
fi

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    sbatch -J "${PWD##*/}" ${EMAIL:+--mail-user="$EMAIL"} -- artis/scripts/artis-kelvin2.sh
else
    # post-processing can remove restart files, so only queue it when no continuation job was submitted
    if ls packets/packets00_0000.out* > /dev/null 2>&1; then
        sbatch -J "exspec_${PWD##*/}" ${EMAIL:+--mail-user="$EMAIL"} -- artis/scripts/exspec-zip-kelvin2.sh
    fi
fi
