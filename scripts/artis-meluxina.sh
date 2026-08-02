#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --ntasks=1024
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --account=p200999
#SBATCH --mail-type=ALL
##SBATCH --mail-user=luke.shingles@gmail.com

module load env/release/2025.1 gompi/2025a zstd git Python

module list

cd $SLURM_SUBMIT_DIR

export MAKEFLAGS="--check-symlink-times --jobs=$(nproc --all)"
export OMPI_CXX=g++
cd artis
make
cd ..

mpicxx --version
echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

# decompress any zipped input files
source ./artis/scripts/exspec-before.sh
hoursleft=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): before srun sn3d. hours left: $hoursleft"
time srun --hint=nomultithread -- ./artis/sn3d -w $hoursleft > out.txt
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
    sbatch --job-name="$SLURM_JOB_NAME" ./artis/scripts/artis-meluxina.sh
fi

if [ -f packets00_0000.out ]; then
    sbatch --job-name="exspec_$SLURM_JOB_NAME" ./artis/scripts/exspec-zip-meluxina.sh
fi
