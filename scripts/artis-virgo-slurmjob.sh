#!/bin/bash
## SLURM META DIRECTIVES HERE DON'T WORK UNDER CENTOS VIRTUAL APPLICATION ENVIRONMENT
## So they are located in artis-virgo-submit.sh as cmd-line parameters to sbatch

export APPTAINER_CONTAINER="/cvmfs/vae.gsi.de/vae26/slurm-25-11/container/user_container-develop.sif"
export APPTAINER_NAME="vae26-user_container"
export APPTAINER_SHARENS=true
export APPTAINER_CONFIGDIR=/tmp/$USER

eval `spack load --first --sh openmpi%gcc`
# GSL is no longer used by artis, but is still loaded so that these scripts work with older versions too
eval `spack load --first --sh gsl%gcc`
eval `spack load --first --sh gcc%gcc`

export LD_LIBRARY_PATH=$(gsl-config --prefix)/lib/:$LD_LIBRARY_PATH
export MAKEFLAGS="--check-symlink-times --jobs=$(nproc --all)"
export OMPI_CXX=g++

cd $SLURM_SUBMIT_DIR

cd artis
make sn3d
cd ..

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

# decompress any zipped input files
source ./artis/scripts/exspec-before.sh

hoursleft=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "ntasks: ${SLURM_NTASKS}"
echo "cpus-per-task: ${SLURM_CPUS_PER_TASK:-1}"
echo "$(date): before srun sn3d. hours left: $hoursleft"

time srun -- ./artis/sn3d -w $hoursleft -o ${SLURM_JOB_ID}.slurm > out.txt

hoursleftafter=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): after srun sn3d finished. hours left: $hoursleftafter"
hourselapsed=$(python3 -c "print($hoursleft - $hoursleftafter)")
echo "hours of runtime: $hourselapsed"
cpuhrs=$(python3 -c "print($SLURM_NTASKS * ${SLURM_CPUS_PER_TASK:-1} * $hourselapsed)")
echo "ntasks: $SLURM_NTASKS -> cpus-per-task: ${SLURM_CPUS_PER_TASK:-1} -> CPU core hrs: $cpuhrs"

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    source ./artis/scripts/artis-virgo-submit.sh
else
    if [ -f packets00_0000.out ]; then
        source ./artis/scripts/exspec-zip-virgo-submit.sh
    fi
fi
