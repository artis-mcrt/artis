#!/bin/bash -x
##SLURM DIRECTIVES HERE DON'T WORK UNDER CENTOS VIRTUAL APPLICATION ENVIRONMENT

cd $SLURM_SUBMIT_DIR

spack load gsl target=$(spack arch -t)

echo "before srun sn3d"
srun -- ./sn3d -w 8 > out.txt
echo "after srun"

mkdir ${SLURM_JOBID}.slurm
./artis/scripts/movefiles.sh ${SLURM_JOBID}.slurm

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    ./artis-virgo-submit.sh
fi

if [ -f packets00_0000.out ]; then
    ./artis/scripts/exspec-gzip-virgo-submit.sh
fi
