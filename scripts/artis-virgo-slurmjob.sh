#!/bin/bash -x
## SLURM META DIRECTIVES HERE DON'T WORK UNDER CENTOS VIRTUAL APPLICATION ENVIRONMENT
## So they are located in artis-virgo-submit.sh as cmd-line parameters to sbatch

cd $SLURM_SUBMIT_DIR

spack load gsl target=$(spack arch -t)

echo "before srun sn3d"
srun -- ./sn3d -w 24 > out.txt
echo "after srun"

mkdir ${SLURM_JOBID}.slurm
./artis/scripts/movefiles.sh ${SLURM_JOBID}.slurm

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    # check if there's a submit script in the submission directory, otherwise use the repository version
    if [ -f ./artis-virgo-submit.sh ]
    then
        ./artis-virgo-submit.sh
    else
        ./artis/scripts/artis-virgo-submit.sh
    fi
fi

if [ -f packets00_0000.out ]; then
    ./artis/scripts/exspec-gzip-virgo-submit.sh
fi
