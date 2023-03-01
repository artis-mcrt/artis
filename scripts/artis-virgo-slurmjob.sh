#!/bin/bash -x
## SLURM META DIRECTIVES HERE DON'T WORK UNDER CENTOS VIRTUAL APPLICATION ENVIRONMENT
## So they are located in artis-virgo-submit.sh as cmd-line parameters to sbatch

cd $SLURM_SUBMIT_DIR

spack load gsl target=$(spack arch -t)

# decompress any zipped input files
source ./artis/scripts/exspec-before.sh

hoursleft=$(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})
echo "$(date): before srun sn3d. hours left: $hoursleft"
time srun -- ./sn3d -w $hoursleft > out.txt
echo "$(date): after srun sn3d finished. hours left: $(python3 ./artis/scripts/slurmjobhoursleft.py ${SLURM_JOB_ID})"
echo "seconds elapsed: $SECONDS"

mkdir ${SLURM_JOB_ID}.slurm
source ./artis/scripts/movefiles.sh ${SLURM_JOB_ID}.slurm

if grep -q "RESTART_NEEDED" "output_0-0.txt"
then
    # check if there's a submit script in the submission directory, otherwise use the repository version
    if [ -f ./artis-virgo-submit.sh ]
    then
        source ./artis-virgo-submit.sh
    else
        source ./artis/scripts/artis-virgo-submit.sh
    fi
fi

if [ -f packets00_0000.out ]; then
    source ./artis/scripts/exspec-gzip-virgo-submit.sh
fi
