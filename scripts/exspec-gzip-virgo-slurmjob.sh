#!/bin/bash -x

export APPTAINER_CONTAINER="/cvmfs/vae.gsi.de/vae24/containers/vae24-user_container_20240418T1037.sif"
export APPTAINER_NAME="vae24-user_container"
export APPTAINER_SHARENS=true
export APPTAINER_CONFIGDIR=/tmp/$USER

eval `spack load --sh gsl%gcc arch=linux-debian11-x86_64`
eval `spack load --sh python arch=linux-debian11-x86_64`

cd $SLURM_SUBMIT_DIR

if [ ! -f emission.out.zst ]; then
  source ./artis/scripts/exspec-before.sh
  ./exspec
else
  echo 'Not running exspec because emission.out.zst was found'
fi

source ./artis/scripts/exspec-after.sh
