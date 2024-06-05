#!/bin/bash -x

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
