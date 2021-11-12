#!/bin/bash -x

spack load gsl target=$(spack arch -t)

cd $SLURM_SUBMIT_DIR

if [ ! -f emission.out.xz ]; then
  source ./artis/scripts/exspec-before.sh
  ./exspec
else
  echo 'Not running exspec because emission.out.xz was found'
fi

source ./artis/scripts/exspec-after.sh
