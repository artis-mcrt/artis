#!/bin/bash -x

spack load gsl target=$(spack arch -t)

cd $SLURM_SUBMIT_DIR

if [ ! -f spec.out ]; then
  ./artis/scripts/exspec-before.sh
  ./exspec
else
  echo 'Not running exspec because spec.out was found'
fi

./artis/scripts/exspec-after.sh
