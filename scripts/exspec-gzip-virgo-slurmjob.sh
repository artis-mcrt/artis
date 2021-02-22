#!/bin/bash -x

spack load gsl target=$(spack arch -t)

cd $SLURM_SUBMIT_DIR

./artis/scripts/exspec-before.sh

./exspec

./artis/scripts/exspec-after.sh
