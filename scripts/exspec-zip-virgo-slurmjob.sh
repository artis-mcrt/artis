#!/bin/bash

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
make exspec
cd ..

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"


source ./artis/scripts/run-exspec-if-needed.sh

source ./artis/scripts/exspec-after.sh
