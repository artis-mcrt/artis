#!/bin/bash

export APPTAINER_CONTAINER="/cvmfs/vae.gsi.de/vae25/containers/user_container-production.sif"
export APPTAINER_NAME="vae25-user_container"
export APPTAINER_SHARENS=true
export APPTAINER_CONFIGDIR=/tmp/$USER

eval `spack load --first --sh openmpi%gcc`
eval `spack load --first --sh gsl%gcc target=x86_64`
eval `spack load --first --sh gcc`
eval `spack load --first --sh zstd`

export LD_LIBRARY_PATH=$(gsl-config --prefix)/lib/:$LD_LIBRARY_PATH
export MAKEFLAGS="--check-symlink-times --jobs=$(nproc --all)"
export OMPI_CXX=g++

cd $SLURM_SUBMIT_DIR

cd artis
make exspec
cd ..

echo "CPU type: $(c++ -march=native -Q --help=target | grep -- '-march=  ' | cut -f3)"

if [[ -f emission.out || -f emission.out.zst ]] && [[ -f exspec.txt ]]; then
  echo 'Not running exspec because emission.out[.zst] and exspec.txt were found'
else
  source ./artis/scripts/exspec-before.sh
  ./artis/exspec
fi

source ./artis/scripts/exspec-after.sh
