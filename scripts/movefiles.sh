#!/usr/bin/env bash

# this script will move the ARTIS output files from the current run into a subfolder
# given by the command line argument
# e.g. ./artis/scripts/movefiles.sh 1234589.slurm

if [ $# -ne 1 ]; then
  echo 1>&2 "Usage: $0 JOB_DIRECTORY"
  exit 3
fi

# mkdir -p $1

# check that the argument is a valid directory or exit
if [ ! -d "$1" ]; then
  echo 1>&2 "JOB_DIRECTORY '$1' is not an existing directory"
  echo 1>&2 "to create, run mkdir '$1'"
  exit 3
fi

mkdir -p $1

#if you want to be able to resume from the end of this run
#cp gridsave.dat packets*.tmp $1

#mv -n will trigger an error if the destination exists
mv -n estimators* $1
# move the files below if the exist, but don't cause an error if they don't exist!
find . -maxdepth 1 -name 'macroatom*' -exec mv -n {} $1 \;
find . -maxdepth 1 -name 'nonthermal*' -exec mv -n {} $1 \;
find . -maxdepth 1 -name 'radfield*' -exec mv -n {} $1 \;
find . -maxdepth 1 -name 'nlte*' -exec mv -n {} $1 \;
mv -n output_*.txt $1
cp $1/output_0-0.txt .
