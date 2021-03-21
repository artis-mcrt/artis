#!/usr/bin/env bash

# this script will move the ARTIS output files from the current run into a subfolder
# given by the command line argument

if [ $# -ne 1 ]; then
  echo 1>&2 "Usage: $0 JOB_DIRECTORY"
  exit 3
fi

#mkdir $1

# check that the argument is a valid directory or exit
[ -d "$1" ] || exit

#cp gridsave.dat $1
mv estimators* $1
# move the files below if the exist, but don't cause an error if they don't exist!
find . -maxdepth 1 -name 'macroatom*' -type f -print0 | xargs -0r mv -t $1
find . -maxdepth 1 -name 'nonthermal*' -type f -print0 | xargs -0r mv -t $1
find . -maxdepth 1 -name 'radfield*' -type f -print0 | xargs -0r mv -t $1
find . -maxdepth 1 -name 'nlte*' -type f -print0 | xargs -0r mv -t $1
mv output_*.txt $1
cp $1/output_0-0.txt .
