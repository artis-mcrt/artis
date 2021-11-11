#!/bin/bash

paths="*.tmp *.out output_*.txt exspec.txt"

# if [[ "$1" == "-d" ]]; then
#   echo 1
# else
#   echo 1>&2 "ADD -d TO CONFIRM DELETION OF THE FOLLOWING FILES"
#   ls $paths 2>/dev/null
# fi

echo "The following files will be deleted"
ls $paths 2>/dev/null

read -p "Are you sure? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]; then
  rm -v $paths
fi