#!/usr/bin/env bash

paths="*.tmp *.out output_*.txt exspec.txt *.out.xz *.out.gz"

# if [[ "$1" == "-d" ]]; then
#   echo 1
# else
#   echo 1>&2 "ADD -d TO CONFIRM DELETION OF THE FOLLOWING FILES"
#   ls $paths 2>/dev/null
# fi

if [ 0 -lt $(ls $paths 2>/dev/null | wc -w) ]; then
  echo "The following ARTIS run files will be deleted:"
  ls $paths 2>/dev/null

  read -p "Are you sure you want to delete these ARTIS run files? " -n 1 -r
  echo    # (optional) move to a new line
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Deleting:"
    rm -v $paths 2>/dev/null
  fi
else
  echo "No ARTIS run files to delete"
fi

