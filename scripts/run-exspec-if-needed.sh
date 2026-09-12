#!/usr/bin/env bash

# sn3d writes emission.out at its last requested timestep
if [[ -f emission.out || -f emission.out.zst ]]; then
  echo 'Not running exspec because emission.out[.zst] was found'
else
  source ./artis/scripts/exspec-before.sh
  # this file is sourced, so exit stops the job script before exspec-after.sh changes the run folder
  if ! ./artis/exspec; then
    echo "exspec failed, so the job stops before exspec-after.sh"
    exit 1
  fi
fi
