#!/usr/bin/env bash

# A job script sources this file, and a sourced file sees the arguments of the job script.
if [ "${BASH_SOURCE:-$0}" = "$0" ] && [ $# -gt 0 ]; then
  echo "Usage: artis/scripts/run-exspec-if-needed.sh with no argument, in the model folder. The script runs exspec if emission.out is absent." >&2
  exit 1
fi

# sn3d writes emission.out at its last requested timestep
if [[ -f emission.out || -f emission.out.zst ]]; then
  echo 'emission.out[.zst] exists, so exspec does not run'
else
  # the job script sources this file, so exit stops the job script before exspec-after.sh changes the model folder
  if ! ./artis/exspec; then
    echo "exspec failed, so the job stops before exspec-after.sh"
    exit 1
  fi
fi
