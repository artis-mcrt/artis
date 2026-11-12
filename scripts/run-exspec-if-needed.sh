#!/usr/bin/env bash

# The record identifies each spectrum, packet, and input file after a complete write.
if python3 ./artis/scripts/check-spectra-complete.py; then
  echo 'The spectrum files are complete, so exspec does not run'
else
  # the job script sources this file, so exit stops the job script before exspec-after.sh changes the model folder
  if ! ./artis/exspec; then
    echo "exspec failed, so the job stops before exspec-after.sh"
    exit 1
  fi
fi
