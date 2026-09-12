#!/usr/bin/env bash

# sn3d writes emission.out at its last requested timestep
if [[ -f emission.out || -f emission.out.zst ]]; then
  echo 'Not running exspec because emission.out[.zst] was found'
else
  source ./artis/scripts/exspec-before.sh
  ./artis/exspec
fi
