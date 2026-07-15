#!/usr/bin/env bash

if [[ -f emission.out || -f emission.out.zst ]] && [[ -f exspec.txt ]]; then
  echo 'Not running exspec because emission.out[.zst] and exspec.txt were found'
else
  source ./artis/scripts/exspec-before.sh
  ./artis/exspec
fi
