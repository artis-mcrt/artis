#!/usr/bin/env bash

if [ -f emission.out ]; then

  # join 3D direction files, if they exist
  ./artis/scripts/mergeangleres.py

  xz -v absorption.out emission*.out || true
  xz -v phixsdata*.txt transitiondata.txt ratecoeff.dat linestat.out || true
  mkdir packets || true
  mv packets*.out* packets/

  # gzip -v --best packets/packets*.out || true
  xz -v -T0 packets/packets*.out || true

  find . -name '*.out' -size +1M -exec xz -v -T0 {} \;
  find . -mindepth 2 -name 'output_*.txt' ! -name "output_0-0.txt" -size +1M -exec xz -v -T0 {} \;
fi