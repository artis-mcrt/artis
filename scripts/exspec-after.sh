#!/usr/bin/env bash

if [[ -f emission.out || -f emission.out.xz ]]; then

  # join 3D direction files, if they exist
  ./artis/scripts/mergeangleres.py

  xz -v absorption.out emission*.out || true
  xz -v phixsdata*.txt transitiondata.txt ratecoeff_v2.dat linestat.out || true
  mkdir packets || true
  mv packets*.out* packets/

  # xz -v -T0 packets/packets*.out || true
  find packets/ -name 'packets*.out' -size +1M -exec xz -v -T0 {} \;

  find . -name '*.out' -size +1M -exec xz -v -T0 {} \;
  find . -mindepth 2 -name 'output_*.txt' ! -name "output_0-0.txt" -size +1M -exec xz -v -T0 {} \;

  # 3D kilonova model.txt and abundances.txt can be huge
  find . -maxdepth 1 -name '*.txt' -size +10M -exec xz -v -T0 {} \;
fi