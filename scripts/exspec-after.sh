#!/usr/bin/env bash

# only compress the files if we successfully ran exspec
if [ -f emission.out* ]; then

  # zstd does decent compression at high speeds
  cmdcompress="zstd -T0 -16 -v --rm -f"

  # fall back to gzip
  # cmdcompress="gzip -v -f"

  # join 3D direction files, if they exist
  python3 ./artis/scripts/mergeangleres.py

  mkdir -p packets
  mv packets*.out* packets/ || true

  # 3D kilonova model.txt and abundances.txt can be huge, so compress txt files
  # do maxdepth 1 first in case job gets killed during run folder compression
  find . -maxdepth 1 -name '*.txt' ! -name "output_0-0.txt" -size +2M -exec $cmdcompress {} \;
  find . -maxdepth 1 -name '*.out' ! -name "slurm-*.out" -size +1M -exec $cmdcompress {} \;
  find . -maxdepth 1 -name 'rateceoff.dat' ! -name "slurm-*.out" -size +1M -exec $cmdcompress {} \;

  find packets/ -name 'packets*.out' -size +1M -exec $cmdcompress {} \;

  find . -maxdepth 2 -name '*.txt' ! -name "output_0-0.txt" -size +2M -exec $cmdcompress {} \;
  find . -maxdepth 2 -name '*.out' ! -name "slurm-*.out" -size +1M -exec $cmdcompress {} \;

  ./artis/scripts/tar_rm_logs.sh

fi

