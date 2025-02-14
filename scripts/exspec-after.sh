#!/usr/bin/env bash

# only compress the files if we successfully ran exspec
if [[ -f emission.out || -f emission.out.zst ]]; then

  # zstd does decent compression at high speeds
  cmdcompress="zstd -T0 -16 -v --rm -f"

  # fall back to gzip
  # cmdcompress="gzip -v -f"

  # join 3D direction files, if they exist
  python3 ./artis/scripts/mergeangleres.py

  mkdir -p packets
  mv packets*.out* packets/ || true

  mkdir -p vpackets
  mv vpackets*.out* vpackets/ || true

  mkdir -p vspecpol
  mv vspecpol*.out* vspecpol/ || true

  # remove empty directories
  find . -maxdepth 1 -type d -empty -delete

  # 3D kilonova model.txt and abundances.txt can be huge, so compress txt files
  # do maxdepth 1 first in case job gets killed during run folder compression
  find . -maxdepth 1 -name '*.txt' ! -name "output_0-0.txt" -size +2M -exec $cmdcompress {} \;
  find . -maxdepth 1 -name '*.out' ! -name "slurm-*.out" -size +1M -exec $cmdcompress {} \;
  find . -maxdepth 1 -name 'ratecoeff.dat' -size +1M -exec $cmdcompress {} \;

  find packets/ -name 'packets*.out' -size +1M -exec $cmdcompress {} \;

  find . -maxdepth 2 -name '*.txt' ! -name "output_0-0.txt" -size +2M -exec $cmdcompress {} \;
  find . -maxdepth 2 -name '*.out' ! -name "slurm-*.out" -size +1M -exec $cmdcompress {} \;

  ./artis/scripts/tar_rm_logs.sh

  mkdir -p speclc_angle_res
  mv *_res_*.out* speclc_angle_res/ || true

  curl -LsSf https://astral.sh/uv/install.sh | sh
  uv tool install artistools

  # convert packets to parquet for fast reading
  uvx artistools lc --frompackets || true

  if [ -f vpkt.txt ]; then
    # convert virtual packets to parquet
    uvx artistools lc --frompackets -plotvspecpol 0 || true
  fi

  # convert estimators to parquet. On JUWELS, you might need to limit the number of processes to 16 in artistools/artistools/configuration.py
  uvx --from artistools -- python3 -c 'import artistools as at; at.estimators.scan_estimators()' || true

fi
