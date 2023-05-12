#!/usr/bin/env bash

# don't compress the files if we didn't successfully run exspec
if [[ -f emission.out || -f emission.out.xz || -f emission.out.lz4 || -f emission.out.zst ]]; then

  if command -v zstd > /dev/null; then
    # zstd does decent compression at high speeds
    cmdcompress="zstd -T0 -16 -v --rm -f"
  elif command -v lz4 > /dev/null; then
    # lz4 is extremely fast, but low compression ratios
    cmdcompress="lz4 -v --best --rm -f"
  else
    # fall back to gzip
    cmdcompress="gzip -v -f"
  fi

  # join 3D direction files, if they exist
  ./artis/scripts/mergeangleres.py

  mkdir -p packets
  mv packets*.out* packets/ || true

  # 3D kilonova model.txt and abundances.txt can be huge, so compress txt files
  # do maxdepth 1 first in case job gets killed during run folder compression
  find . -maxdepth 1 -name '*.txt' ! -name "output_0-0.txt" -size +2M -exec $cmdcompress {} \;
  find . -maxdepth 1 -name '*.out' -size +1M -exec $cmdcompress {} \;

  find packets/ -name 'packets*.out' -size +1M -exec $cmdcompress {} \;

  find . -maxdepth 2 -name '*.txt' ! -name "output_0-0.txt" -size +2M -exec $cmdcompress {} \;
  find . -maxdepth 2 -name '*.out' -size +1M -exec $cmdcompress {} \;

fi

