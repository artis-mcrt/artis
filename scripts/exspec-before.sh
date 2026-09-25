#!/usr/bin/env bash

# exspec-after.sh moves the packet files into packets/, so move them back before exspec.
if [[ ! -f packets00_0000.out && ! -f packets00_0000.out.zst && (-f packets/packets00_0000.out || -f packets/packets00_0000.out.zst) ]]; then
  mv packets/packets*.out* .
fi

# A build with libzstd reads the .zst files directly, so the script decompresses only the .xz
# files. A build without libzstd stops with a message that names the .zst file. Decompress the
# .zst files by hand for such a build, e.g. with: find . -maxdepth 1 -name '*.zst' -exec zstd -d --rm {} \;

find . -maxdepth 1 -name 'packets**.out.xz' -exec xz -d -v -T0 {} \;
find . -maxdepth 1 -name '*.txt.xz' -exec xz -d -v -T0 {} \;

if [[ -d data/ ]]; then
    find data/ -name '*.txt.xz' -exec xz -d -v -T0 {} \;
fi

if [[ -d artis/data/ ]]; then
    find artis/data/ -name '*.txt.xz' -exec xz -d -v -T0 {} \;
fi
