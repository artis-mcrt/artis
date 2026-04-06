#!/usr/bin/env bash

# packet output files outside the artis folder, so move them back to run exspec
if [[ ! -f packets00_0000.out && ! -f packets00_0000.out.zst && (-f packets/packets00_0000.out || -f packets/packets00_0000.out.zst) ]]; then
  mv packets/packets*.out* .
fi

find . -maxdepth 1 -name 'packets**.out.zst' -exec zstd -d -v -T0 --rm {} \;
find . -maxdepth 1 -name 'packets**.out.xz' -exec xz -d -v -T0 {} \;

#unzip, e.g. phixsdata_v2.txt.xz transitiondata.txt.xz ratecoeff.dat.xz

find . -maxdepth 1 -name '*.dat.zst' -exec zstd -d -v -T0 --rm {} \;
find . -maxdepth 1 -name '*.dat.xz' -exec xz -d -v -T0 {} \;

find . -maxdepth 1 -name '*.txt.zst' -exec zstd -d -v -T0 --rm {} \;
find . -maxdepth 1 -name '*.txt.xz' -exec xz -d -v -T0 {} \;

# if data folder exists, unzip files there too
if [[ -d data/ ]]; then
    find data/ -name '*.txt.zst' -exec zstd -d -v -T0 --rm {} \;
    find data/ -name '*.txt.xz' -exec xz -d -v -T0 {} \;
fi

if [[ -d artis/data/ ]]; then
    find artis/data/ -name '*.txt.zst' -exec zstd -d -v -T0 --rm {} \;
    find artis/data/ -name '*.txt.xz' -exec xz -d -v -T0 {} \;
fi
