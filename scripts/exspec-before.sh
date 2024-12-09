#!/usr/bin/env bash

# packet output files outside the artis folder, so move them back to run exspec
if [[ ! -f packets00_0000.out && ! -f packets00_0000.out.zst && -f packets/packets00_0000.out.zst ]]; then
  mv packets/packets*.out* .
fi

find . -maxdepth 1 -name 'packets**.out.zst' -exec zstd -d -v -T0 --rm {} \;
find . -maxdepth 1 -name 'packets**.out.gz' -exec gzip -d -v {} \;
find . -maxdepth 1 -name 'packets**.out.xz' -exec xz -d -v -T0 {} \;

#unzip, e.g. phixsdata_v2.txt.xz transitiondata.txt.xz ratecoeff.dat.xz

find . -maxdepth 1 -name '*.dat.zst' -exec zstd -d -v -T0 --rm {} \;
find . -maxdepth 1 -name '*.dat.xz' -exec xz -d -v -T0 {} \;
find . -maxdepth 1 -name '*.dat.gz' -exec gzip -d -v {} \;

find . -maxdepth 1 -name '*.txt.zst' -exec zstd -d -v -T0 --rm {} \;
find . -maxdepth 1 -name '*.txt.xz' -exec xz -d -v -T0 {} \;
find . -maxdepth 1 -name '*.txt.gz' -exec gzip -d -v {} \;
