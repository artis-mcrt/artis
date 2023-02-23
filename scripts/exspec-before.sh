#!/usr/bin/env bash

# files outside the artis folder, so move them back
if [ ! -f packets00_0000.out* ]; then
  mv packets/packets*.out* .
fi

# gzip -d packets*.out.gz || true
# xz -d -T 0 packets*.out.xz || true

#gzip -d -v phixsdata_v2.txt.gz transitiondata.txt.gz ratecoeff_v2.dat.gz || true
#xz -d -v -T 0 phixsdata_v2.txt.xz transitiondata.txt.xz ratecoeff_v2.dat.xz || true

find . -maxdepth 1 -name '*.gz' -size +10M -exec gzip -d -v {} \;
find . -maxdepth 1 -name '*.xz' -size +10M -exec xz -d -v -T0 {} \;
