# files not in the artis folder, so move them
if [ ! -f packets00_0000.out* ]; then
  mv packets/packets*.out* .
fi

gzip -d packets*.out.gz || true
xz -d -T 0 packets*.out.xz || true
gzip -d -v phixsdata_v2.txt.gz transitiondata.txt.gz ratecoeff.dat.gz || true
xz -d -v -T 0 phixsdata_v2.txt.xz transitiondata.txt.xz ratecoeff.dat.xz || true
