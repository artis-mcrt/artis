if [ -f spec.out ]; then
  xz -v absorption.out emission*.out || true
  xz -v phixsdata_v2.txt transitiondata.txt ratecoeff.dat linestat.out || true
  mkdir packets || true
  mv packets*.out* packets/
  # xz -v -2 -T 0 packets/packets*.out || true
  gzip -v --best packets/packets*.out || true
fi