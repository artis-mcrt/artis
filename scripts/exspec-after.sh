if [ -f spec.out ]; then
  xz -v absorption.out emission*.out || true
  xz -v phixsdata_v2.txt transitiondata.txt ratecoeff.dat linestat.out || true
  mkdir packets || true
  mv packets*.out* packets/
  xz -v -T 0 packets*.out packets/packets*.out || true
fi