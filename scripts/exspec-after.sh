if [ -f spec.out ]; then
  xz -v absorption.out emission*.out || true
  xz -v phixsdata_v2.txt transitiondata.txt ratecoeff.dat linestat.out || true
  mkdir packets || true
  mv packets*.out* packets/

  # gzip -v --best packets/packets*.out || true
  xz -v -T0 packets/packets*.out || true

  find . -name '*.out' -size +1M -exec xz -v -T0 {} \;
  find . -name 'output*.txt' -size +1M -exec xz -v -T0 {} \;
fi