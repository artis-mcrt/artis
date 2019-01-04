if [ -f spec.out ]; then
  gzip -v absorption.out emission*.out || true
  gzip -v phixsdata_v2.txt transitiondata.txt ratecoeff.dat linestat.out || true
  gzip -v packets*.out || true
  mkdir packets || true
  mv packets*.out.gz packets/
fi
