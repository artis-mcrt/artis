if [ -f spec.out ]; then
  gzip -v absorption.out emission*.out || true
  gzip -v packets*.out || true
  gzip -v phixsdata.txt transitiondata.txt ratecoeff.dat linestate.out || true
  mkdir packets || true
  mv packets*.out.gz packets/
fi
