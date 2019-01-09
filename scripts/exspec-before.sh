if [! -f packets00_0000.out* ]; then
  mv packets/packets*.out* .
fi

gunzip packets*.out.gz || true
gunzip -v phixsdata_v2.txt.gz transitiondata.txt.gz ratecoeff.dat.gz || true
