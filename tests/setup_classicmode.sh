#!/usr/bin/env zsh

set -x

rsync -av classicmode_inputfiles/ classicmode_testrun/

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory classicmode_testrun/

cp ../data/* classicmode_testrun/

cp ../artisoptions_classic.h classicmode_testrun/artisoptions.h

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 15000;/g' classicmode_testrun/artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' classicmode_testrun/artisoptions.h

sed -i'' -e 's/constexpr bool VPKT_ON.*/constexpr bool VPKT_ON = true;,g' classicmode_testrun/artisoptions.h

set +x
