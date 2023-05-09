#!/usr/bin/env zsh

set -x

rsync -av classicmode_3d_inputfiles/ classicmode_3d_testrun/

xz -d -T0 -v classicmode_3d_testrun/*.xz

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory classicmode_3d_testrun/

cp ../data/* classicmode_3d_testrun/

cp ../artisoptions_classic.h classicmode_3d_testrun/artisoptions.h

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 15000;/g' classicmode_3d_testrun/artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' classicmode_3d_testrun/artisoptions.h


set +x
