#!/usr/bin/env zsh

set -x

rsync -av kilonova_inputfiles/ kilonova_testrun/

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory kilonova_testrun/

cp ../data/* kilonova_testrun/

cp ../artisoptions_kilonova_lte.h kilonova_testrun/artisoptions.h

cd kilonova_testrun

xz -dvk -T0 *.xz

sed -i'' -e 's/#define MPKTS.*/#define MPKTS 5000/g' artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' artisoptions.h

cd -

set +x
