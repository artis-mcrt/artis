#!/usr/bin/env zsh

set -x

mkdir -p kilonova_testrun

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory kilonova_testrun/

rsync -av kilonova_inputfiles/ kilonova_testrun/

cp ../data/* kilonova_testrun/

cp ../artisoptions_kilonova_lte.h kilonova_testrun/artisoptions.h

cd kilonova_testrun

xz -dvk -T0 *.xz

sed -i'' -e 's/#define MPKTS.*/#define MPKTS 40000/g' artisoptions.h

sed -i'' -e 's/constexpr int TABLESIZE.*/constexpr int TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 1000.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 20000.;/g' artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' artisoptions.h

cd -

set +x
