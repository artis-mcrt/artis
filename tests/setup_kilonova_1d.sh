#!/usr/bin/env zsh

set -x

runfolder=kilonova_1d_testrun

mkdir -p $runfolder

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O -L https://github.com/artis-mcrt/artis/releases/download/v2026.5.15/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory $runfolder/

# same input files as the other test run
rsync -av kilonova_1d_inputfiles/ $runfolder/

ln -s ../../data/ $runfolder

cp ../artisoptions_kilonova_lte.h $runfolder/artisoptions.h

cd $runfolder

xz -f -d -v -T0 *.xz

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 80000;/g' artisoptions.h

sed -i'' -e 's/constexpr int TABLESIZE.*/constexpr int TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 1000.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 20000.;/g' artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_EMISSIONABSORPTION_SPEC_AT_END.*/constexpr bool WRITE_EMISSIONABSORPTION_SPEC_AT_END = true;/g' artisoptions.h

cd -

set +x
