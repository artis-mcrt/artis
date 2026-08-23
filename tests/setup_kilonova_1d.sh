#!/usr/bin/env zsh

set -x

runfolder=kilonova_1d_testrun

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O -L https://github.com/artis-mcrt/artis/releases/download/v2026.5.15/atomicdata_feconi.tar.xz; fi

mkdir -p $runfolder

cd $runfolder

tar -xf ../atomicdata_feconi.tar.xz --directory ./

rsync -av ../kilonova_1d_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_kilonova_lte.h artisoptions.h

xz -f -d -v -T0 *.xz

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 80000;/g' artisoptions.h

sed -i'' -e 's/constexpr int TABLESIZE.*/constexpr int TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 1000.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 20000.;/g' artisoptions.h


cd -

set +x
