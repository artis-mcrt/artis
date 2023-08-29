#!/usr/bin/env zsh

set -x

runfolder=classicmode_testrun

mkdir -p $runfolder

if [ ! -f atomicdata_classic.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_classic.tar.xz; fi

tar -xf atomicdata_classic.tar.xz --directory $runfolder/

rsync -av classicmode_inputfiles/ $runfolder/

cp ../data/* $runfolder/

cp ../artisoptions_classic.h $runfolder/artisoptions.h

cd $runfolder

xz -dvk -T0 *.xz

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 80000;/g' artisoptions.h

sed -i'' -e 's/constexpr int TABLESIZE.*/constexpr int TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 1000.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 20000.;/g' artisoptions.h

cd -

set +x
