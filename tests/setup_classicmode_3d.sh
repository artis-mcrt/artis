#!/usr/bin/env zsh

set -x

runfolder=classicmode_3d_testrun

mkdir -p $runfolder

if [ ! -f atomicdata_classic.tar.xz ]; then curl --insecure -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_classic.tar.xz; fi

tar -xf atomicdata_classic.tar.xz --directory $runfolder/

rsync -av classicmode_3d_inputfiles/ $runfolder/

ln -s ../../data/ $runfolder

cp ../artisoptions_classic.h $runfolder/artisoptions.h

cd $runfolder

xz -dv -T0 *.xz

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 15000;/g' artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_EMISSIONABSORPTION_SPEC_AT_END.*/constexpr bool WRITE_EMISSIONABSORPTION_SPEC_AT_END = true;/g' artisoptions.h

cd -

set +x
