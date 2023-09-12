#!/usr/bin/env zsh

set -x

runfolder=classicmode_3d_testrun

mkdir -p $runfolder

rsync -av classicmode_3d_inputfiles/ $runfolder/

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory $runfolder/

cp ../data/* $runfolder/

cp ../artisoptions_classic.h $runfolder/artisoptions.h

cd $runfolder

xz -dv -T0 *.xz

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 15000;/g' artisoptions.h

sed -i'' -e 's/constexpr int GRID_TYPE.*/constexpr int GRID_TYPE = GRID_CARTESIAN3D;/g' artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' artisoptions.h

cd -

set +x
