#!/usr/bin/env zsh

set -x

runfolder=classicmode_3d_testrun

mkdir -p $runfolder

if [ ! -f atomicdata_classic.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_classic.tar.xz; fi

tar -xf atomicdata_classic.tar.xz --directory $runfolder/

rsync -av classicmode_3d_inputfiles/ $runfolder/

cp ../data/* $runfolder/

cp ../artisoptions_classic.h $runfolder/artisoptions.h

cd $runfolder

xz -dv -T0 *.xz

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 15000;/g' artisoptions.h

sed -i'' -e 's/constexpr auto GRID_TYPE.*/constexpr auto GRID_TYPE = GridType::CARTESIAN3D;/g' artisoptions.h

sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_X.*/constexpr int CUBOID_NCOORDGRID_X = 10;/g' artisoptions.h
sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_Y.*/constexpr int CUBOID_NCOORDGRID_Y = 10;/g' artisoptions.h
sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_Z.*/constexpr int CUBOID_NCOORDGRID_Z = 10;/g' artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' artisoptions.h

cd -

set +x
