#!/usr/bin/env zsh

set -x

runfolder=classicmode_1d_3dgrid_testrun

mkdir -p $runfolder

if [ ! -f atomicdata_classic.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_classic.tar.xz; fi

tar -xf atomicdata_classic.tar.xz --directory $runfolder/

rsync -av classicmode_1d_3dgrid_inputfiles/ $runfolder/

cp ../data/* $runfolder/

cp ../artisoptions_classic.h $runfolder/artisoptions.h

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 15000;/g' $runfolder/artisoptions.h

sed -i'' -e 's/constexpr int GRID_TYPE.*/constexpr int GRID_TYPE = GRID_CARTESIAN3D;/g' artisoptions.h

sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_X.*/constexpr int CUBOID_NCOORDGRID_X = 100;/g' classicmode_3d_testrun/artisoptions.h
sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_Y.*/constexpr int CUBOID_NCOORDGRID_Y = 100;/g' classicmode_3d_testrun/artisoptions.h
sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_Z.*/constexpr int CUBOID_NCOORDGRID_Z = 100;/g' classicmode_3d_testrun/artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' $runfolder/artisoptions.h

sed -i'' -e 's/constexpr bool VPKT_ON.*/constexpr bool VPKT_ON = true;/g' $runfolder/artisoptions.h

set +x
