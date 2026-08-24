#!/usr/bin/env zsh

set -x

runfolder=classicmode_1d_3dgrid_testrun

if [ ! -f atomicdata_classic.tar.xz ]; then curl -O -L https://github.com/artis-mcrt/artis/releases/download/v2026.5.15/atomicdata_classic.tar.xz; fi

mkdir -p $runfolder

cd $runfolder

tar -xf ../atomicdata_classic.tar.xz --directory ./

rsync -av ../classicmode_1d_3dgrid_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_classic.h artisoptions.h

sed -i.bak -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 15000;/g' artisoptions.h

sed -i.bak -e 's/constexpr std::optional<GridType> GRID_TYPE_OVERRIDE.*/constexpr std::optional<GridType> GRID_TYPE_OVERRIDE = GridType::CARTESIAN3D;/g' artisoptions.h

sed -i.bak -e 's/constexpr int CUBOID_NCOORDGRID_X.*/constexpr int CUBOID_NCOORDGRID_X = 100;/g' artisoptions.h
sed -i.bak -e 's/constexpr int CUBOID_NCOORDGRID_Y.*/constexpr int CUBOID_NCOORDGRID_Y = 100;/g' artisoptions.h
sed -i.bak -e 's/constexpr int CUBOID_NCOORDGRID_Z.*/constexpr int CUBOID_NCOORDGRID_Z = 100;/g' artisoptions.h

sed -i.bak -e 's/constexpr bool VPKT_ON.*/constexpr bool VPKT_ON = true;/g' artisoptions.h
sed -i.bak -e 's/constexpr bool VPKT_WRITE_CONTRIBS.*/constexpr bool VPKT_WRITE_CONTRIBS = true;/g' artisoptions.h

rm -f artisoptions.h.bak

cd -

set +x
