#!/usr/bin/env zsh

set -x

runfolder=nebular_1d_3dgrid_testrun

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O -L https://github.com/artis-mcrt/artis/releases/download/v2026.5.15/atomicdata_feconi.tar.xz; fi

mkdir -p $runfolder

cd $runfolder

rsync -av ../nebular_1d_3dgrid_inputfiles/ ./

tar -xf ../atomicdata_feconi.tar.xz --directory .

ln -s ../../ artis

cp artis/artisoptions_nltenebular.h artisoptions.h

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 1000000;/g' artisoptions.h

sed -i'' -e 's/constexpr std::optional<GridType> GRID_TYPE_OVERRIDE.*/constexpr std::optional<GridType> GRID_TYPE_OVERRIDE = GridType::CARTESIAN3D;/g' artisoptions.h

sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_X.*/constexpr int CUBOID_NCOORDGRID_X = 50;/g' artisoptions.h
sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_Y.*/constexpr int CUBOID_NCOORDGRID_Y = 50;/g' artisoptions.h
sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_Z.*/constexpr int CUBOID_NCOORDGRID_Z = 50;/g' artisoptions.h

sed -i'' -e 's/constexpr int TABLESIZE.*/constexpr int TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 2000.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 10000.;/g' artisoptions.h

sed -i'' -e 's/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP.*/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 7;/g' artisoptions.h

sed -i'' -e 's/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP.*/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 7;/g' artisoptions.h


cd -

set +x
