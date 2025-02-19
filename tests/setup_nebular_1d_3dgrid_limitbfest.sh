#!/usr/bin/env zsh

set -x

runfolder=nebular_1d_3dgrid_limitbfest_testrun

rsync -av nebular_1d_3dgrid_inputfiles/ nebular_1d_3dgrid_limitbfest_testrun/

rsync --ignore-times -av nebular_1d_3dgrid_limitbfest_inputfiles/ nebular_1d_3dgrid_limitbfest_testrun/

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory nebular_1d_3dgrid_limitbfest_testrun/

cp ../data/* nebular_1d_3dgrid_limitbfest_testrun/

cp ../artisoptions_nltenebular.h nebular_1d_3dgrid_limitbfest_testrun/artisoptions.h

cd nebular_1d_3dgrid_limitbfest_testrun

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 1000000;/g' artisoptions.h

sed -i'' -e 's/constexpr auto GRID_TYPE.*/constexpr auto GRID_TYPE = GridType::CARTESIAN3D;/g' artisoptions.h

sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_X.*/constexpr int CUBOID_NCOORDGRID_X = 50;/g' artisoptions.h
sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_Y.*/constexpr int CUBOID_NCOORDGRID_Y = 50;/g' artisoptions.h
sed -i'' -e 's/constexpr int CUBOID_NCOORDGRID_Z.*/constexpr int CUBOID_NCOORDGRID_Z = 50;/g' artisoptions.h

sed -i'' -e 's/constexpr int TABLESIZE.*/constexpr int TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 2000.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 10000.;/g' artisoptions.h

sed -i'' -e 's/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP.*/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 7;/g' artisoptions.h

sed -i'' -e 's/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP.*/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 7;/g' artisoptions.h

sed -i'' -e 's/constexpr bool SF_AUGER_CONTRIBUTION_ON.*/constexpr bool SF_AUGER_CONTRIBUTION_ON = false;/g' artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' artisoptions.h

sed -i'' -e 's/  \/\/ return LEVEL_IS_NLTE.*/return LEVEL_IS_NLTE(element_z, ionstage, level);/g' artisoptions.h

cd -

set +x
