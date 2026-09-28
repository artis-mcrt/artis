#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

runfolder=nebular_1d_3dgrid_limitbfest_testrun

getatomicdata atomicdata_feconi.tar.xz

mkdir -p $runfolder

cd $runfolder

rsync -av ../nebular_1d_3dgrid_inputfiles/ ./

rsync --ignore-times -av ../nebular_1d_3dgrid_limitbfest_inputfiles/ ./

tar -xf ../atomicdata_feconi.tar.xz --directory .

ln -s ../../ artis

cp artis/artisoptions_nltenebular.h artisoptions.h

sedopt "constexpr std::int64_t NUM_PACKETS.*" "constexpr std::int64_t NUM_PACKETS = 4'000'000;"

sedopt 'constexpr std::optional<GridType> GRID_TYPE_OVERRIDE.*' 'constexpr std::optional<GridType> GRID_TYPE_OVERRIDE = GridType::CARTESIAN3D;'

sedopt 'constexpr int CUBOID_NCOORDGRID_X.*' 'constexpr int CUBOID_NCOORDGRID_X = 50;'
sedopt 'constexpr int CUBOID_NCOORDGRID_Y.*' 'constexpr int CUBOID_NCOORDGRID_Y = 50;'
sedopt 'constexpr int CUBOID_NCOORDGRID_Z.*' 'constexpr int CUBOID_NCOORDGRID_Z = 50;'

sedopt 'constexpr int RATECOEFF_TABLESIZE.*' 'constexpr int RATECOEFF_TABLESIZE = 20;'
sedopt 'constexpr double MINTEMP.*' 'constexpr double MINTEMP = 2000.;'
sedopt 'constexpr double MAXTEMP.*' 'constexpr double MAXTEMP = 10000.;'

sedopt 'constexpr int FIRST_NLTE_RADFIELD_TIMESTEP.*' 'constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 7;'

sedopt 'constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP.*' 'constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 7;'

sedopt 'constexpr bool SF_AUGER_CONTRIBUTION_ON.*' 'constexpr bool SF_AUGER_CONTRIBUTION_ON = false;'

sedopt 'constexpr bool LEVEL_HAS_BFEST.*' 'constexpr bool LEVEL_HAS_BFEST(int element_z, int ionstage, int level) { return level <= ION_NLEVELS_EXCITED_NLTE(element_z, ionstage); }'

rm -f artisoptions.h.bak

cd -

set +x
