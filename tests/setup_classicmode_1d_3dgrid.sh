#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

runfolder=classicmode_1d_3dgrid_testrun

getatomicdata atomicdata_classic.tar.xz

mkdir -p $runfolder

cd $runfolder

tar -xf ../atomicdata_classic.tar.xz --directory ./

rsync -av ../classicmode_1d_3dgrid_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_classic.h artisoptions.h

sedopt "constexpr std::int64_t NUM_PACKETS.*" "constexpr std::int64_t NUM_PACKETS = 60'000;"

sedopt 'constexpr std::optional<GridType> GRID_TYPE_OVERRIDE.*' 'constexpr std::optional<GridType> GRID_TYPE_OVERRIDE = GridType::CARTESIAN3D;'

sedopt 'constexpr int CUBOID_NCOORDGRID_X.*' 'constexpr int CUBOID_NCOORDGRID_X = 100;'
sedopt 'constexpr int CUBOID_NCOORDGRID_Y.*' 'constexpr int CUBOID_NCOORDGRID_Y = 100;'
sedopt 'constexpr int CUBOID_NCOORDGRID_Z.*' 'constexpr int CUBOID_NCOORDGRID_Z = 100;'

sedopt 'constexpr bool VPKT_ON.*' 'constexpr bool VPKT_ON = true;'
sedopt 'constexpr bool VPKT_WRITE_CONTRIBS.*' 'constexpr bool VPKT_WRITE_CONTRIBS = true;'

sedopt 'constexpr bool KEEP_ESCAPED_GAMMAS.*' 'constexpr bool KEEP_ESCAPED_GAMMAS = true;'

rm -f artisoptions.h.bak

cd -

set +x
