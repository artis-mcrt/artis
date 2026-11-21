#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

modelpath=kilonova_3dspherical_3dgrid_testrun

getatomicdata atomicdata_feconi.tar.xz

mkdir -p $modelpath

cd $modelpath

tar -xf ../atomicdata_feconi.tar.xz --directory ./

# same input files as the other test run
rsync -av ../kilonova_3dspherical_inputfiles/ ./

# the checksum files of this test
rsync -av --ignore-times ../kilonova_3dspherical_3dgrid_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_kilonova_lte.h artisoptions.h

setoption NUM_PACKETS '3.2e5'

# map the 3D spherical model to a 3D Cartesian propagation grid. Each model cell then holds the midpoint of more
# than one propagation cell
setoption GRID_TYPE_OVERRIDE 'GridType::CARTESIAN3D'
setoption CUBOID_NCOORDGRID_X '24'
setoption CUBOID_NCOORDGRID_Y '24'
setoption CUBOID_NCOORDGRID_Z '24'

setoption RATECOEFF_TABLESIZE '20'
setoption MINTEMP '1000.'
setoption MAXTEMP '20000.'

setoption KEEP_ESCAPED_GAMMAS 'true'

setoption USE_LUT_PHOTOION 'true'

cd -

set +x
