#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

modelpath=classicmode_1d_3dgrid_thermalinput_testrun

getatomicdata atomicdata_classic.tar.xz

mkdir -p $modelpath

cd $modelpath

tar -xf ../atomicdata_classic.tar.xz --directory ./

# same input files as the other test run
rsync -av ../classicmode_1d_3dgrid_inputfiles/ ./

# energydistribution.txt and energyrate.txt of this test, and the checksum files
rsync -av --ignore-times ../classicmode_1d_3dgrid_thermalinput_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_classic.h artisoptions.h

setoption NUM_PACKETS '6e4'

setoption GRID_TYPE_OVERRIDE 'GridType::CARTESIAN3D'

setoption INITIAL_PACKETS_ON 'false'
setoption THERMAL_ENERGY_INPUT_FROM_FILES 'true'

cd -

set +x
