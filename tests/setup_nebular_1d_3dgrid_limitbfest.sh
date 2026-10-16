#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

modelpath=nebular_1d_3dgrid_limitbfest_testrun

getatomicdata atomicdata_feconi.tar.xz

mkdir -p $modelpath

cd $modelpath

rsync -av ../nebular_1d_3dgrid_inputfiles/ ./

rsync --ignore-times -av ../nebular_1d_3dgrid_limitbfest_inputfiles/ ./

tar -xf ../atomicdata_feconi.tar.xz --directory .

ln -s ../../ artis

cp artis/artisoptions_nltenebular.h artisoptions.h

setoption NUM_PACKETS '4e6'

setoption GRID_TYPE_OVERRIDE 'GridType::CARTESIAN3D'

setoption RATECOEFF_TABLESIZE '20'
setoption MINTEMP '2000.'
setoption MAXTEMP '10000.'

setoption MULTIBIN_RADFIELD_FIRST_TIMESTEP '7'

setoption DETAILED_BF_ESTIMATORS_FIRST_TIMESTEP '7'

setoption SF_AUGER_CONTRIBUTION_ON 'false'

setoption DETAILED_BF_ESTIMATORS_LEVEL_SELECTED '{ return level <= ION_NLEVELS_EXCITED_NLTE(element_z, ionstage); }'

cd -

set +x
