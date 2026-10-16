#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

modelpath=nebular_1d_3dgrid_testrun

getatomicdata atomicdata_feconi.tar.xz

mkdir -p $modelpath

cd $modelpath

rsync -av ../nebular_1d_3dgrid_inputfiles/ ./

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

setoption WRITE_ION_HEATING_COOLING_RATES 'true'

cd -

set +x
