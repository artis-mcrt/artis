#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

modelpath=kilonova_2d_testrun

getatomicdata atomicdata_feconi.tar.xz

mkdir -p $modelpath

cd $modelpath

tar -xf ../atomicdata_feconi.tar.xz --directory ./

rsync -av ../kilonova_2d_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_kilonova_lte.h artisoptions.h

setoption NUM_PACKETS '3.2e5'

setoption RATECOEFF_TABLESIZE '20'
setoption MINTEMP '1000.'
setoption MAXTEMP '20000.'

setoption KEEP_ESCAPED_GAMMAS 'true'

setoption USE_LUT_PHOTOION 'true'

cd -

set +x
