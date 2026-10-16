#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

modelpath=classicmode_1d_3dgrid_testrun

getatomicdata atomicdata_classic.tar.xz

mkdir -p $modelpath

cd $modelpath

tar -xf ../atomicdata_classic.tar.xz --directory ./

rsync -av ../classicmode_1d_3dgrid_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_classic.h artisoptions.h

setoption NUM_PACKETS '6e4'

setoption GRID_TYPE_OVERRIDE 'GridType::CARTESIAN3D'

setoption VPKT_ON 'true'
setoption VPKT_WRITE_CONTRIBS 'true'

setoption KEEP_ESCAPED_GAMMAS 'true'

cd -

set +x
