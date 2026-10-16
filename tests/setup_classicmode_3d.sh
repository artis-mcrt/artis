#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

modelpath=classicmode_3d_testrun

getatomicdata atomicdata_classic.tar.xz

mkdir -p $modelpath

cd $modelpath

tar -xf ../atomicdata_classic.tar.xz --directory ./

rsync -av ../classicmode_3d_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_classic.h artisoptions.h

setoption NUM_PACKETS '6e4'

setoption VPKT_ON 'true'
setoption VPKT_WRITE_CONTRIBS 'true'
setoption VPKT_USE_BINNED_OPACITIES 'true'

setoption KEEP_ESCAPED_GAMMAS 'true'

cd -

set +x
