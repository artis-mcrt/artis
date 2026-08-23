#!/usr/bin/env zsh

set -x

runfolder=classicmode_3d_testrun

if [ ! -f atomicdata_classic.tar.xz ]; then curl -O -L https://github.com/artis-mcrt/artis/releases/download/v2026.5.15/atomicdata_classic.tar.xz; fi

mkdir -p $runfolder

cd $runfolder

tar -xf ../atomicdata_classic.tar.xz --directory ./

rsync -av ../classicmode_3d_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_classic.h artisoptions.h

xz -f -d -v -T0 *.xz

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 15000;/g' artisoptions.h


sed -i'' -e 's/constexpr bool VPKT_ON.*/constexpr bool VPKT_ON = true;/g' artisoptions.h
sed -i'' -e 's/constexpr bool VPKT_WRITE_CONTRIBS.*/constexpr bool VPKT_WRITE_CONTRIBS = true;/g' artisoptions.h
sed -i'' -e 's/constexpr bool VPKT_USE_EXPANSION_OPACITIES.*/constexpr bool VPKT_USE_EXPANSION_OPACITIES = true;/g' artisoptions.h

cd -

set +x
