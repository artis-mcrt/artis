#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

runfolder=classicmode_3d_testrun

getatomicdata atomicdata_classic.tar.xz

mkdir -p $runfolder

cd $runfolder

tar -xf ../atomicdata_classic.tar.xz --directory ./

rsync -av ../classicmode_3d_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_classic.h artisoptions.h

xz -f -d -v -T0 *.xz

sedopt "constexpr std::int64_t NUM_PACKETS.*" "constexpr std::int64_t NUM_PACKETS = 60'000;"

sedopt 'constexpr bool VPKT_ON.*' 'constexpr bool VPKT_ON = true;'
sedopt 'constexpr bool VPKT_WRITE_CONTRIBS.*' 'constexpr bool VPKT_WRITE_CONTRIBS = true;'
sedopt 'constexpr bool VPKT_USE_EXPANSION_OPACITIES.*' 'constexpr bool VPKT_USE_EXPANSION_OPACITIES = true;'

sedopt 'constexpr bool KEEP_ESCAPED_GAMMAS.*' 'constexpr bool KEEP_ESCAPED_GAMMAS = true;'

rm -f artisoptions.h.bak

cd -

set +x
