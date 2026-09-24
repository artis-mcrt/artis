#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

runfolder=kilonova_1d_testrun

getatomicdata atomicdata_feconi.tar.xz

mkdir -p $runfolder

cd $runfolder

tar -xf ../atomicdata_feconi.tar.xz --directory ./

rsync -av ../kilonova_1d_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_kilonova_lte.h artisoptions.h

xz -f -d -v -T0 *.xz

sedopt "constexpr std::int64_t NUM_PACKETS.*" "constexpr std::int64_t NUM_PACKETS = 320'000;"

sedopt 'constexpr int RATECOEFF_TABLESIZE.*' 'constexpr int RATECOEFF_TABLESIZE = 20;'
sedopt 'constexpr double MINTEMP.*' 'constexpr double MINTEMP = 1000.;'
sedopt 'constexpr double MAXTEMP.*' 'constexpr double MAXTEMP = 20000.;'

sedopt 'constexpr bool KEEP_ESCAPED_GAMMAS.*' 'constexpr bool KEEP_ESCAPED_GAMMAS = true;'

sedopt 'constexpr bool USE_LUT_PHOTOION.*' 'constexpr bool USE_LUT_PHOTOION = true;'

rm -f artisoptions.h.bak

cd -

set +x
