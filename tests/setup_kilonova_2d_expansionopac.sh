#!/usr/bin/env zsh

set -x

source ./setupfuncs.sh

runfolder=kilonova_2d_expansionopac_testrun

getatomicdata atomicdata_feconi.tar.xz

mkdir -p $runfolder

cd $runfolder

tar -xf ../atomicdata_feconi.tar.xz --directory ./

# same input files as the other test run
rsync -av ../kilonova_2d_inputfiles/ ./

# the phixsdata_v2.txt of this test, and the checksum files
rsync -av --ignore-times ../kilonova_2d_expansionopac_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_kilonova_lte.h artisoptions.h

xz -f -d -v -T0 *.xz

sedopt "constexpr std::int64_t NUM_PACKETS.*" "constexpr std::int64_t NUM_PACKETS = 320'000;"

sedopt 'constexpr int RATECOEFF_TABLESIZE.*' 'constexpr int RATECOEFF_TABLESIZE = 20;'
sedopt 'constexpr double MINTEMP.*' 'constexpr double MINTEMP = 1000.;'
sedopt 'constexpr double MAXTEMP.*' 'constexpr double MAXTEMP = 20000.;'

sedopt 'constexpr bool RPKT_USE_EXPANSION_OPACITIES.*' 'constexpr bool RPKT_USE_EXPANSION_OPACITIES = true;'
sedopt 'constexpr std::optional<float> RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.*' 'constexpr std::optional<float> RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY = 1.;'

rm -f artisoptions.h.bak

cd -

set +x
