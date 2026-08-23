#!/usr/bin/env zsh

set -x

runfolder=kilonova_2d_expansionopac_testrun

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O -L https://github.com/artis-mcrt/artis/releases/download/v2026.5.15/atomicdata_feconi.tar.xz; fi

mkdir -p $runfolder

cd $runfolder

tar -xf ../atomicdata_feconi.tar.xz --directory ./

# same input files as the other test run
rsync -av ../kilonova_2d_inputfiles/ ./

# for the checksum files
rsync -av --ignore-times ../kilonova_2d_expansionopac_inputfiles/ ./

ln -s ../../ artis

cp artis/artisoptions_kilonova_lte.h artisoptions.h

xz -f -d -v -T0 *.xz

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 80000;/g' artisoptions.h

sed -i'' -e 's/constexpr int RATECOEFF_TABLESIZE.*/constexpr int RATECOEFF_TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 1000.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 20000.;/g' artisoptions.h


sed -i'' -e 's/constexpr bool RPKT_USE_EXPANSION_OPACITIES.*/constexpr bool RPKT_USE_EXPANSION_OPACITIES = true;/g' artisoptions.h
sed -i'' -e 's/constexpr std::optional<float> RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.*/constexpr std::optional<float> RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY = 1.;/g' artisoptions.h

cd -

set +x
