#!/usr/bin/env zsh

set -x

runfolder=kilonova_2d_2dgrid_barnesthermalisation_testrun

mkdir -p $runfolder

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory $runfolder/

# same input files as the other test run
rsync -av kilonova_2d_3dgrid_inputfiles/ $runfolder/

# for the checksum files
rsync -av --ignore-times kilonova_2d_2dgrid_barnesthermalisation_inputfiles/ $runfolder/

cp ../data/* $runfolder/

cp ../artisoptions_kilonova_lte.h $runfolder/artisoptions.h

cd $runfolder

xz -dv -T0 *.xz

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 80000;/g' artisoptions.h

sed -i'' -e 's/constexpr int GRID_TYPE.*/constexpr int GRID_TYPE = GRID_CYLINDRICAL2D;/g' artisoptions.h

sed -i'' -e 's/constexpr int TABLESIZE.*/constexpr int TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 1000.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 20000.;/g' artisoptions.h

sed -i'' -e 's/constexpr auto PARTICLE_THERMALISATION_SCHEME.*/constexpr auto PARTICLE_THERMALISATION_SCHEME = ThermalisationScheme::BARNES;/g' artisoptions.h

sed -i'' -e 's/constexpr auto GAMMA_THERMALISATION_SCHEME.*/constexpr auto GAMMA_THERMALISATION_SCHEME = ThermalisationScheme::BARNES;/g' artisoptions.h

cd -

set +x
