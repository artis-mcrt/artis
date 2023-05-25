#!/usr/bin/env zsh

set -x

rsync -av nebularonezone_inputfiles/ nebularonezone_testrun/

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory nebularonezone_testrun/

cp ../data/* nebularonezone_testrun/

cp ../artisoptions_nltenebular.h nebularonezone_testrun/artisoptions.h

cd nebularonezone_testrun

sed -i'' -e 's/constexpr int MPKTS.*/constexpr int MPKTS = 1000000;/g' artisoptions.h

sed -i'' -e 's/constexpr int TABLESIZE.*/constexpr int TABLESIZE = 20;/g' artisoptions.h
sed -i'' -e 's/constexpr double MINTEMP.*/constexpr double MINTEMP = 2000.;/g' artisoptions.h
sed -i'' -e 's/constexpr double MAXTEMP.*/constexpr double MAXTEMP = 10000.;/g' artisoptions.h

sed -i'' -e 's/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP.*/constexpr int FIRST_NLTE_RADFIELD_TIMESTEP = 7;/g' artisoptions.h

sed -i'' -e 's/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP.*/constexpr int DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP = 7;/g' artisoptions.h

sed -i'' -e 's/constexpr bool SF_AUGER_CONTRIBUTION_ON.*/constexpr bool SF_AUGER_CONTRIBUTION_ON = false;/g' artisoptions.h

sed -i'' -e 's/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/constexpr bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' artisoptions.h

cd -

set +x
