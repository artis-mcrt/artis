#!/usr/bin/env zsh

set -x

rsync -av nebularonezone_inputfiles/ nebularonezone_testrun/

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory nebularonezone_testrun/

cp ../data/* nebularonezone_testrun/

cp ../artisoptions_nltenebular.h nebularonezone_testrun/artisoptions.h

cd nebularonezone_testrun

# sed -i 's/#define MPKTS.*/#define MPKTS 15000/g' artisoptions.h

sed -i 's/#define TABLESIZE.*/#define TABLESIZE 20/g' artisoptions.h
sed -i 's/#define MINTEMP.*/#define MINTEMP 2000./g' artisoptions.h
sed -i 's/#define MAXTEMP.*/#define MAXTEMP 10000./g' artisoptions.h

sed -i 's/static const int FIRST_NLTE_RADFIELD_TIMESTEP.*/static const int FIRST_NLTE_RADFIELD_TIMESTEP = 7;/g' artisoptions.h

sed -i 's/#define DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP.*/#define DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP 7/g' artisoptions.h

sed -i 's/#define SF_AUGER_CONTRIBUTION_ON.*/#define SF_AUGER_CONTRIBUTION_ON false/g' artisoptions.h

sed -i 's/static bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC.*/static bool WRITE_PARTIAL_EMISSIONABSORPTIONSPEC = true;/g' artisoptions.h

cd -

set +x
