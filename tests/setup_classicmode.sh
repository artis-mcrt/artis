#!/usr/bin/env zsh

set -x

rsync -av classicmode_inputfiles/ classicmode_testrun/

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory classicmode_testrun/

cp ../data/* classicmode_testrun/
# cp ../artisoptions_classic.h classicmode_testrun/artisoptions.h

sed 's/#define MPKTS.*/#define MPKTS 15000/g' ../artisoptions_classic.h > classicmode_testrun/artisoptions.h

set +x
