#!/usr/bin/env zsh

set -x

rsync -av nebularonezone_inputfiles/ nebularonezone_testrun/

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory nebularonezone_testrun/

cp ../data/* nebularonezone_testrun/

set +x
