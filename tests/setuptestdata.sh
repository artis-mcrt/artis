#!/usr/bin/env zsh

set -x

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory nebularonezone/

cp ../data/* nebularonezone/

if [ ! -f nebularonezone_reference_20211109.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/nebularonezone_reference_20211109.tar.xz; fi

tar -xf nebularonezone_reference_20211109.tar.xz
tar -xf atomicdata_feconi.tar.xz --directory nebularonezone_reference/
cp ../data/* nebularonezone_reference/

set +x
