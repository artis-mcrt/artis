#!/usr/bin/env zsh

set -x

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://psweb.mp.qub.ac.uk/artis/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory nebularonezone/

tar -xf atomicdata_feconi.tar.xz --directory nebularonezone_reference/

cp ../data/* nebularonezone/
cp ../data/* nebularonezone_reference/

if [ ! -f nebularonezone_reference_20210324.tar.xz ]; then curl -O https://psweb.mp.qub.ac.uk/artis/artis/nebularonezone_reference_20210324.tar.xz; fi

tar -xf nebularonezone_reference_20210324.tar.xz

set +x