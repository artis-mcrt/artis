#!/usr/bin/env zsh

set -x

rsync -av classicmode_inputfiles/ classicmode_testrun/

if [ ! -f atomicdata_feconi.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/atomicdata_feconi.tar.xz; fi

tar -xf atomicdata_feconi.tar.xz --directory classicmode_testrun/

cp ../data/* classicmode_testrun/
# cp ../artisoptions_classic.h classicmode_testrun/artisoptions.h
sed 's/#define MPKTS.*/#define MPKTS 5000/g' ../artisoptions_classic.h > classicmode_testrun/artisoptions.h

#if [ ! -f classicmode_reference_20211126.tar.xz ]; then curl -O https://theory.gsi.de/~lshingle/artis_http_public/artis/classicmode_reference_20211126.tar.xz; fi

#tar -xf classicmode_reference_20211126.tar.xz
#tar -xf atomicdata_feconi.tar.xz --directory classicmode_reference/
#cp ../data/* classicmode_reference/
#if [ ! -f classicmode_reference/input.txt ]; then cp classicmode_inputfiles/input-newrun.txt classicmode_reference/input.txt; fi

set +x
