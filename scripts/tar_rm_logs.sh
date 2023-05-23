#!/usr/bin/env bash

tmpdir=$(mktemp -d /tmp/XXXXXX)

function trap_ctrlc() {
   if [[ -d "$tmpdir" ]]; then
     rm -rf $tmpdir
     echo "\nCtrl-C caught...deleted temp dir: $tmpdir"
   fi

   exit 2
}

trap "trap_ctrlc" 2

if [ -f logfiles.tar* ]; then
    echo "logfiles.tar* already exists! Not overwriting"
    exit 1
else
    tar -cvf $tmpdir/logfiles.tar */output_*.txt* && zstd -v -T0 -15 $tmpdir/logfiles.tar && mv -v $tmpdir/logfiles.tar.zst . && rm -rfv $tmpdir && find . -mindepth 2 -name "output_*.txt*" ! -name "output_0-0.txt*" -delete
fi

