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

if [[ -f logfiles.tar.zst ]]; then
    echo "logfiles.tar.zst already exists! Not overwriting"
    exit 1
else
    tar --zstd -cvf $tmpdir/logfiles.tar.zst */output_*.txt* && mv $tmpdir/logfiles.tar.zst . && rm -rf $tmpdir && find . -mindepth 2 -name "output_*.txt*" ! -name "output_0-0.txt*" -delete

fi

