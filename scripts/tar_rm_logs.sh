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

if compgen -G 'logfiles.tar*' > /dev/null; then
    echo "logfiles.tar* already exists! Not overwriting"
    exit 1
else
    # the artis folder holds the code and the data files of the code, so do not change them
    find . -mindepth 2 ! -path './artis/*' -name "output_*.txt*" -print > $tmpdir/filelist.txt
    tar -cvf $tmpdir/logfiles.tar --files-from $tmpdir/filelist.txt && zstd -v -T0 -15 $tmpdir/logfiles.tar && mv -v $tmpdir/logfiles.tar.zst . && rm -rfv $tmpdir && find . -mindepth 2 ! -path './artis/*' -name "output_*.txt*" ! -name "output_0-0.txt*" -delete
fi
