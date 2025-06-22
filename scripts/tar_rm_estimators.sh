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

find . -type d -name "*.slurm" -print0 | while IFS= read -r -d '' dir; do
    if [[ -d "$dir" ]]; then
        echo "runfolder: $dir"
        cd "$dir"
        if [ -f estimators_allranks.tar* ]; then
            echo "  estimators_allranks.tar*  exists already!"
            ls -lh estimators_allranks.tar*
        else
            if [ -e estimbatch00_*.parquet* -a -e estimators_0001.out* ]; then
                find . -mindepth 0 -name "estimators_*.out*" -print | sort > $tmpdir/estimatorfilelist.txt
                echo "  Creating tarball of estimators_allranks.tar.zst"
                tar -cf $tmpdir/estimators_allranks.tar --files-from $tmpdir/estimatorfilelist.txt && zstd -v -T0 -15 $tmpdir/estimators_allranks.tar && mv -v $tmpdir/estimators_allranks.tar.zst . && rm -f $tmpdir/* && find . -mindepth 0 -name "estimators_*.out*" ! -name "estimators_0000.out*" -delete
                # copy the timestamp of the first estimator file to the tarball
                #touch -m -r estimators_0000.out* estimators_allranks.tar.zst
            else
                echo "  no parquet batches or estimators_0001.out* file found, skipping tar creation"
            fi
        fi
        cd - > /dev/null
    fi
done
