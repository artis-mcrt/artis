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

find . -type d \( -name "*.slurm" -o -name "job_from_ts*" \) -print0 | while IFS= read -r -d '' dir; do
    if [[ -d "$dir" ]]; then
        echo "runfolder: $dir"
        cd "$dir"
        if [ -f estimators_allranks.tar* ]; then
            echo "  estimators_allranks.tar*  exists already!"
            ls -lh estimators_allranks.tar*
        else
            # A file of all ranks that is older than a file of a rank does not hold the later text, e.g. when
            # combine_estimator_files.py ran while the job ran. The files of the ranks must then stay, because
            # artistools reads the file of all ranks when the files of the ranks are absent.
            stale_allranks_file=""
            for allranksfile in estimators_allranks.out estimators_allranks.out.zst estimators_allranks.out.gz estimators_allranks.out.xz; do
                if [ -f "$allranksfile" ] && [ -n "$(find . -maxdepth 1 -name 'estimators_[0-9]*.out*' -newer "$allranksfile" -print -quit)" ]; then
                    stale_allranks_file=$allranksfile
                fi
            done
            # A parquet cache that is older than a file of a rank does not hold the later timesteps, e.g. when the job
            # wrote more timesteps after the conversion. The files of the ranks must then stay for a new conversion.
            stale_cache=""
            for cachefile in estimators_allranks.out.parquet* estimbatch*.parquet*; do
                if [ -f "$cachefile" ] && [ -n "$(find . -maxdepth 1 -name 'estimators_[0-9]*.out*' -newer "$cachefile" -print -quit)" ]; then
                    stale_cache=$cachefile
                fi
            done
            if [ -n "$stale_allranks_file" ]; then
                echo "  $stale_allranks_file is older than a file of a rank. Combine the files of the ranks again. The script keeps the files."
            elif [ -n "$stale_cache" ]; then
                echo "  $stale_cache is older than a file of a rank. Read the estimators with artistools again. The script keeps the files."
            # artistools writes the cache estimators_allranks.out.parquet, and an earlier version wrote estimbatch*.parquet*
            elif (compgen -G 'estimators_allranks.out.parquet*' > /dev/null || compgen -G 'estimbatch00_*.parquet*' > /dev/null) && compgen -G 'estimators_[0-9]*.out*' > /dev/null; then
                find . -mindepth 0 -name "estimators_[0-9]*.out*" -print | sort > $tmpdir/estimatorfilelist.txt
                echo "  Creating tarball of estimators_allranks.tar"
                tar -cf $tmpdir/estimators_allranks.tar --files-from $tmpdir/estimatorfilelist.txt && mv -v $tmpdir/estimators_allranks.tar . && rm -f $tmpdir/* && find . -mindepth 0 -name "estimators_[0-9]*.out*" -delete
            else
                echo "  The folder has no parquet cache or no estimator file of a rank, so the script makes no tar file."
            fi
        fi
        cd - > /dev/null
    fi
done
