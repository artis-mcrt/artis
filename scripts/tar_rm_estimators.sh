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
            # The files of the ranks go into the archive only when a file of all ranks holds the same text. artistools
            # then reads that file for a new conversion. combine_estimator_files.py makes the file. A file of all ranks
            # that is older than a file of a rank does not hold the later text, e.g. when the script ran during the job.
            allranksfile=""
            for candidate in estimators_allranks.out estimators_allranks.out.zst estimators_allranks.out.gz estimators_allranks.out.xz; do
                if [ -f "$candidate" ]; then
                    allranksfile=$candidate
                    break
                fi
            done
            if ! compgen -G 'estimators_[0-9]*.out*' > /dev/null; then
                echo "  The folder has no estimator file of a rank, so the script makes no tar file."
            elif [ -z "$allranksfile" ]; then
                echo "  The folder has no estimators_allranks.out file. Run combine_estimator_files.py first. The script keeps the files."
            elif [ -n "$(find . -maxdepth 1 -name 'estimators_[0-9]*.out*' -newer "$allranksfile" -print -quit)" ]; then
                echo "  $allranksfile is older than a file of a rank. Combine the files of the ranks again. The script keeps the files."
            else
                find . -mindepth 0 -name "estimators_[0-9]*.out*" -print | sort > $tmpdir/estimatorfilelist.txt
                echo "  Creating tarball of estimators_allranks.tar"
                tar -cf $tmpdir/estimators_allranks.tar --files-from $tmpdir/estimatorfilelist.txt && mv -v $tmpdir/estimators_allranks.tar . && rm -f $tmpdir/* && find . -mindepth 0 -name "estimators_[0-9]*.out*" -delete
            fi
        fi
        cd - > /dev/null
    fi
done
