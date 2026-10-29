#!/usr/bin/env bash
# Run an end-to-end test. The test job of .github/workflows/ci.yml calls this script for each step.

set -eo pipefail

testsdir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repopath=$(dirname "$testsdir")
testname=$1
step=${2:-all}
modelpath="$testsdir/${testname}_testrun"
testmode=${TESTMODE:-OFF}
mpirun_command=(mpirun --oversubscribe --mca osc '^ucx' --mca mpi_yield_when_idle 1)

# the sanitizers of TESTMODE=ON need these settings
if [ "$testmode" = ON ]; then
    export ASAN_OPTIONS=${ASAN_OPTIONS:-detect_stack_use_after_return=1:detect_leaks=0}
    export OMPI_MCA_memory='^patcher'
fi

usage() {
    cat >&2 <<EOF
Usage: tests/runtest.sh <testname> [step]

Without a step, the script runs these steps in order. It stops at the first step that fails.
  setup            remove tests/<testname>_testrun and make it again with tests/setup_<testname>.sh
  build            build sn3d and exspec from the artisoptions.h of the model folder, and copy them into it
  job0             start a new simulation
  checksums_job0   compare the checksums of the output files with results_md5_job0.txt (only directly after job0)
  resume           resume the simulation from the restart files
  exspec           remove the restart files, then run exspec and scripts/mergeangleres.py
  checksums_final  compare the checksums of the output files with results_md5_final.txt

TESTMODE=ON in the environment adds the sanitizers and the extra assertions to the build.

Tests:
EOF
    for setupscript in "$testsdir"/setup_*.sh; do
        setupscript=${setupscript##*/setup_}
        echo "  ${setupscript%.sh}" >&2
    done
    exit 1
}

if [ -z "$testname" ] || [ ! -f "$testsdir/setup_${testname}.sh" ]; then
    usage
fi

if [ "$step" != all ] && [ "$step" != setup ] && [ ! -d "$modelpath" ]; then
    echo "[error] The model folder $modelpath is absent. Run the setup step first." >&2
    exit 1
fi

case $step in
    all)
        for step_name in setup build job0 checksums_job0 resume exspec checksums_final; do
            if ! "${BASH_SOURCE[0]}" "$testname" "$step_name"; then
                echo "[error] The step $step_name of the test $testname failed." >&2
                exit 1
            fi
        done
        echo "All steps of the test $testname succeeded. The model folder is $modelpath."
        ;;

    setup)
        rm -rf "$modelpath"
        (cd "$testsdir" && source "./setup_${testname}.sh")
        ;;

    build)
        cd "$repopath"
        saved_options=artisoptions.h.saved_by_runtest
        if [ -L "$saved_options" ] || [ -e "$saved_options" ]; then
            echo "[error] An earlier test did not restore $saved_options. Rename it to artisoptions.h." >&2
            exit 1
        fi
        # mv keeps the timestamp, so a later build of your own options does not compile everything again
        if [ -L artisoptions.h ] || [ -e artisoptions.h ]; then
            mv artisoptions.h "$saved_options"
            trap 'mv -f "$saved_options" artisoptions.h' EXIT
        else
            trap 'rm -f artisoptions.h' EXIT
        fi
        ln -s "tests/${testname}_testrun/artisoptions.h" artisoptions.h
        make --check-symlink-times --jobs="$(getconf _NPROCESSORS_ONLN)" REPRODUCIBLE=ON ZSTD=ON TESTMODE="$testmode" \
            MAX_NODE_SIZE=2 FASTMATH=OFF sn3d exspec
        cp -L sn3d exspec "$modelpath/"
        # the objects hold the test options, so a later build in this build folder must compile them again
        touch -t 200001010000 "$(dirname "$(readlink sn3d)")"/*.o
        ;;

    job0)
        cd "$modelpath"
        # without input.txt, sn3d starts a new simulation from input-newrun.txt
        rm -f input.txt
        time "${mpirun_command[@]}" -np 4 ./sn3d
        ;;

    resume)
        cd "$modelpath"
        cp input-resume.txt input.txt
        time "${mpirun_command[@]}" -np 4 ./sn3d
        ;;

    exspec)
        cd "$modelpath"
        rm -f ./*.tmp
        ls -l packets/
        time "${mpirun_command[@]}" -np 2 ./exspec
        python3 "$repopath/scripts/mergeangleres.py"
        ;;

    checksums_job0 | checksums_final)
        cd "$modelpath"
        checksum_set=${step#checksums_}
        checksumfile=results_md5_${checksum_set}.txt
        # The file in the model folder exists only if md5sum wrote the full list. The setup step copies the
        # reference file into the model folder, so this removal comes first.
        rm -f "$checksumfile" "$checksumfile.part"
        if [ "$checksum_set" = job0 ]; then
            jobfolder=job_from_ts0000
        else
            # the job folder with the highest start timestep is the one of the resume run
            jobfolder=$(ls -d job_from_ts* | grep -v ts0000 | tail -n 1)
        fi
        reffile="$testsdir/${testname}_inputfiles/$checksumfile"
        # The checksums are of the decompressed content under the plain name, in the output format of md5sum.
        shopt -s nullglob
        for zstfile in *.out.zst "$jobfolder"/*.out.zst packets/*.out.zst vpackets/*.out.zst vspecpol/*.out.zst \
            vpkt_grid/*.out.zst speclc_angle_res/*.zst; do
            zstdcat "$zstfile" | md5sum | sed "s|  -$|  ${zstfile%.zst}|"
        done | tee "$checksumfile.part"
        shopt -u nullglob
        mv "$checksumfile.part" "$checksumfile"
        if [ ! -s "$reffile" ]; then
            echo "${GITHUB_ACTIONS:+::warning::}The repository has no $checksumfile for $testname. The step compares nothing."
        elif diff "$reffile" "$checksumfile"; then
            echo "The checksums agree with ${reffile#"$repopath"/}."
        else
            echo "[error] The checksums differ from ${reffile#"$repopath"/}. In the diff, '<' marks the reference and '>' marks this run." >&2
            exit 1
        fi
        ;;

    *)
        usage
        ;;
esac
