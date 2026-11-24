#!/usr/bin/env bash
# Run an end-to-end test. The test job of .github/workflows/ci.yml calls this script for each step.

set -eo pipefail

# the checksum steps need the same order of the file names in each locale
export LC_ALL=C

testsdir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repopath=$(dirname "$testsdir")
testname=$1
step=${2-all}
modelpath="$testsdir/${testname}_testrun"
testmode=${TESTMODE:-OFF}
mpirun_command=(mpirun --oversubscribe --mca osc '^ucx' --mca mpi_yield_when_idle 1)
steps=(setup build inputcheck formatcheck job0 checksums_job0 resume exspec checksums_final)

# the sanitizers of TESTMODE=ON need these settings
if [ "$testmode" = ON ]; then
    export ASAN_OPTIONS=${ASAN_OPTIONS:-detect_stack_use_after_return=1:detect_leaks=0}
    export OMPI_MCA_memory='^patcher'
fi

test_names() {
    for setupscript in "$testsdir"/setup_*.sh; do
        setupscript=${setupscript##*/setup_}
        echo "${setupscript%.sh}"
    done
}

usage() {
    cat <<EOF
Usage: tests/runtest.sh <testname> [step]
       tests/runtest.sh --list
       tests/runtest.sh --help

Without a step, the script runs these steps in order. It stops at the first step that fails, except a checksum
step that finds different checksums. The other steps then run, and at the end the script lists the files that
differ and exits with an error status.
  setup            remove tests/<testname>_testrun and make it again with tests/setup_<testname>.sh
  build            build sn3d and exspec from the artisoptions.h of the model folder, and copy them into it
  inputcheck       check the input of the new simulation with sn3d -c
  formatcheck      check the reader of model.txt of format version 2 with sn3d -c, for the small models in
                   tests/<testname>_inputfiles/formatcheck_models. A test with no such folder checks nothing
  job0             start a new simulation
  checksums_job0   compare the checksums of the output files with results_md5_job0.txt (only directly after job0)
  resume           check the input and the restart files with sn3d -c, then resume the simulation
  exspec           remove the restart files, then run exspec and scripts/mergeangleres.py
  checksums_final  compare the checksums of the output files with results_md5_final.txt

TESTMODE=ON in the environment adds the sanitizers and the extra assertions to the build.

Tests:
EOF
    test_names | sed 's/^/  /'
}

case $testname in
    -h | --help)
        usage
        exit 0
        ;;
    --list)
        test_names
        exit 0
        ;;
    "")
        usage >&2
        exit 1
        ;;
esac

if [ ! -f "$testsdir/setup_${testname}.sh" ]; then
    echo "[error] The script has no test '$testname'. Run tests/runtest.sh --list to see the tests." >&2
    exit 1
fi

step_is_known=false
for step_name in all "${steps[@]}"; do
    if [ "$step" = "$step_name" ]; then
        step_is_known=true
    fi
done
if [ "$step_is_known" = false ]; then
    echo "[error] The script has no step '$step'. Run tests/runtest.sh --help to see the steps." >&2
    exit 1
fi

if [ "$step" != all ] && [ "$step" != setup ] && [ ! -d "$modelpath" ]; then
    echo "[error] The model folder $modelpath is absent. Run the setup step first." >&2
    exit 1
fi

# Run the input check of sn3d with $1 ranks. The check must change no file and no folder of the model folder outside
# its own folder inputcheck.
run_inputcheck() {
    list_files_and_folders() { find . -path ./inputcheck -prune -o -print | sort; }
    local files_and_folders_before files_changed time_before
    files_and_folders_before=$(list_files_and_folders)
    time_before=$(mktemp)
    sleep 1  # some file systems store the time of a change in whole seconds
    "${mpirun_command[@]}" -np "$1" ./sn3d -c
    files_changed=$(find . -path ./inputcheck -prune -o ! -path . -newer "$time_before" -print)
    rm -f "$time_before"
    if [ "$files_and_folders_before" != "$(list_files_and_folders)" ] || [ -n "$files_changed" ]; then
        echo "[error] The input check changed the model folder: $files_changed" >&2
        exit 1
    fi
    rm -rf inputcheck
}

case $step in
    all)
        checksum_steps_with_differences=()
        failed_step=
        for step_name in "${steps[@]}"; do
            step_status=0
            "$testsdir/runtest.sh" "$testname" "$step_name" || step_status=$?
            if [ "$step_status" -eq 0 ]; then
                continue
            fi
            # A different compiler or CPU type can give different checksums, so the other steps still run. A checksum
            # step gives the status 1 with its checksum file only for a difference of the checksums.
            if [[ $step_name == checksums_* ]] && [ "$step_status" -eq 1 ] &&
                [ -f "$modelpath/results_md5_${step_name#checksums_}.txt" ]; then
                echo "[error] The step $step_name of the test $testname failed." >&2
                checksum_steps_with_differences+=("$step_name")
            else
                failed_step=$step_name
                break
            fi
        done
        for step_name in "${checksum_steps_with_differences[@]}"; do
            checksumfile=results_md5_${step_name#checksums_}.txt
            echo "[error] These output files differ from tests/${testname}_inputfiles/$checksumfile. Each is changed, missing, or extra:" >&2
            # diff exits with status 1 when the files differ
            diff "$testsdir/${testname}_inputfiles/$checksumfile" "$modelpath/$checksumfile" |
                sed -n 's/^[<>] [0-9a-f]*  /  /p' | sort -u >&2 || true
        done
        if [ -n "$failed_step" ]; then
            echo "[error] The step $failed_step of the test $testname failed." >&2
            exit 1
        fi
        if [ ${#checksum_steps_with_differences[@]} -gt 0 ]; then
            echo "[error] The test $testname failed only at the comparison of the checksums. The model folder is $modelpath." >&2
            exit 1
        fi
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
            echo "[error] $saved_options exists. Another build step runs now, or an earlier build step did not restore it. If no build step runs, rename it to artisoptions.h." >&2
            exit 1
        fi
        # the build folder of the test options has its own name, so your own build compiles nothing again
        if [ -L artisoptions.h ] || [ -e artisoptions.h ]; then
            mv artisoptions.h "$saved_options"
            trap 'mv -f "$saved_options" artisoptions.h' EXIT
        else
            trap 'rm -f artisoptions.h' EXIT
        fi
        ln -s "tests/${testname}_testrun/artisoptions.h" artisoptions.h
        make --jobs="$(getconf _NPROCESSORS_ONLN)" REPRODUCIBLE=ON ZSTD=ON TESTMODE="$testmode" \
            MAX_NODE_SIZE=2 FASTMATH=OFF OPENMP=OFF STDPAR=OFF GPU=OFF OPTIMIZE=ON PGO= sn3d exspec
        cp -L sn3d exspec "$modelpath/"
        ;;

    inputcheck)
        cd "$modelpath"
        # without input.txt, the check reads input-newrun.txt
        rm -f input.txt
        run_inputcheck 2
        ;;

    formatcheck)
        checkmodels="$testsdir/${testname}_inputfiles/formatcheck_models"
        if [ ! -d "$checkmodels" ]; then
            echo "The test $testname has no folder formatcheck_models. The step checks nothing."
            exit 0
        fi
        # each case runs in its own folder outside the model folder, thus the model folder gets no new file
        checkroot=$(mktemp -d)
        trap 'rm -rf "$checkroot"' EXIT

        # Run sn3d -c for one case. $1 is the name of the case, and $2 is its model.txt. $3 is "complete" for a model
        # that sn3d must read, or the text of the stop message for a model that sn3d must stop for. Each other
        # argument is a text that the output must also have. An argument "file:name" makes the file name in the case
        # folder, with the content 0.
        run_formatcheck_case() {
            local casename=$1 modelfile=$2 expected=$3
            shift 3
            local casefolder="$checkroot/$casename" path status=0 text
            mkdir "$casefolder"
            for path in "$modelpath"/*; do
                case ${path##*/} in
                    model.txt* | input.txt | inputcheck | formatcheck_models) ;;
                    *) ln -s "$path" "$casefolder/" ;;
                esac
            done
            cp "$modelfile" "$casefolder/model.txt"
            local texts=()
            for text in "$@"; do
                if [ "${text#file:}" != "$text" ]; then
                    echo 0 > "$casefolder/${text#file:}"
                else
                    texts+=("$text")
                fi
            done
            (cd "$casefolder" && "${mpirun_command[@]}" -np 2 ./sn3d -c) > "$casefolder/sn3d_check.log" 2>&1 \
                || status=$?
            if [ "$expected" = complete ]; then
                texts+=("input check: complete")
                if [ "$status" != 0 ]; then
                    tail -n 20 "$casefolder/sn3d_check.log" >&2
                    echo "[error] formatcheck $casename: sn3d -c stopped, but it must read this model." >&2
                    exit 1
                fi
            else
                texts+=("$expected")
                if [ "$status" = 0 ]; then
                    echo "[error] formatcheck $casename: sn3d -c read this model, but it must stop." >&2
                    exit 1
                fi
            fi
            for text in "${texts[@]}"; do
                if ! grep -rqF -- "$text" "$casefolder/sn3d_check.log" "$casefolder/inputcheck"; then
                    tail -n 20 "$casefolder/sn3d_check.log" >&2
                    echo "[error] formatcheck $casename: the output of sn3d -c does not have the text '$text'." >&2
                    exit 1
                fi
            done
            echo "formatcheck $casename: correct"
        }

        model1d="$checkmodels/spherical1d_nofegroup_model.txt"
        model2d="$checkmodels/cylindrical2d_model.txt"
        run_formatcheck_case spherical1d_nofegroup "$model1d" complete "model.txt has no column X_Fegroup"
        run_formatcheck_case cylindrical2d "$model2d" complete

        # damaged copies of the 2D model, for which sn3d must stop
        celllinecount=$(grep -c '^[0-9]' "$model2d")
        sed "s/^cell_line_count .*/cell_line_count $((celllinecount - 1))/" "$model2d" > "$checkroot/wrongcount.txt"
        run_formatcheck_case wrongcount "$checkroot/wrongcount.txt" \
            "the key cell_line_count gives $((celllinecount - 1)) cell lines, but the file has $celllinecount"
        # a file that lost its last three lines
        sed -e '$d' "$model2d" | sed -e '$d' | sed -e '$d' > "$checkroot/cutlines.txt"
        run_formatcheck_case cutlines "$checkroot/cutlines.txt" \
            "the key cell_line_count gives $celllinecount cell lines, but the file has $((celllinecount - 3))"
        # a cell line that lost its end
        sed -e '$s/ [^ ]*$//' "$model2d" > "$checkroot/cutvalue.txt"
        run_formatcheck_case cutvalue "$checkroot/cutvalue.txt" "has no readable value for the column"
        # a version line with text after the version
        sed '1s/^format_version 2$/format_version 2 junk/' "$model2d" > "$checkroot/textafterversion.txt"
        run_formatcheck_case textafterversion "$checkroot/textafterversion.txt" "has the format version '2 junk'"
        # a cell line with text after its values
        sed -e '$s/$/ junk/' "$model2d" > "$checkroot/textaftervalues.txt"
        run_formatcheck_case textaftervalues "$checkroot/textaftervalues.txt" "after the values of the columns"
        # a line of a cell with a density of zero that holds a negative mass fraction
        sed -e '$s/^\([0-9]* \)[^ ]* [^ ]* /\10 -0.5 /' "$model2d" > "$checkroot/emptycellvalue.txt"
        run_formatcheck_case emptycellvalue "$checkroot/emptycellvalue.txt" \
            "has the value -0.5 in the column X_Fegroup, which is outside the limits of that column"
        # a key cells_per_axis that does not agree with the lists of edges
        sed 's/^cells_per_axis 5 10$/cells_per_axis 5 9/' "$model2d" > "$checkroot/wrongcellsperaxis.txt"
        run_formatcheck_case wrongcellsperaxis "$checkroot/wrongcellsperaxis.txt" \
            "the key cells_per_axis gives the cell counts '5 9', but the cell edges give '5 10'"
        run_formatcheck_case yefile "$model2d" "Ye.txt is present, but model.txt has format version 2" file:Ye.txt
        run_formatcheck_case abundancefile "$model2d" \
            "abundances.txt is present, but model.txt has format version 2" file:abundances.txt
        # a mirrored model has the masses of its full twin
        get_masses_of_case() {
            sed -E 's/^[0-9TZ:-]+ //' "$checkroot/$1/inputcheck/output_0-0.txt" |
                grep -E '^(Total input model mass|Total grid-mapped mass|  [0-9A-Z][A-Za-z0-9-]*: [0-9])'
        }
        require_same_masses() {
            if [ -z "$(get_masses_of_case "$1")" ] || [ "$(get_masses_of_case "$1")" != "$(get_masses_of_case "$2")" ]; then
                echo "[error] formatcheck: the cases $1 and $2 have different masses." >&2
                exit 1
            fi
            echo "formatcheck $1 and $2: the same masses"
        }
        # The mirrored model holds the cells with z >= 0 of the 2D model, which has 5 x 10 cells. The key
        # cell_line_count is optional, thus the models from awk have no such key
        modelmirror="$checkroot/cylindrical2d_mirror.txt"
        awk '/^cells_per_axis/ { print "cells_per_axis 5 5"; print "mirror_plane z=0"; next }
            /^vel_z_edges_kmps/ { printf "vel_z_edges_kmps"; for (i = NF - 5; i <= NF; i++) printf " %s", $i; print ""; next }
            /^cell_line_count/ || (/^[0-9]/ && $1 < 25) { next }
            /^[0-9]/ { $1 -= 25 }
            { print }' "$model2d" > "$modelmirror"
        # the twin has each cell of the mirrored model also at its mirror position, in the sequence of the cell ids
        modeltwin="$checkroot/cylindrical2d_mirrortwin.txt"
        {
            awk '/^cells_per_axis/ { print "cells_per_axis 5 10"; next }
                /^mirror_plane/ || /^[0-9]/ { next }
                /^vel_z_edges_kmps/ { printf "vel_z_edges_kmps"; for (i = NF; i > 2; i--) printf " -%s", $i; for (i = 2; i <= NF; i++) printf " %s", $i; print ""; next }
                { print }' "$modelmirror"
            awk '/^[0-9]/ { id = $1; $1 = (id % 5) + 5 * (4 - int(id / 5)); print; $1 = (id % 5) + 5 * (5 + int(id / 5)); print }' \
                "$modelmirror" | sort -n
        } > "$modeltwin"
        run_formatcheck_case cylindrical2d_mirror "$modelmirror" complete "model.txt has mirror symmetry in the plane z = 0"
        run_formatcheck_case cylindrical2d_mirrortwin "$modeltwin" complete
        require_same_masses cylindrical2d_mirror cylindrical2d_mirrortwin
        # Each of the 5 shells of the 1D model in two polar cells from the plane z = 0 to the +z axis. The polar cells
        # have different widths in cos(theta), and the second one has two times the density. Thus a mirror cell at
        # an incorrect polar index gives a different mass
        modelspherical2dmirror="$checkroot/spherical2d_mirror.txt"
        {
            awk '/^coordinates/ { print "coordinates spherical2d"; next }
                /^cells_per_axis/ { print "cells_per_axis 5 2"; print "mirror_plane z=0"; next }
                /^vel_r_edges_kmps/ { print; print "costheta_edges 0.0 0.3 1.0"; next }
                /^cell_line_count/ { next }
                { print }' "$model1d"
            awk '/^[0-9]/ { $1 += 5; $2 *= 2; print }' "$model1d"
        } > "$modelspherical2dmirror"
        # the twin has the polar cells from the -z axis to the +z axis
        modelspherical2dtwin="$checkroot/spherical2d_mirrortwin.txt"
        {
            awk '/^cells_per_axis/ { print "cells_per_axis 5 4"; next }
                /^mirror_plane/ || /^[0-9]/ { next }
                /^costheta_edges/ { print "costheta_edges -1.0 -0.3 0.0 0.3 1.0"; next }
                { print }' "$modelspherical2dmirror"
            awk '/^[0-9]/ { id = $1; $1 = (id % 5) + 5 * (1 - int(id / 5)); print; $1 = (id % 5) + 5 * (2 + int(id / 5)); print }' \
                "$modelspherical2dmirror" | sort -n
        } > "$modelspherical2dtwin"
        run_formatcheck_case spherical2d_mirror "$modelspherical2dmirror" complete
        run_formatcheck_case spherical2d_mirrortwin "$modelspherical2dtwin" complete
        require_same_masses spherical2d_mirror spherical2d_mirrortwin
        sed 's/^mirror_plane z=0$/mirror_plane z=1/' "$modelmirror" > "$checkroot/mirrorvalue.txt"
        run_formatcheck_case mirrorvalue "$checkroot/mirrorvalue.txt" "the key mirror_plane must have the value z=0"
        # BSD sed gives a \n in a replacement as the letter n, thus awk inserts the line
        awk '/^cell_line_count/ { print "mirror_plane z=0" } { print }' "$model1d" > "$checkroot/mirror1d.txt"
        run_formatcheck_case mirror1d "$checkroot/mirror1d.txt" "a spherical1d model takes no key mirror_plane"
        # a polar cell at the z axis that is too thin for the position tolerance of a cell boundary
        sed 's/^costheta_edges 0.0 0.3 1.0$/costheta_edges 0.0 0.999999995 1.0/' "$modelspherical2dmirror" \
            > "$checkroot/thinpolarcell.txt"
        run_formatcheck_case thinpolarcell "$checkroot/thinpolarcell.txt" "are the edges of a cell that is too thin"
        # the polar axis of a mirrored spherical model starts at the plane z = 0
        sed 's/^costheta_edges 0.0 0.3 1.0$/costheta_edges -1.0 0.3 1.0/' "$modelspherical2dmirror" \
            > "$checkroot/mirrorfullpolaraxis.txt"
        run_formatcheck_case mirrorfullpolaraxis "$checkroot/mirrorfullpolaraxis.txt" \
            "the first value of costheta_edges is -1, but it must be 0"
        # the z axis of a mirrored cylindrical model starts at the plane z = 0
        awk '/^cell_line_count/ { print "mirror_plane z=0" } { print }' "$model2d" > "$checkroot/mirrorfullaxis.txt"
        run_formatcheck_case mirrorfullaxis "$checkroot/mirrorfullaxis.txt" \
            "the values of vel_z_edges_kmps must have equal steps from 0 to"
        # two radial edges of the 1D model that differ in km/s and give the same velocity in cm/s
        sed 's/^vel_r_edges_kmps 0.0 /vel_r_edges_kmps 0.0 1000.0000000000002 1000.0000000000003 /' "$model1d" \
            > "$checkroot/equaledges.txt"
        run_formatcheck_case equaledges "$checkroot/equaledges.txt" "The values must also increase in cm/s"
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
        # the check of a resumed simulation needs the rank count of the simulation
        run_inputcheck 4
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
            # diff gives the status 1 for a difference and 2 for an error
            diff_status=$?
            if [ "$diff_status" -eq 1 ]; then
                echo "[error] The checksums differ from ${reffile#"$repopath"/}. In the diff, '<' marks the reference and '>' marks this run." >&2
            fi
            exit "$diff_status"
        fi
        ;;
esac
