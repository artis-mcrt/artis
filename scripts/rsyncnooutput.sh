#!/usr/bin/env bash
rsync -av --exclude='speclc_angle_res' --exclude='exspec*.txt' --exclude="*.pdf" --exclude="out*.txt" --exclude='*.slurm' --exclude='job_from_ts*' --exclude='*.out*' --exclude='packets*' --exclude='logfiles*' --exclude='*.tmp' "$@"
