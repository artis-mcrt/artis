#!/usr/bin/env bash
rsync -av --exclude='speclc_angle_res' --exclude='exspec*.txt' --exclude="*.pdf" --exclude="out*.txt" --exclude='*.slurm' --exclude='*.out*' --exclude='packets*' --exclude='logfiles*' --exclude='*.tmp' "$@"
