#!/bin/bash
# 1920 cores and 4GB per core for 3D LTE kilonova models (and maybe 3D NLTE Type Ia?)
sbatch -J $(basename $(exec pwd)) --ntasks=1920 --mem-per-cpu=4096MB --partition=long --time=24:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh
