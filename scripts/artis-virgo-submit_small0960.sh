#!/bin/bash
# 960 cores and 1.5GB per core for simple 1D models (maybe 3D Type Ia without full NLTE?)
sbatch -J $(basename $(exec pwd)) --ntasks=960 --mem-per-cpu=1536M --partition=long --time=24:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh
