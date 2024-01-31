#!/bin/bash -x
# 1920 cores and 2GB per core for 3D LTE kilonova models (and maybe 3D NLTE Type Ia?)
sbatch -J $(basename $(exec pwd)) --ntasks=1920 --mem-per-cpu=2000MB --partition=long --time=24:00:00 --constraint=amd,epyc,7713 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh

# 960 cores and 1.5GB per core for simple 1D models (maybe 3D Type Ia without full NLTE?)
#sbatch -J $(basename $(exec pwd)) --ntasks=960 --mem-per-cpu=1536M --partition=long --time=24:00:00 --constraint=amd --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh
