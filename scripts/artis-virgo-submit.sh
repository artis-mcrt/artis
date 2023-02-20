#!/bin/bash
#--ntasks=960 --mem-per-cpu=1.5GB is probably fine for 1D models
#--ntasks=1920 --mem-per-cpu=4GB for 3D kilonova models
sbatch -J $(basename $(exec pwd)) --ntasks=960 --mem-per-cpu=1.5G --partition=long --time=24:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh
