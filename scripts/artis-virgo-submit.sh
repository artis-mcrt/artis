#!/bin/bash -x
sbatch -J $(basename $(exec pwd)) --ntasks=1920 --ntasks-per-node=128 --mem-per-cpu=2000MB --partition=long --time=24:00:00 --constraint=amd,epyc,7713 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh
