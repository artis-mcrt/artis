#!/bin/bash
sbatch -J $(basename $(exec pwd)) --ntasks=960 --ntasks-per-core=1 --constraint=amd --partition=long --time=24:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de --no-requeue -- artis/scripts/artis-virgo-slurmjob.sh
