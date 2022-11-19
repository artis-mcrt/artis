#!/bin/bash
sbatch -J ex_$(basename $(exec pwd)) --ntasks=8 --ntasks-per-core=1 --partition=long --constraint=amd --time=48:00:00 --mail-type=ALL --mail-user=${USER}@gsi.de -- artis/scripts/exspec-gzip-virgo-slurmjob.sh
